#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Module 7: genome-scale metabolic model (GEM) generation with gapseq 1.4.0
 * and MEMOTE 0.17.0.
 *
 * CLI contracts are verified against the pinned tool help:
 *   gapseq doall <genome> [medium] [Bacteria|Archaea]
 *   gapseq adapt -m <model> -w <growth> -c <weights> -g <genes> -b <blast>
 *   memote report snapshot --filename <html> <model>
 *   memote run --ignore-git --filename <json> <model>
 */

// ============================================================================
// PARAMETERS
// ============================================================================

params.samwise_dir = java.nio.file.Paths
    .get((params.samwise_dir ?: params.working_dir ?: projectDir).toString())
    .toAbsolutePath()
    .normalize()
    .toString()
params.working_dir = java.nio.file.Paths
    .get((params.working_dir ?: params.samwise_dir).toString())
    .toAbsolutePath()
    .normalize()
    .toString()
params.input_manifest = null
params.protein_fasta_dir = null
params.media_csv = null
params.template_organism = 'Bacteria'
params.run_gapseq_adapt = false
params.adaptation_compounds = null // Retired: use --adapt_manifest.
params.adapt_manifest = null
params.run_memote = true
params.memote_mode = 'snapshot' // snapshot | run | both
params.threads = null
params.tool_env_dir = null
params.gapseq_env_dir = null
params.auto_install = true
params.conda_pkgs_dir = null
params.gapseq_extra_args = '' // Retired: the pinned doall interface is positional.
params.memote_extra_args = ''
params.publish_gems_mode = 'copy'
params.publish_reports_mode = 'copy'
params.results_dir = params.working_dir
params.outdir = "${params.results_dir}/module_7_gems"
params.module6_output_dir = "${params.results_dir}/module_6_magannotate"

// Pinned package specifications. These are intentionally not configurable.
params.gapseq_package = 'gapseq=1.4.0'
params.memote_package = 'memote=0.17.0'

params.media_minimal = "${params.samwise_dir}/background/media/gapseq_M9_glucose_aerobic.csv"
params.media_comprehensive = "${params.samwise_dir}/background/media/gapseq_all_nutrients.csv"

// ============================================================================
// HELPERS AND VALIDATION
// ============================================================================

def absPath(value) {
    def text = value == null ? '' : value.toString().trim()
    if (!text || text in ['null', 'NA']) {
        return ''
    }
    return java.nio.file.Paths.get(text).toAbsolutePath().normalize().toString()
}

def requireRegularFile(String label, String pathText) {
    def candidate = new File(pathText)
    if (!candidate.isFile() || candidate.length() == 0L) {
        error("${label} is missing, not a regular file, or empty: ${pathText}")
    }
}

def parseMagIdsFromManifest(String pathText) {
    def ids = [] as Set
    def rows = new File(pathText).readLines()
    def headerIndex = rows.findIndexOf { row -> row?.trim() }
    if (headerIndex < 0) {
        error("Module 6 manifest is empty: ${pathText}")
    }
    def headers = rows[headerIndex].split(/\t/, -1)*.trim()
    def magIndex = headers.indexOf('mag_id')
    if (magIndex < 0) {
        error("Module 6 manifest lacks required 'mag_id' column: ${pathText}")
    }
    rows.drop(headerIndex + 1).each { line ->
        if (!line?.trim()) return
        def columns = line.split(/\t/, -1)
        if (magIndex < columns.length && columns[magIndex].trim()) {
            ids << columns[magIndex].trim()
        }
    }
    if (!ids) {
        error("Module 6 manifest contains no non-empty mag_id values: ${pathText}")
    }
    return ids
}

def parseAdaptManifestTSV(String pathText, Set upstreamMagIds) {
    def values = [:]
    def rows = new File(pathText).readLines()
    def headerIndex = rows.findIndexOf { row -> row?.trim() }
    if (headerIndex < 0) {
        error("Adapt manifest is empty: ${pathText}")
    }

    def headers = rows[headerIndex].split(/\t/, -1)*.trim()
    def magIndex = headers.indexOf('mag_id')
    def compoundsIndex = headers.indexOf('adapt_compounds')
    if (magIndex < 0 || compoundsIndex < 0) {
        error("Adapt manifest requires tab-separated headers: mag_id and adapt_compounds")
    }

    rows.drop(headerIndex + 1).eachWithIndex { line, offset ->
        if (!line?.trim()) return
        def columns = line.split(/\t/, -1)
        def lineNumber = headerIndex + offset + 2
        if (magIndex >= columns.length || compoundsIndex >= columns.length) {
            error("Adapt manifest line ${lineNumber} is missing a required field")
        }
        def magId = columns[magIndex].trim()
        def compounds = columns[compoundsIndex].trim()
        if (!magId || !compounds) {
            error("Adapt manifest line ${lineNumber} requires non-empty mag_id and adapt_compounds")
        }
        if (values.containsKey(magId)) {
            error("Adapt manifest contains duplicate mag_id: ${magId}")
        }
        if (!upstreamMagIds.contains(magId)) {
            error("Adapt manifest mag_id is absent from the Module 6 manifest: ${magId}")
        }
        compounds.split(/,/).each { token ->
            if (!(token.trim() ==~ /^cpd\d{5}:((TRUE)|(FALSE))$/)) {
                error("Invalid adapt compound '${token}' for ${magId}; expected cpd#####:(TRUE|FALSE)")
            }
        }
        values[magId] = compounds
    }

    if (!values) {
        error("Adapt manifest contains no data rows: ${pathText}")
    }
    return values
}

// ============================================================================
// PREPARE MODULE 6 INPUTS
// ============================================================================

process PREPARE_GAPSEQ_INPUTS {
    tag 'prepare_gapseq_inputs'

    publishDir "${params.outdir}/inputs", mode: 'copy', pattern: 'gapseq_input_manifest.tsv'
    publishDir "${params.outdir}/inputs", mode: 'copy', pattern: 'gapseq_inputs_stats.tsv'
    publishDir "${params.outdir}/inputs", mode: "${params.publish_gems_mode}", pattern: 'protein_fastas'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'prepare_gapseq_inputs.log'

    input:
    val eggnog_genepred_fasta
    val eggnog_input_manifest
    path setup_status

    output:
    path 'gapseq_input_manifest.tsv', emit: manifest
    path 'gapseq_inputs_stats.tsv', emit: stats
    path 'protein_fastas', emit: protein_fastas_dir
    path 'prepare_gapseq_inputs.log', emit: log

    script:
    """
    set -euo pipefail
    mkdir -p protein_fastas
    LOG=prepare_gapseq_inputs.log
    echo "PREPARE_GAPSEQ_INPUTS started: \$(date)" > "\$LOG"

    python3 - "${eggnog_genepred_fasta}" "${eggnog_input_manifest}" "\$LOG" <<'PY'
import csv
import sys
from pathlib import Path

fasta_path = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
log_path = Path(sys.argv[3])
output_dir = Path("protein_fastas")


def log(message):
    with log_path.open("a") as handle:
        print(message, file=handle)

if not fasta_path.is_file() or fasta_path.stat().st_size == 0:
    raise SystemExit(f"ERROR: Unified eggNOG predicted-protein FASTA missing or empty: {fasta_path}")
if not manifest_path.is_file() or manifest_path.stat().st_size == 0:
    raise SystemExit(f"ERROR: Module 6 manifest missing or empty: {manifest_path}")

with manifest_path.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if not reader.fieldnames or "mag_id" not in reader.fieldnames:
        raise SystemExit("ERROR: Module 6 manifest requires a mag_id column")
    mag_ids = {row["mag_id"].strip() for row in reader if row.get("mag_id", "").strip()}

if not mag_ids:
    raise SystemExit("ERROR: Module 6 manifest contains no MAG IDs")

counts = {mag_id: 0 for mag_id in mag_ids}
handles = {}
unmapped = 0
records = 0


def write_record(header, sequence):
    global unmapped, records
    if not header or not sequence:
        return
    mag_id = header.split("|", 1)[0].split()[0]
    if mag_id not in counts:
        unmapped += 1
        return
    handle = handles.get(mag_id)
    if handle is None:
        handle = (output_dir / f"{mag_id}.faa").open("w")
        handles[mag_id] = handle
    handle.write(f">{header}\\n{sequence}\\n")
    counts[mag_id] += 1
    records += 1

try:
    header = None
    sequence_parts = []
    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                write_record(header, "".join(sequence_parts))
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line)
    write_record(header, "".join(sequence_parts))
finally:
    for handle in handles.values():
        handle.close()

with Path("gapseq_input_manifest.tsv").open("w") as handle:
    print("mag_id\tprotein_fasta", file=handle)
    for mag_id in sorted(mag_ids):
        if counts[mag_id]:
            print(f"{mag_id}\t{(output_dir / f'{mag_id}.faa').resolve()}", file=handle)
        else:
            log(f"WARNING: No predicted proteins mapped to MAG: {mag_id}")

with Path("gapseq_inputs_stats.tsv").open("w") as handle:
    print("metric\tvalue", file=handle)
    print(f"total_mags\t{len(mag_ids)}", file=handle)
    print(f"mags_with_sequences\t{sum(value > 0 for value in counts.values())}", file=handle)
    print(f"total_sequences_mapped\t{records}", file=handle)
    print(f"sequences_unmapped\t{unmapped}", file=handle)

if not records:
    raise SystemExit("ERROR: No predicted proteins matched Module 6 MAG IDs")
log(f"Prepared {sum(value > 0 for value in counts.values())} per-MAG .faa protein FASTAs from {records} proteins")
PY

    echo "PREPARE_GAPSEQ_INPUTS completed: \$(date)" >> "\$LOG"
    """

    stub:
    """
    mkdir -p protein_fastas
    printf 'mag_id\\tprotein_fasta\\nMAG_STUB\\t%s/protein_fastas/MAG_STUB.faa\\n' "\$PWD" > gapseq_input_manifest.tsv
    printf '>MAG_STUB_protein\\nMSTUBSEQ\\n' > protein_fastas/MAG_STUB.faa
    printf 'metric\\tvalue\\ntotal_mags\\t1\\nmags_with_sequences\\t1\\ntotal_sequences_mapped\\t1\\nsequences_unmapped\\t0\\n' > gapseq_inputs_stats.tsv
    echo 'PREPARE_GAPSEQ_INPUTS stub' > prepare_gapseq_inputs.log
    """
}

// ============================================================================
// PINNED TOOL ENVIRONMENT
// ============================================================================

process SETUP_GAPSEQ {
    tag 'setup_gapseq_memote'

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: 'gapseq_setup_status.env'

    output:
    path 'gapseq_setup_status.env', emit: status

    script:
    def managedEnv = !params.gapseq_env_dir && !params.tool_env_dir
    def envDir = params.gapseq_env_dir ? absPath(params.gapseq_env_dir) : (params.tool_env_dir ? "${absPath(params.tool_env_dir)}/gapseq" : "${absPath(params.outdir)}/conda_envs/gapseq")
    def packageCache = params.conda_pkgs_dir ? absPath(params.conda_pkgs_dir) : "${absPath(params.outdir)}/conda_pkgs/gapseq"

    """
    set -euo pipefail
    export PYTHONNOUSERSITE=1
    unset PYTHONPATH
    STATUS=gapseq_setup_status.env
    GAPSEQ_ENV="${envDir}"
    CONDA_PKGS_DIRS="${packageCache}"
    MANAGED_ENV="${managedEnv}"
    export CONDA_PKGS_DIRS
    mkdir -p "\$CONDA_PKGS_DIRS"

    echo "GAPSEQ_ENV=\$GAPSEQ_ENV" > "\$STATUS"
    echo "GAPSEQ_PACKAGE=${params.gapseq_package}" >> "\$STATUS"
    echo "MEMOTE_PACKAGE=${params.memote_package}" >> "\$STATUS"
    echo "MANAGED_BY_SAMWISE=\$MANAGED_ENV" >> "\$STATUS"
    echo "CONDA_PKGS_DIRS=\$CONDA_PKGS_DIRS" >> "\$STATUS"

    find_installer() {
        if command -v mamba >/dev/null 2>&1; then echo mamba
        elif command -v conda >/dev/null 2>&1; then echo conda
        else echo ''
        fi
    }

    validate_env() {
        local prefix="\$1"
        [[ -x "\$prefix/bin/gapseq" ]] || return 1
        [[ -x "\$prefix/bin/memote" ]] || return 1
        export PATH="\$prefix/bin:\$PATH"
        local gapseq_version
        local gapseq_status
        local memote_version
        local memote_status
        if gapseq_version="\$("\$prefix/bin/gapseq" -v 2>&1)"; then
            :
        else
            gapseq_status=\$?
            {
                echo "ERROR: gapseq version check failed with exit status \$gapseq_status."
                printf '%s\\n' "\$gapseq_version"
            } >> "\$STATUS"
            printf '%s\\n' "\$gapseq_version" >&2
            return "\$gapseq_status"
        fi
        if memote_version="\$("\$prefix/bin/memote" --version 2>&1)"; then
            :
        else
            memote_status=\$?
            {
                echo "ERROR: MEMOTE version check failed with exit status \$memote_status."
                printf '%s\\n' "\$memote_version"
            } >> "\$STATUS"
            printf '%s\\n' "\$memote_version" >&2
            return "\$memote_status"
        fi
        printf '%s\\n' "\$gapseq_version" | grep -Eq '(^|[^0-9])1\\.4\\.0([^0-9]|\$)' || return 1
        printf '%s\\n' "\$memote_version" | grep -Eq '(^|[^0-9])0\\.17\\.0([^0-9]|\$)' || return 1
    }

    if [[ -d "\$GAPSEQ_ENV" ]] && ! validate_env "\$GAPSEQ_ENV"; then
        if [[ "\$MANAGED_ENV" == true ]]; then
            echo 'Removing stale SAMWISE-managed pinned environment.' >> "\$STATUS"
            rm -rf "\$GAPSEQ_ENV"
        else
            echo 'ERROR: User-managed environment is missing pinned gapseq 1.4.0 and/or MEMOTE 0.17.0; it will not be modified.' >> "\$STATUS"
            exit 1
        fi
    fi

    if [[ ! -d "\$GAPSEQ_ENV" ]]; then
        if [[ "\$MANAGED_ENV" != true ]]; then
            echo 'ERROR: User-managed environment does not exist; SAMWISE will not create it.' >> "\$STATUS"
            exit 1
        fi
        if [[ "${params.auto_install}" != true ]]; then
            echo 'ERROR: Pinned environment is absent and --auto_install is false.' >> "\$STATUS"
            exit 1
        fi
        INSTALLER="\$(find_installer)"
        if [[ -z "\$INSTALLER" ]]; then
            echo 'ERROR: Neither mamba nor conda is available.' >> "\$STATUS"
            exit 1
        fi
        mkdir -p "\$(dirname "\$GAPSEQ_ENV")"
        "\$INSTALLER" create -y -p "\$GAPSEQ_ENV" --override-channels -c conda-forge -c bioconda \\
            '${params.gapseq_package}' '${params.memote_package}' 'python>=3.8' 'libsbml' >> "\$STATUS" 2>&1
    fi

    if ! validate_env "\$GAPSEQ_ENV"; then
        echo 'ERROR: Environment validation failed; required versions are gapseq 1.4.0 and MEMOTE 0.17.0.' >> "\$STATUS"
        exit 1
    fi

    export PATH="\$GAPSEQ_ENV/bin:\$PATH"
    echo "GAPSEQ_VERSION=\$(gapseq -v 2>&1 | tr '\\n' ' ')" >> "\$STATUS"
    echo "MEMOTE_VERSION=\$(memote --version 2>&1 | tr '\\n' ' ')" >> "\$STATUS"
    echo "SETUP_COMPLETED=\$(date -u +%FT%TZ)" >> "\$STATUS"
    """

    stub:
    """
    printf 'GAPSEQ_ENV=%s\\nGAPSEQ_VERSION=1.4.0\\nMEMOTE_VERSION=0.17.0\\nMANAGED_BY_SAMWISE=true\\n' "\$PWD" > gapseq_setup_status.env
    """
}

// ============================================================================
// GAPSEQ MODEL CONSTRUCTION
// ============================================================================

process RUN_GAPSEQ_DOALL {
    tag "gapseq_doall_${mag_id}"

    publishDir "${params.outdir}/gems/doall", mode: "${params.publish_gems_mode}"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    cpus { params.threads ? (params.threads as int) : 1 }

    input:
    tuple val(mag_id), path(protein_fasta)
    val media_csv
    val template_organism
    path setup_status

    output:
    tuple val(mag_id), path("${mag_id}_gapseq_doall"), val('doall'), emit: gems
    path "${mag_id}_gapseq_doall.log", emit: log
    path "${mag_id}_doall_metadata.tsv", emit: metadata

    script:
    def outputDir = "${mag_id}_gapseq_doall"
    def logFile = "${mag_id}_gapseq_doall.log"

    """
    set -euo pipefail
    LOG="${logFile}"
    ENV_DIR="\$(grep '^GAPSEQ_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    export PATH="\$ENV_DIR/bin:\$PATH"
    echo "RUN_GAPSEQ_DOALL started for ${mag_id}: \$(date)" > "\$LOG"
    echo "gapseq doall ${protein_fasta} ${media_csv} ${template_organism}" >> "\$LOG"

    mkdir -p "${outputDir}"
    cd "${outputDir}"
    gapseq doall "../${protein_fasta}" "${media_csv}" "${template_organism}" >> "../\$LOG" 2>&1
    cd ..

    require_one() {
        local label="\$1"
        shift
        local -a matches=("\$@")
        if [[ "\${#matches[@]}" -ne 1 || ! -f "\${matches[0]}" ]]; then
            echo "ERROR: Expected exactly one \$label; found \${#matches[@]}" >> "\$LOG"
            printf '%s\\n' "\${matches[@]}" >> "\$LOG"
            exit 1
        fi
        printf '%s' "\${matches[0]}"
    }

    MODEL_RDS_MATCHES=()
    while IFS= read -r match; do MODEL_RDS_MATCHES+=("\$match"); done < <(find "${outputDir}" -type f -name '*.RDS' ! -name '*-draft.RDS' ! -name '*-rxnWeights.RDS' ! -name '*-rxnXgenes.RDS' | sort)
    MODEL_XML_MATCHES=()
    while IFS= read -r match; do MODEL_XML_MATCHES+=("\$match"); done < <(find "${outputDir}" -type f -name '*.xml' ! -name '*-draft.xml' | sort)
    WEIGHTS_MATCHES=()
    while IFS= read -r match; do WEIGHTS_MATCHES+=("\$match"); done < <(find "${outputDir}" -type f -name '*-rxnWeights.RDS' | sort)
    GENES_MATCHES=()
    while IFS= read -r match; do GENES_MATCHES+=("\$match"); done < <(find "${outputDir}" -type f -name '*-rxnXgenes.RDS' | sort)
    REACTIONS_MATCHES=()
    while IFS= read -r match; do REACTIONS_MATCHES+=("\$match"); done < <(find "${outputDir}" -type f -name '*-all-Reactions.tbl' | sort)

    MODEL_RDS="\$(require_one 'model RDS' "\${MODEL_RDS_MATCHES[@]}")"
    MODEL_XML="\$(require_one 'model SBML XML' "\${MODEL_XML_MATCHES[@]}")"
    WEIGHTS="\$(require_one 'rxnWeights RDS' "\${WEIGHTS_MATCHES[@]}")"
    GENES="\$(require_one 'rxnXgenes RDS' "\${GENES_MATCHES[@]}")"
    REACTIONS="\$(require_one 'all-Reactions table' "\${REACTIONS_MATCHES[@]}")"

    published_dir="${params.outdir}/gems/doall/${outputDir}"
    rel_path() { printf '%s' "\${1#${outputDir}/}"; }
    {
        printf 'mag_id\\tstage\\tmodel_rds\\tmodel_xml\\trxn_weights\\trxn_genes\\treactions_tbl\\tlog\\n'
        printf '%s\\tdoall\\t%s/%s\\t%s/%s\\t%s/%s\\t%s/%s\\t%s/%s\\t%s/%s\\n' \\
            '${mag_id}' "\$published_dir" "\$(rel_path "\$MODEL_RDS")" "\$published_dir" "\$(rel_path "\$MODEL_XML")" \\
            "\$published_dir" "\$(rel_path "\$WEIGHTS")" "\$published_dir" "\$(rel_path "\$GENES")" \\
            "\$published_dir" "\$(rel_path "\$REACTIONS")" '${params.outdir}/logs' "\$LOG"
    } > "${mag_id}_doall_metadata.tsv"

    # Ensure files are fully written and visible on NFS
    sync
    sleep 2
    
    # Verify critical files exist and are non-empty
    for file in "\$MODEL_RDS" "\$MODEL_XML" "\$WEIGHTS" "\$GENES" "\$REACTIONS"; do
        if [[ ! -f "\$file" ]] || [[ ! -s "\$file" ]]; then
            echo "ERROR: Required output file missing or empty: \$file" >> "\$LOG"
            ls -lh "\$(dirname "\$file")" >> "\$LOG" 2>&1
            exit 1
        fi
    done
    echo "Validated all output files exist and are non-empty" >> "\$LOG"

    echo "RUN_GAPSEQ_DOALL completed for ${mag_id}: \$(date)" >> "\$LOG"
    """

    stub:
    def outputDir = "${mag_id}_gapseq_doall"
    def logFile = "${mag_id}_gapseq_doall.log"
    """
    mkdir -p "${outputDir}"
    printf 'stub' > "${outputDir}/${mag_id}.RDS"
    printf '<sbml/>' > "${outputDir}/${mag_id}.xml"
    printf 'stub' > "${outputDir}/${mag_id}-rxnWeights.RDS"
    printf 'stub' > "${outputDir}/${mag_id}-rxnXgenes.RDS"
    printf 'stub' > "${outputDir}/${mag_id}-all-Reactions.tbl"
    echo 'RUN_GAPSEQ_DOALL stub' > "${logFile}"
    printf 'mag_id\\tstage\\tmodel_rds\\tmodel_xml\\trxn_weights\\trxn_genes\\treactions_tbl\\tlog\\n%s\\tdoall\\t%s/gems/doall/${outputDir}/${mag_id}.RDS\\t%s/gems/doall/${outputDir}/${mag_id}.xml\\t%s/gems/doall/${outputDir}/${mag_id}-rxnWeights.RDS\\t%s/gems/doall/${outputDir}/${mag_id}-rxnXgenes.RDS\\t%s/gems/doall/${outputDir}/${mag_id}-all-Reactions.tbl\\t%s/logs/${logFile}\\n' '${mag_id}' '${params.outdir}' '${params.outdir}' '${params.outdir}' '${params.outdir}' '${params.outdir}' '${params.outdir}' > "${mag_id}_doall_metadata.tsv"
    """
}

// ============================================================================
// MAG-SPECIFIC ADAPTATION
// ============================================================================

process RUN_GAPSEQ_ADAPT {
    tag "gapseq_adapt_${mag_id}"

    publishDir "${params.outdir}/gems/adapt", mode: "${params.publish_gems_mode}"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(mag_id), path(doall_dir), val(adapt_compounds)
    path setup_status

    output:
    tuple val(mag_id), path("${mag_id}_gapseq_adapt"), val('adapt'), emit: gems
    path "${mag_id}_gapseq_adapt.log", emit: log
    path "${mag_id}_adapt_metadata.tsv", emit: metadata

    script:
    def outputDir = "${mag_id}_gapseq_adapt"
    def logFile = "${mag_id}_gapseq_adapt.log"

    """
    set -euo pipefail
    LOG="${logFile}"
    ENV_DIR="\$(grep '^GAPSEQ_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    export PATH="\$ENV_DIR/bin:\$PATH"
    echo "RUN_GAPSEQ_ADAPT started for ${mag_id}: \$(date)" > "\$LOG"

    require_one() {
        local label="\$1"
        shift
        local -a matches=("\$@")
        if [[ "\${#matches[@]}" -ne 1 || ! -f "\${matches[0]}" ]]; then
            echo "ERROR: Expected exactly one \$label; found \${#matches[@]}" >> "\$LOG"
            printf '%s\\n' "\${matches[@]}" >> "\$LOG"
            exit 1
        fi
        printf '%s' "\${matches[0]}"
    }

    MODEL_MATCHES=()
    while IFS= read -r match; do MODEL_MATCHES+=("\$match"); done < <(find -L "${doall_dir}" -type f -name '*.RDS' ! -name '*-draft.RDS' ! -name '*-rxnWeights.RDS' ! -name '*-rxnXgenes.RDS' | sort)
    WEIGHTS_MATCHES=()
    while IFS= read -r match; do WEIGHTS_MATCHES+=("\$match"); done < <(find -L "${doall_dir}" -type f -name '*-rxnWeights.RDS' | sort)
    GENES_MATCHES=()
    while IFS= read -r match; do GENES_MATCHES+=("\$match"); done < <(find -L "${doall_dir}" -type f -name '*-rxnXgenes.RDS' | sort)
    REACTIONS_MATCHES=()
    while IFS= read -r match; do REACTIONS_MATCHES+=("\$match"); done < <(find -L "${doall_dir}" -type f -name '*-all-Reactions.tbl' | sort)
    MODEL_RDS="\$(require_one 'model RDS' "\${MODEL_MATCHES[@]}")"
    WEIGHTS="\$(require_one 'rxnWeights RDS' "\${WEIGHTS_MATCHES[@]}")"
    GENES="\$(require_one 'rxnXgenes RDS' "\${GENES_MATCHES[@]}")"
    REACTIONS="\$(require_one 'all-Reactions table' "\${REACTIONS_MATCHES[@]}")"

    # Also find the original XML for fallback
    MODEL_XML_MATCHES=()
    while IFS= read -r match; do MODEL_XML_MATCHES+=("\$match"); done < <(find -L "${doall_dir}" -type f -name '*.xml' ! -name '*-draft.xml' | sort)
    MODEL_XML="\$(require_one 'model XML' "\${MODEL_XML_MATCHES[@]}")"

    mkdir -p "${outputDir}"
    echo "gapseq adapt -m \$MODEL_RDS -w ${adapt_compounds} -c \$WEIGHTS -g \$GENES -b \$REACTIONS -f ${outputDir}" >> "\$LOG"
    gapseq adapt -m "\$MODEL_RDS" -w "${adapt_compounds}" -c "\$WEIGHTS" -g "\$GENES" -b "\$REACTIONS" -f "${outputDir}" >> "\$LOG" 2>&1

    # Check if adapt actually produced new files
    ADAPTED_RDS_MATCHES=()
    while IFS= read -r match; do ADAPTED_RDS_MATCHES+=("\$match"); done < <(find "${outputDir}" -type f -name '*-adapt.RDS' | sort)
    ADAPTED_XML_MATCHES=()
    while IFS= read -r match; do ADAPTED_XML_MATCHES+=("\$match"); done < <(find "${outputDir}" -type f -name '*-adapt.xml' | sort)
    
    ADAPTATION_STATUS="success"
    if [[ "\${#ADAPTED_RDS_MATCHES[@]}" -eq 1 && "\${#ADAPTED_XML_MATCHES[@]}" -eq 1 ]]; then
        # Adapt produced new files - use them
        ADAPTED_RDS="\${ADAPTED_RDS_MATCHES[0]}"
        ADAPTED_XML="\${ADAPTED_XML_MATCHES[0]}"
        echo "Adaptation produced new model files" >> "\$LOG"
    elif [[ "\${#ADAPTED_RDS_MATCHES[@]}" -eq 0 && "\${#ADAPTED_XML_MATCHES[@]}" -eq 0 ]]; then
        # No adapt files - model was already growing, copy originals to output dir
        echo "No adaptation files produced (model already growing); using original doall model" >> "\$LOG"
        cp "\$MODEL_RDS" "${outputDir}/"
        cp "\$MODEL_XML" "${outputDir}/"
        ADAPTED_RDS="${outputDir}/\$(basename "\$MODEL_RDS")"
        ADAPTED_XML="${outputDir}/\$(basename "\$MODEL_XML")"
        ADAPTATION_STATUS="no_changes"
    else
        # Unexpected state - partial output
        echo "ERROR: Expected either 0 or 1 adapted files, found RDS=\${#ADAPTED_RDS_MATCHES[@]} XML=\${#ADAPTED_XML_MATCHES[@]}" >> "\$LOG"
        printf '%s\\n' "\${ADAPTED_RDS_MATCHES[@]}" "\${ADAPTED_XML_MATCHES[@]}" >> "\$LOG"
        exit 1
    fi

    published_dir="${params.outdir}/gems/adapt/${outputDir}"
    rel_path() { printf '%s' "\${1#${outputDir}/}"; }
    {
        printf 'mag_id\\tstage\\tmodel_rds\\tmodel_xml\\tadaptation_status\\tlog\\n'
        printf '%s\\tadapt\\t%s/%s\\t%s/%s\\t%s\\t%s/%s\\n' '${mag_id}' \\
            "\$published_dir" "\$(rel_path "\$ADAPTED_RDS")" "\$published_dir" "\$(rel_path "\$ADAPTED_XML")" \\
            "\$ADAPTATION_STATUS" '${params.outdir}/logs' "\$LOG"
    } > "${mag_id}_adapt_metadata.tsv"
    echo "RUN_GAPSEQ_ADAPT completed for ${mag_id}: \$(date)" >> "\$LOG"
    """

    stub:
    def outputDir = "${mag_id}_gapseq_adapt"
    def logFile = "${mag_id}_gapseq_adapt.log"
    """
    mkdir -p "${outputDir}"
    printf 'stub' > "${outputDir}/${mag_id}-adapt.RDS"
    printf '<sbml/>' > "${outputDir}/${mag_id}-adapt.xml"
    echo 'RUN_GAPSEQ_ADAPT stub' > "${logFile}"
    printf 'mag_id\\tstage\\tmodel_rds\\tmodel_xml\\tlog\\n%s\\tadapt\\t%s/gems/adapt/${outputDir}/${mag_id}-adapt.RDS\\t%s/gems/adapt/${outputDir}/${mag_id}-adapt.xml\\t%s/logs/${logFile}\\n' '${mag_id}' '${params.outdir}' '${params.outdir}' '${params.outdir}' > "${mag_id}_adapt_metadata.tsv"
    """
}

// ============================================================================
// MEMOTE VALIDATION
// ============================================================================

process RUN_MEMOTE_SNAPSHOT {
    tag "memote_snapshot_${mag_id}"

    publishDir "${params.outdir}/reports/snapshot", mode: "${params.publish_reports_mode}", pattern: '*.html'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(mag_id), path(gems_dir), val(model_stage)
    path setup_status

    output:
    path "${mag_id}_memote_snapshot.html", emit: report
    path "${mag_id}_memote_snapshot.log", emit: log
    path "${mag_id}_snapshot_metadata.tsv", emit: metadata

    script:
    def reportFile = "${mag_id}_memote_snapshot.html"
    def logFile = "${mag_id}_memote_snapshot.log"

    """
    set -euo pipefail
    LOG="${logFile}"
    ENV_DIR="\$(grep '^GAPSEQ_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    export PATH="\$ENV_DIR/bin:\$PATH"
    export PYTHONNOUSERSITE=1
    unset PYTHONPATH
    XML_MATCHES=()
    if [[ "${model_stage}" == 'doall' ]]; then
        while IFS= read -r match; do XML_MATCHES+=("\$match"); done < <(find -L "${gems_dir}" -type f -name '*.xml' ! -name '*-draft.xml' | sort)
    elif [[ "${model_stage}" == 'adapt' ]]; then
        # Look for ANY final XML file (could be -adapt.xml or copied original .xml)
        while IFS= read -r match; do XML_MATCHES+=("\$match"); done < <(find -L "${gems_dir}" -type f -name '*.xml' ! -name '*-draft.xml' | sort)
    else
        echo "ERROR: Unsupported model stage: ${model_stage}" > "\$LOG"
        exit 1
    fi
    if [[ "\${#XML_MATCHES[@]}" -ne 1 ]]; then
        echo "ERROR: Expected exactly one final XML model; found \${#XML_MATCHES[@]}" > "\$LOG"
        exit 1
    fi
    MODEL_XML="\${XML_MATCHES[0]}"
    if [[ ! -f "\$MODEL_XML" ]] || [[ ! -s "\$MODEL_XML" ]]; then
        echo "ERROR: Model XML file is missing or empty: \$MODEL_XML" >> "\$LOG"
        exit 1
    fi
    echo "memote report snapshot --filename ${reportFile} \$MODEL_XML" > "\$LOG"
    memote report snapshot --filename "${reportFile}" ${params.memote_extra_args} "\$MODEL_XML" >> "\$LOG" 2>&1
    printf 'mag_id\\tstage\\tsnapshot_html\\tlog\\n%s\\t%s\\t%s/reports/snapshot/${reportFile}\\t%s/logs/${logFile}\\n' '${mag_id}' '${model_stage}' '${params.outdir}' '${params.outdir}' > "${mag_id}_snapshot_metadata.tsv"
    """

    stub:
    def reportFile = "${mag_id}_memote_snapshot.html"
    def logFile = "${mag_id}_memote_snapshot.log"
    """
    echo '<html></html>' > "${reportFile}"
    echo 'RUN_MEMOTE_SNAPSHOT stub' > "${logFile}"
    printf 'mag_id\\tstage\\tsnapshot_html\\tlog\\n%s\\t%s\\t%s/reports/snapshot/${reportFile}\\t%s/logs/${logFile}\\n' '${mag_id}' '${model_stage}' '${params.outdir}' '${params.outdir}' > "${mag_id}_snapshot_metadata.tsv"
    """
}

process RUN_MEMOTE_RUN {
    tag "memote_run_${mag_id}"

    publishDir "${params.outdir}/reports/run", mode: "${params.publish_reports_mode}", pattern: '*.json'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.log'

    input:
    tuple val(mag_id), path(gems_dir), val(model_stage)
    path setup_status

    output:
    path "${mag_id}_memote_run.json", emit: report
    path "${mag_id}_memote_run.log", emit: log
    path "${mag_id}_run_metadata.tsv", emit: metadata

    script:
    def reportFile = "${mag_id}_memote_run.json"
    def logFile = "${mag_id}_memote_run.log"

    """
    set -euo pipefail
    LOG="${logFile}"
    ENV_DIR="\$(grep '^GAPSEQ_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    export PATH="\$ENV_DIR/bin:\$PATH"
    export PYTHONNOUSERSITE=1
    unset PYTHONPATH
    XML_MATCHES=()
    if [[ "${model_stage}" == 'doall' ]]; then
        while IFS= read -r match; do XML_MATCHES+=("\$match"); done < <(find -L "${gems_dir}" -type f -name '*.xml' ! -name '*-draft.xml' | sort)
    elif [[ "${model_stage}" == 'adapt' ]]; then
     # Look for ANY final XML file (could be -adapt.xml or copied original .xml)
        while IFS= read -r match; do XML_MATCHES+=("\$match"); done < <(find -L "${gems_dir}" -type f -name '*.xml' ! -name '*-draft.xml' | sort)
    else
        echo "ERROR: Unsupported model stage: ${model_stage}" > "\$LOG"
        exit 1
    fi
    if [[ "\${#XML_MATCHES[@]}" -ne 1 ]]; then
        echo "ERROR: Expected exactly one final XML model; found \${#XML_MATCHES[@]}" > "\$LOG"
        exit 1
    fi
    MODEL_XML="\${XML_MATCHES[0]}"
    if [[ ! -f "\$MODEL_XML" ]] || [[ ! -s "\$MODEL_XML" ]]; then
        echo "ERROR: Model XML file is missing or empty: \$MODEL_XML" >> "\$LOG"
        exit 1
    fi
    echo "memote run --ignore-git --filename ${reportFile} \$MODEL_XML" > "\$LOG"
    memote run --ignore-git --filename "${reportFile}" ${params.memote_extra_args} "\$MODEL_XML" >> "\$LOG" 2>&1
    printf 'mag_id\\tstage\\trun_json\\tlog\\n%s\\t%s\\t%s/reports/run/${reportFile}\\t%s/logs/${logFile}\\n' '${mag_id}' '${model_stage}' '${params.outdir}' '${params.outdir}' > "${mag_id}_run_metadata.tsv"
    """

    stub:
    def reportFile = "${mag_id}_memote_run.json"
    def logFile = "${mag_id}_memote_run.log"
    """
    echo '{}' > "${reportFile}"
    echo 'RUN_MEMOTE_RUN stub' > "${logFile}"
    printf 'mag_id\\tstage\\trun_json\\tlog\\n%s\\t%s\\t%s/reports/run/${reportFile}\\t%s/logs/${logFile}\\n' '${mag_id}' '${model_stage}' '${params.outdir}' '${params.outdir}' > "${mag_id}_run_metadata.tsv"
    """
}

// ============================================================================
// FINAL SAMWISE MANIFEST
// ============================================================================

process WRITE_OUTPUT_MANIFEST {
    tag 'write_module_7_manifest'

    publishDir "${params.outdir}/summary", mode: 'copy'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: 'write_output_manifest.log'

    input:
    path input_manifest
    path doall_metadata
    path adapt_metadata
    path snapshot_metadata
    path run_metadata

    output:
    path 'module_7_gems_manifest.tsv', emit: manifest
    path 'module_7_gems_summary.tsv', emit: summary
    path 'write_output_manifest.log', emit: log

    script:
    """
    set -euo pipefail
    LOG=write_output_manifest.log
    echo "WRITE_OUTPUT_MANIFEST started: \$(date)" > "\$LOG"

    python3 - "${input_manifest}" "${params.outdir}" <<'PY'
import csv
import sys
from pathlib import Path

input_manifest = Path(sys.argv[1])
outdir = sys.argv[2]

def load_rows(pattern):
    rows = {}
    for path in sorted(Path('.').glob(pattern)):
        # Skip empty marker files or directories
        if not path.is_file() or path.stat().st_size == 0:
            continue
        with path.open(newline='') as handle:
            for row in csv.DictReader(handle, delimiter='\\t'):
                mag_id = row.get('mag_id', '').strip()
                if not mag_id:
                    raise SystemExit(f"ERROR: Metadata file lacks mag_id values: {path}")
                if mag_id in rows:
                    raise SystemExit(f"ERROR: Duplicate metadata for MAG {mag_id} in {pattern}")
                rows[mag_id] = row
    return rows

with input_manifest.open(newline='') as handle:
    source_rows = list(csv.DictReader(handle, delimiter='\\t'))

mag_ids = [row['mag_id'].strip() for row in source_rows if row.get('mag_id', '').strip()]
if not mag_ids:
    raise SystemExit('ERROR: Prepared gapseq input manifest contains no MAGs')
if len(set(mag_ids)) != len(mag_ids):
    raise SystemExit('ERROR: Prepared gapseq input manifest contains duplicate MAG IDs')

doall = load_rows('*_doall_metadata.tsv')
adapt = load_rows('*_adapt_metadata.tsv')
snapshot = load_rows('*_snapshot_metadata.tsv')
run = load_rows('*_run_metadata.tsv')

missing_doall = sorted(set(mag_ids) - set(doall))
if missing_doall:
    raise SystemExit(f"ERROR: Missing doall metadata for MAGs: {', '.join(missing_doall)}")

headers = [
    'mag_id', 'protein_fasta',
    'doall_model_rds', 'doall_model_xml', 'doall_rxn_weights', 'doall_rxn_genes',
    'doall_reactions_tbl', 'doall_log',
    'adapted', 'adapt_model_rds', 'adapt_model_xml', 'adapt_log',
    'final_model_stage', 'final_model_rds', 'final_model_xml',
    'memote_snapshot_html', 'memote_snapshot_log',
    'memote_run_json', 'memote_run_log',
]

with Path('module_7_gems_manifest.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=headers, delimiter='\\t', lineterminator='\\n')
    writer.writeheader()
    for mag_id in sorted(mag_ids):
        d = doall[mag_id]
        a = adapt.get(mag_id, {})
        final = a if a else d
        s = snapshot.get(mag_id, {})
        r = run.get(mag_id, {})
        writer.writerow({
            'mag_id': mag_id,
            'protein_fasta': f"{outdir}/inputs/protein_fastas/{mag_id}.faa",
            'doall_model_rds': d.get('model_rds', ''),
            'doall_model_xml': d.get('model_xml', ''),
            'doall_rxn_weights': d.get('rxn_weights', ''),
            'doall_rxn_genes': d.get('rxn_genes', ''),
            'doall_reactions_tbl': d.get('reactions_tbl', ''),
            'doall_log': d.get('log', ''),
            'adapted': 'true' if a else 'false',
            'adapt_model_rds': a.get('model_rds', ''),
            'adapt_model_xml': a.get('model_xml', ''),
            'adapt_log': a.get('log', ''),
            'final_model_stage': final.get('stage', ''),
            'final_model_rds': final.get('model_rds', ''),
            'final_model_xml': final.get('model_xml', ''),
            'memote_snapshot_html': s.get('snapshot_html', ''),
            'memote_snapshot_log': s.get('log', ''),
            'memote_run_json': r.get('run_json', ''),
            'memote_run_log': r.get('log', ''),
        })

with Path('module_7_gems_summary.tsv').open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\\t', lineterminator='\\n')
    writer.writerow(['metric', 'value'])
    writer.writerow(['total_mags', len(mag_ids)])
    writer.writerow(['adapted_mags', len(adapt)])
    writer.writerow(['doall_only_mags', len(mag_ids) - len(adapt)])
    writer.writerow(['memote_snapshot_validated_mags', len(snapshot)])
    writer.writerow(['memote_run_validated_mags', len(run)])
PY

    echo "WRITE_OUTPUT_MANIFEST completed: \$(date)" >> "\$LOG"
    """

    stub:
    """
    printf 'mag_id\\tprotein_fasta\\tdoall_model_rds\\tdoall_model_xml\\tdoall_rxn_weights\\tdoall_rxn_genes\\tdoall_reactions_tbl\\tdoall_log\\tadapted\\tadapt_model_rds\\tadapt_model_xml\\tadapt_log\\tfinal_model_stage\\tfinal_model_rds\\tfinal_model_xml\\tmemote_snapshot_html\\tmemote_snapshot_log\\tmemote_run_json\\tmemote_run_log\\n' > module_7_gems_manifest.tsv
    printf 'metric\\tvalue\\ntotal_mags\\t0\\nadapted_mags\\t0\\ndoall_only_mags\\t0\\nmemote_snapshot_validated_mags\\t0\\nmemote_run_validated_mags\\t0\\n' > module_7_gems_summary.tsv
    echo 'WRITE_OUTPUT_MANIFEST stub' > write_output_manifest.log
    """
}

// ============================================================================
// WORKFLOW
// ============================================================================

workflow {
    def eggnogGenepredFasta = params.protein_fasta_dir ? absPath(params.protein_fasta_dir) : "${absPath(params.module6_output_dir)}/eggnog/samwise_eggnog.emapper.genepred.fasta"
    def eggnogInputManifest = params.input_manifest ? absPath(params.input_manifest) : "${absPath(params.module6_output_dir)}/summary/eggnog_input_manifest.tsv"
    def selectedMedia = params.media_csv ? absPath(params.media_csv) : absPath(params.media_comprehensive)
    def template = params.template_organism?.toString()?.trim()
    def runAdapt = params.run_gapseq_adapt.toString().toBoolean()
    def runMemote = params.run_memote.toString().toBoolean()
    def memoteMode = params.memote_mode?.toString()?.trim()?.toLowerCase()

    requireRegularFile('Unified Module 6 eggNOG predicted-protein FASTA', eggnogGenepredFasta)
    requireRegularFile('Module 6 eggNOG input manifest', eggnogInputManifest)
    requireRegularFile('Gapfilling media CSV', selectedMedia)
    if (!(template in ['Bacteria', 'Archaea'])) {
        error("template_organism must be exactly Bacteria or Archaea; received: ${template}")
    }
    if (!(memoteMode in ['snapshot', 'run', 'both'])) {
        error("memote_mode must be snapshot, run, or both; received: ${params.memote_mode}")
    }
    if (params.gapseq_extra_args?.toString()?.trim()) {
        error('gapseq_extra_args is retired because gapseq 1.4.0 doall uses a fixed positional interface.')
    }
    if (params.adaptation_compounds?.toString()?.trim()) {
        error('adaptation_compounds is retired; provide MAG-specific values through --adapt_manifest.')
    }

    def upstreamMagIds = parseMagIdsFromManifest(eggnogInputManifest)
    def adaptManifestPath = params.adapt_manifest ? absPath(params.adapt_manifest) : ''
    def adaptValues = [:]
    if (runAdapt) {
        if (!adaptManifestPath) {
            error('run_gapseq_adapt is true but no --adapt_manifest was supplied.')
        }
        requireRegularFile('Adapt manifest', adaptManifestPath)
        adaptValues = parseAdaptManifestTSV(adaptManifestPath, upstreamMagIds)
    } else if (adaptManifestPath) {
        log.warn('An adapt manifest was supplied but --run_gapseq_adapt is false; adaptation will be skipped.')
    }

    log.info("Module 7 parameter resolution: samwise_dir=${params.samwise_dir}; working_dir=${params.working_dir}; results_dir=${params.results_dir}")
    log.info("Module 7 output and publishing: outdir=${params.outdir}; gems_mode=${params.publish_gems_mode}; reports_mode=${params.publish_reports_mode}")
    log.info("Module 7 outdir: ${params.outdir}")
    log.info("Pinned packages: ${params.gapseq_package}; ${params.memote_package}")
    log.info("Prepared protein FASTAs are always written with the .faa suffix.")
    log.info("Media CSV: ${selectedMedia}")
    log.info("Template organism: ${template}")
    log.info("Adapt enabled: ${runAdapt}; manifest entries: ${adaptValues.size()}")
    log.info("MEMOTE enabled: ${runMemote}; mode: ${memoteMode}")

    // Set up and validate the pinned tool environment before preparing FASTA inputs.
    SETUP_GAPSEQ()
    def setupStatus = SETUP_GAPSEQ.out.status
    PREPARE_GAPSEQ_INPUTS(
        channel.value(eggnogGenepredFasta),
        channel.value(eggnogInputManifest),
        setupStatus,
    )

    def preparedInputs = PREPARE_GAPSEQ_INPUTS.out.manifest
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.mag_id, file(row.protein_fasta)) }

    RUN_GAPSEQ_DOALL(
        preparedInputs,
        channel.value(selectedMedia),
        channel.value(template),
        setupStatus,
    )

    def doallGems = RUN_GAPSEQ_DOALL.out.gems
    def adaptInputs = runAdapt \
        ? doallGems.filter { item -> adaptValues.containsKey(item[0]) }.map { item -> tuple(item[0], item[1], adaptValues[item[0]]) } \
        : channel.empty()
    def passThroughGems = runAdapt \
        ? doallGems.filter { item -> !adaptValues.containsKey(item[0]) } \
        : doallGems

    // Always invoke with a possibly empty channel, so downstream output channels exist in every mode.
    RUN_GAPSEQ_ADAPT(adaptInputs, setupStatus)
    def finalGems = passThroughGems.mix(RUN_GAPSEQ_ADAPT.out.gems)
    // Declare variables BEFORE the if/else blocks
    def snapshotInputs
    def runInputs
    
    if (runMemote && memoteMode == 'both') {
        finalGems.multiMap { mag_id, gems_dir, stage ->
            snapshot: tuple(mag_id, gems_dir, stage)
            run: tuple(mag_id, gems_dir, stage)
        }.set { memoteSplit }
        snapshotInputs = memoteSplit.snapshot
        runInputs = memoteSplit.run
    } else if (runMemote && memoteMode == 'snapshot') {
        snapshotInputs = finalGems
        runInputs = channel.empty()
    } else if (runMemote && memoteMode == 'run') {
        snapshotInputs = channel.empty()
        runInputs = finalGems
    } else {
        snapshotInputs = channel.empty()
        runInputs = channel.empty()
    }

    RUN_MEMOTE_SNAPSHOT(snapshotInputs, setupStatus)
    RUN_MEMOTE_RUN(runInputs, setupStatus)

    WRITE_OUTPUT_MANIFEST(
        PREPARE_GAPSEQ_INPUTS.out.manifest,
        RUN_GAPSEQ_DOALL.out.metadata.collect().ifEmpty([]),
        RUN_GAPSEQ_ADAPT.out.metadata.collect().ifEmpty([]),
        RUN_MEMOTE_SNAPSHOT.out.metadata.collect().ifEmpty([]),
        RUN_MEMOTE_RUN.out.metadata.collect().ifEmpty([]),
    )
}
