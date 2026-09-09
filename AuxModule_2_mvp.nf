#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
/*
* AuxModule_2_mvp.nf
*
* Runs selected MVP modules against assemblies from Module 2 and,
* optionally, Module 2b.
*
* Example:
*
*   nextflow run AuxModule_2_mvp.nf \
*       --working_dir /path/to/results \
*       --mvp_modules "0,1,2,3,6" \
*       --threads 24 \
*       --memory_gb 192 \
*       --genomad_db_path /path/to/databases/genomad_db \
*       --checkv_db_path /path/to/databases/checkv_db \
*       --install_databases false
*
* Supported MVP module selections:
*
*   0,1,2,3,4,5,6,7,99,100
*
* Modules are always executed in canonical numerical order, regardless of
* the order supplied to --mvp_modules.
*
* MVP writes directly into:
*
*   <working_dir>/AuxModule_2_mvp/
*
* This preserves MVP's required directory names and allows completed
* outputs to remain available if a later MVP module fails.
*/

params.samwise_dir = params.samwise_dir ?: params.working_dir ?: projectDir
params.working_dir = params.working_dir ?: params.samwise_dir
params.output_dir = null
params.assembly_manifest = null
params.trimmed_manifest = null
params.include_individual_assemblies = true
params.include_coassemblies = true
params.include_subtractive_assemblies = true
params.coassembly_manifest = null
params.coassembly_trimmed_manifest = null
params.subtractive_assembly_manifest = null
params.subtractive_assembly_dir = null
params.subtractive_unmapped_reads_dir = null
params.bbmap_version = "40.02"
params.bbmap_env_dir = null
/*
* MVP module selection.
*
* Examples:
*
*   --mvp_modules "0,1,2,3,6"
*   --mvp_modules "4,5,100"
*   --mvp_modules "6,100"
*
* If earlier modules are omitted, their outputs must already exist in
* params.outdir from a previous run.
*/
params.mvp_modules = "0,1,2,3,4,5,100"
/*
* MVP installation.
*/
params.mvp_version = "1.1.5"
params.auto_install = true
params.tool_env_dir = null
/*
* geNomad and CheckV databases.
*/
params.install_databases = false
params.mvp_database_dir = null
params.genomad_db_path = null
params.checkv_db_path = null
/*
* Resources.
*/
params.threads = 4
params.memory_gb = 0
/*
* MVP Module 00 options.
*/
params.skip_check_errors = false
/*
* MVP Module 01 options.
*/
params.min_seq_size = 0
params.genomad_relaxed = false
params.genomad_conservative = false
params.skip_modify_headers = false
/*
* MVP Module 02 options.
*/
params.viral_min_genes = 1
params.host_viral_genes_ratio = 1
/*
* MVP Module 03 options.
*/
params.min_ani = 95
params.min_tcov = 85
params.min_qcov = 0
params.read_type = "short"
params.unfiltered_protein_file = false
/*
* MVP Module 04 options.
*/
params.interleaved = true
params.delete_mapping_intermediates = true
/*
* MVP Module 05 options.
*/
params.covered_fraction = null
params.normalization = "RPKM"
params.filtration = "conservative"
/*
* MVP Module 06 options.
*/
params.functional_fasta_files = "representative"
params.phrogs_evalue = 0.01
params.phrogs_score = 60
params.pfam_evalue = 0.01
params.pfam_score = 50
params.functional_ads = true
params.ads_evalue = 0.01
params.ads_score = 60
params.functional_rdrp = false
params.rdrp_evalue = 0.001
params.rdrp_score = 50
params.functional_dram = true
params.delete_functional_intermediates = true
/*
* MVP Module 07 options.
*/
params.binning_sample_group = null
params.read_mapping_sample_group = null
params.keep_bam = false
params.delete_binning_intermediates = false
/*
* MVP Module 99 options.
*
* Required only when 99 is included in --mvp_modules.
*
* Valid steps:
*   setup_metadata
*   prep_submission
*
* A template is required for prep_submission.
*/
params.miuvig_identifier = null
params.miuvig_step = null
params.miuvig_template = null
params.force = false
/*
* Convert a path to an absolute normalized path.
*/
def absPath(value) {
    def text = value == null
        ? ""
        : value.toString().trim()
    if (!text || text == "null" || text == "NA") {
        return ""
    }
    return java.nio.file.Paths
        .get(text)
        .toAbsolutePath()
        .normalize()
        .toString()
}
/*
* Derived directories.
*/
params.results_dir = params.working_dir
    ? absPath(params.working_dir)
    : (params.output_dir
        ? absPath(params.output_dir)
        : absPath("."))
params.module1_outdir = "${params.results_dir}/module_1_readtrimming"
params.module2_outdir = "${params.results_dir}/module_2_readassembly"
params.module2b_outdir = "${params.results_dir}/module_2b_coassembly"
params.module5_outdir = "${params.results_dir}/module_5_subassembly"
params.outdir = "${params.results_dir}/AuxModule_2_mvp"
params.generated_interleaved_reads_dir = "${params.outdir}/interleaved_reads"
/*
* Parse a comma-, semicolon-, or whitespace-separated module list.
*/
def parseMvpModules(value) {
    def aliases = ["00": "0", "01": "1", "02": "2", "03": "3", "04": "4", "05": "5", "06": "6", "07": "7", "099": "99"]
    return value
        .toString()
        .split(/[,;\s]+/)
        .collect { item -> item.trim() }
        .findAll { item -> item }
        .collect { item -> aliases.containsKey(item) ? aliases[item] : item }
        .unique()
}
workflow {
    def supported_modules = ["0", "1", "2", "3", "4", "5", "6", "7", "99", "100"]
    def canonical_order = ["0", "1", "2", "3", "4", "5", "6", "7", "99", "100"]

    def selected_modules = parseMvpModules(params.mvp_modules)

    if (!selected_modules) {
        error(
            """
            No MVP modules were selected.
            Example:
              --mvp_modules "0,1,2,3,6"
            """.stripIndent()
        )
    }

    def invalid_modules = selected_modules.findAll { module ->
        !(module in supported_modules)
    }

    if (invalid_modules) {
        error(
            """
            Unsupported MVP module selection:
              ${invalid_modules.join(', ')}
            Supported modules:
              ${supported_modules.join(', ')}
            """.stripIndent()
        )
    }

    def ordered_modules = canonical_order.findAll { module ->
        module in selected_modules
    }

    if (!(params.read_type.toString() in ["short", "long"])) {
        error(
            """
            Invalid --read-type value: ${params.read_type}
            Supported values:
              short
              long
            """.stripIndent()
        )
    }

    if (!(params.normalization.toString() in ["RPKM", "FPKM"])) {
        error(
            """
            Invalid --normalization value: ${params.normalization}
            Supported values:
              RPKM
              FPKM
            """.stripIndent()
        )
    }

    if (!(params.filtration.toString() in ["relaxed", "conservative"])) {
        error(
            """
            Invalid --filtration value: ${params.filtration}
            Supported values:
              relaxed
              conservative
            """.stripIndent()
        )
    }

    if (params.genomad_relaxed.toString().toBoolean() && params.genomad_conservative.toString().toBoolean()) {
        error(
            """
            --genomad_relaxed and --genomad_conservative cannot both be enabled.
            """.stripIndent()
        )
    }

    if ("99" in selected_modules) {
        if (!params.miuvig_identifier) {
            error(
                """
                MVP Module 99 was selected, but no MIUViG identifier was supplied.
                Provide:
                  --miuvig_identifier IDENTIFIER
                """.stripIndent()
            )
        }

        if (!params.miuvig_step || !(params.miuvig_step.toString() in ["setup_metadata", "prep_submission"])) {
            error(
                """
                MVP Module 99 requires either:
                  --miuvig_step setup_metadata
                or:
                  --miuvig_step prep_submission
                """.stripIndent()
            )
        }

        if (params.miuvig_step.toString() == "prep_submission" && !params.miuvig_template) {
            error(
                """
                MVP Module 99 prep_submission requires:
                  --miuvig_template /path/to/template.tsv
                """.stripIndent()
            )
        }
    }

    def include_individual = params.include_individual_assemblies
        .toString()
        .toBoolean()

    def include_coassemblies = params.include_coassemblies
        .toString()
        .toBoolean()

    def include_subtractive = params.include_subtractive_assemblies
        .toString()
        .toBoolean()

    if (!include_individual && !include_coassemblies && !include_subtractive) {
        error(
            """
            No assembly classes were enabled.
            Enable at least one of:
              --include_individual_assemblies true
              --include_coassemblies true
              --include_subtractive_assemblies true
            """.stripIndent()
        )
    }

    def assembly_manifest_file = params.assembly_manifest
        ? absPath(params.assembly_manifest)
        : "${params.module2_outdir}/summary/assembly_manifest.tsv"

    def trimmed_manifest_file = params.trimmed_manifest
        ? absPath(params.trimmed_manifest)
        : "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    def coassembly_manifest_file = params.coassembly_manifest
        ? absPath(params.coassembly_manifest)
        : "${params.module2b_outdir}/summary/assembly_manifest.tsv"

    def coassembly_reads_file = params.coassembly_trimmed_manifest
        ? absPath(params.coassembly_trimmed_manifest)
        : "${params.module2b_outdir}/summary/coassembly_trimmed_manifest.tsv"

    def subtractive_manifest_file = params.subtractive_assembly_manifest
        ? absPath(params.subtractive_assembly_manifest)
        : "${params.module5_outdir}/summary/subtractive_assembly_summary_manifest.tsv"

    def subtractive_assembly_dir = params.subtractive_assembly_dir
        ? absPath(params.subtractive_assembly_dir)
        : "${params.module5_outdir}/assemblies"

    def subtractive_unmapped_reads_dir = params.subtractive_unmapped_reads_dir
        ? absPath(params.subtractive_unmapped_reads_dir)
        : "${params.module5_outdir}/unmapped_reads"

    log.info("Auxiliary Module 2 MVP output directory: ${params.outdir}")
    log.info("Selected MVP modules: ${ordered_modules.join(',')}")
    log.info("MVP modules will execute in this order: ${ordered_modules.join(' -> ')}")

    SETUP_AUXMODULE2_MVP()

    DETERMINE_BBTOOLS_REQUIREMENT(
        channel.value(include_individual),
        channel.value(include_coassemblies),
        channel.value(include_subtractive),
        channel.value(trimmed_manifest_file),
        channel.value(coassembly_reads_file),
        channel.value(subtractive_manifest_file),
        channel.value(subtractive_unmapped_reads_dir),
        channel.value(params.results_dir),
        SETUP_AUXMODULE2_MVP.out.status,
    )

    SETUP_BBTOOLS(
        DETERMINE_BBTOOLS_REQUIREMENT.out.status
    )

    PREPARE_MVP_READS(
        channel.value(include_individual),
        channel.value(include_coassemblies),
        channel.value(include_subtractive),
        channel.value(trimmed_manifest_file),
        channel.value(coassembly_reads_file),
        channel.value(subtractive_manifest_file),
        channel.value(subtractive_unmapped_reads_dir),
        channel.value(params.results_dir),
        SETUP_AUXMODULE2_MVP.out.status,
        SETUP_BBTOOLS.out.status,
    )

    PREPARE_MVP_METADATA(
        channel.value(include_individual),
        channel.value(include_coassemblies),
        channel.value(include_subtractive),
        channel.value(assembly_manifest_file),
        PREPARE_MVP_READS.out.reads_manifest,
        channel.value(coassembly_manifest_file),
        channel.value(subtractive_manifest_file),
        channel.value(subtractive_assembly_dir),
        channel.value(params.results_dir),
        SETUP_AUXMODULE2_MVP.out.status,
    )

    RUN_MVP(
        PREPARE_MVP_METADATA.out.metadata,
        PREPARE_MVP_METADATA.out.input_summary,
        SETUP_AUXMODULE2_MVP.out.status,
    )
}

/*
* Inspect enabled read manifests to determine whether any paired read sets
* must be interleaved before MVP metadata can be generated.
*/
process DETERMINE_BBTOOLS_REQUIREMENT {
    tag "determine_bbtools_requirement"
    publishDir "${params.outdir}/summary", mode: "copy", pattern: "bbtools_requirement.env"
    publishDir "${params.outdir}/logs", mode: "copy", pattern: "determine_bbtools_requirement.log"

    input:
    val include_individual
    val include_coassemblies
    val include_subtractive
    val trimmed_manifest
    val coassembly_trimmed_manifest
    val subtractive_manifest
    val subtractive_unmapped_reads_dir
    val base_dir
    path mvp_tools_status

    output:
    path "bbtools_requirement.env", emit: status
    path "determine_bbtools_requirement.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG="determine_bbtools_requirement.log"
    STATUS="bbtools_requirement.env"

    TOOL_ENV="\$(
        grep '^TOOL_ENV=' "${mvp_tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"
    PYTHON_EXECUTABLE="\$(
        grep '^PYTHON_EXECUTABLE=' "${mvp_tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    if [[ -z "\$PYTHON_EXECUTABLE" || ! -x "\$PYTHON_EXECUTABLE" ]]; then
        echo "ERROR: Python executable was not found in the MVP environment." >&2        
        cat "${mvp_tools_status}" >&2 || true
        exit 1
    fi

    echo "BBTools requirement analysis started: \$(date)" > "\$LOG"
    echo "Path resolution base directory: ${base_dir}" >> "\$LOG"
    echo "Include individual assemblies: ${include_individual}" >> "\$LOG"
    echo "Include coassemblies: ${include_coassemblies}" >> "\$LOG"
    echo "Include subtractive assemblies: ${include_subtractive}" >> "\$LOG"
    echo "Module 1 trimmed-read manifest: ${trimmed_manifest}" >> "\$LOG"
    echo "Module 2b trimmed-read manifest: ${coassembly_trimmed_manifest}" >> "\$LOG"
    echo "Module 5 subtractive assembly manifest: ${subtractive_manifest}" >> "\$LOG"
    echo "Module 5 unmapped-read directory: ${subtractive_unmapped_reads_dir}" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    "\$PYTHON_EXECUTABLE" - \\
        "${include_individual}" \\
        "${include_coassemblies}" \\
        "${include_subtractive}" \\
        "${trimmed_manifest}" \\
        "${coassembly_trimmed_manifest}" \\
        "${subtractive_manifest}" \\
        "${subtractive_unmapped_reads_dir}" \\
        "${base_dir}" \\
        "\$STATUS" \\
        "\$LOG" <<'PY'
import csv
import sys
from pathlib import Path
(
    include_individual,
    include_coassemblies,
    include_subtractive,
    trimmed_manifest,
    coassembly_trimmed_manifest,
    subtractive_manifest,
    subtractive_unmapped_reads_dir,
    base_dir,
    status_file,
    log_file,
) = sys.argv[1:]

include_individual = include_individual.lower() == "true"
include_coassemblies = include_coassemblies.lower() == "true"
include_subtractive = include_subtractive.lower() == "true"

base_dir = Path(base_dir).resolve()
status_file = Path(status_file)
log_file = Path(log_file)

def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)

def resolve_path(value):
    value = str(value or "").strip()
    if not value:
        return ""
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return str(path.resolve())

def read_tsv(path_text, required_columns=None, required=True):
    required_columns = required_columns or []
    path = Path(resolve_path(path_text))
    if not path.is_file():
        if required:
            raise RuntimeError(f"Required file does not exist: {path}")
        return []
    if path.stat().st_size == 0:
        if required:
            raise RuntimeError(f"Required file is empty: {path}")
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=chr(9))
        observed_columns = set(reader.fieldnames or [])
        missing = set(required_columns) - observed_columns
        if missing:
            raise RuntimeError(
                f"File {path} is missing required columns: {sorted(missing)}"
            )
        return [
            row
            for row in reader
            if any(str(value or "").strip() for value in row.values())
        ]

def file_exists_nonempty(path_text):
    path = Path(resolve_path(path_text))
    return path.is_file() and path.stat().st_size > 0

def first_existing_column(rows, candidates):
    if not rows:
        return None
    available = set(rows[0].keys())
    for candidate in candidates:
        if candidate in available:
            return candidate
    return None

def detect_read_layout(row, label, sample_label, fallback_interleaved=""):
    layout = str(row.get("layout", "")).strip().lower()

    paired_candidates = [
        ("read1", "read2"),
        ("r1", "r2"),
        ("forward", "reverse"),
        ("fq1", "fq2"),
    ]
    interleaved_candidates = [
        "interleaved",
        "read_path",
        "reads",
        "reads_path",
        "unmapped_reads",
        "unmapped_read_path",
    ]

    def get_pair():
        for c1, c2 in paired_candidates:
            v1 = str(row.get(c1, "")).strip()
            v2 = str(row.get(c2, "")).strip()
            if v1 or v2:
                if not (v1 and v2):
                    raise RuntimeError(
                        f"{label} sample {sample_label} has incomplete paired-read fields "
                        f"({c1}={v1!r}, {c2}={v2!r})"
                    )
                return v1, v2
        return "", ""

    def get_interleaved():
        for column in interleaved_candidates:
            value = str(row.get(column, "")).strip()
            if value:
                return value
        return fallback_interleaved

    if layout:
        if layout == "paired":
            r1, r2 = get_pair()
            if not (r1 and r2):
                raise RuntimeError(
                    f"{label} sample {sample_label} declares layout=paired but no complete read pair was found"
                )
            return "paired"
        if layout == "interleaved":
            inter = get_interleaved()
            if not inter:
                raise RuntimeError(
                    f"{label} sample {sample_label} declares layout=interleaved but no interleaved read path was found"
                )
            return "interleaved"
        raise RuntimeError(
            f"{label} sample {sample_label} has unsupported layout value: {layout!r}"
        )

    r1, r2 = get_pair()
    if r1 and r2:
        return "paired"

    inter = get_interleaved()
    if inter:
        return "interleaved"

    raise RuntimeError(
        f"Could not determine read layout for {label} sample {sample_label}"
    )

needs_bbtools = False
inspected_counts = {
    "individual": 0,
    "coassembly": 0,
    "subtractive": 0,
}

if include_individual:
    rows = read_tsv(trimmed_manifest, required_columns=["sample_id"], required=True)
    for row in rows:
        sample_id = str(row.get("sample_id", "")).strip() or "UNKNOWN"
        layout = detect_read_layout(
            row,
            "Module 1 trimmed-read manifest",
            sample_id,
        )
        inspected_counts["individual"] += 1
        if layout == "paired":
            needs_bbtools = True

if include_coassemblies:
    if not file_exists_nonempty(coassembly_trimmed_manifest):
        log(
            "Module 2b coassembly trimmed-read manifest was not found or was empty. "
            "Proceeding without coassemblies."
        )
    else:
        rows = read_tsv(coassembly_trimmed_manifest, required_columns=[], required=True)
        id_column = first_existing_column(
            rows,
            ["sample_id", "coassembly_id", "coassembly", "sample"],
        )
        if not id_column:
            raise RuntimeError(
                "Could not identify the coassembly read sample ID column"
            )
        for row in rows:
            sample_id = str(row.get(id_column, "")).strip() or "UNKNOWN"
            layout = detect_read_layout(
                row,
                "Module 2b coassembly trimmed-read manifest",
                sample_id,
            )
            inspected_counts["coassembly"] += 1
            if layout == "paired":
                needs_bbtools = True

if include_subtractive:
    if not file_exists_nonempty(subtractive_manifest):
        log(
            "Module 5 subtractive assembly manifest was not found or was empty. "
            "Proceeding without subtractive assemblies."
        )
    else:
        rows = read_tsv(subtractive_manifest, required_columns=[], required=True)
        sample_column = first_existing_column(
            rows,
            ["sample_id", "safe_sample_id", "sample"],
        )
        if not sample_column:
            raise RuntimeError(
                "Could not identify the subtractive sample ID column"
            )
        subtractive_dir = Path(resolve_path(subtractive_unmapped_reads_dir))
        for row in rows:
            sample_id = str(row.get(sample_column, "")).strip() or "UNKNOWN"
            safe_sample_id = str(row.get("safe_sample_id", "")).strip() or sample_id
            fallback_interleaved = str(
                (subtractive_dir / f"{safe_sample_id}.unmapped_interleaved.fastq.gz").resolve()
            )
            layout = detect_read_layout(
                row,
                "Module 5 subtractive assembly manifest",
                sample_id,
                fallback_interleaved=fallback_interleaved,
            )
            inspected_counts["subtractive"] += 1
            if layout == "paired":
                needs_bbtools = True

with status_file.open("w") as handle:
    print(f"NEEDS_BBTOOLS={'true' if needs_bbtools else 'false'}", file=handle)
    print(f"INDIVIDUAL_ROWS={inspected_counts['individual']}", file=handle)
    print(f"COASSEMBLY_ROWS={inspected_counts['coassembly']}", file=handle)
    print(f"SUBTRACTIVE_ROWS={inspected_counts['subtractive']}", file=handle)
    print("STATUS=complete", file=handle)

log("----------------------------------------")
log(f"Individual read rows inspected: {inspected_counts['individual']}")
log(f"Coassembly read rows inspected: {inspected_counts['coassembly']}")
log(f"Subtractive read rows inspected: {inspected_counts['subtractive']}")
log(f"BBTools required: {needs_bbtools}")
PY

    echo "BBTools requirement analysis finished: \$(date)" >> "\$LOG"
    """
}

/*
* Prepare BBMap / BBTools only when at least one enabled input class
* contains paired reads that must be interleaved for MVP.
*/
process SETUP_BBTOOLS {
    tag "setup_bbmap_bbtools"
    cache false
    publishDir "${params.outdir}/setup", mode: "copy", pattern: "bbtools_status.env"

    input:
    path requirement_status

    output:
    path "bbtools_status.env", emit: status

    script:
    def env_dir = params.bbmap_env_dir
        ? absPath(params.bbmap_env_dir)
        : absPath("${params.outdir}/conda_envs/bbmap")
    """
    set -euo pipefail
    STATUS_FILE="bbtools_status.env"
    TOOL_ENV="${env_dir}"
    NEEDS_BBTOOLS="\$(
        grep '^NEEDS_BBTOOLS=' "${requirement_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"
    echo "BBMap / BBTools setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested BBMap version: ${params.bbmap_version}" >> "\$STATUS_FILE"
    echo "Requested Conda environment: \$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Interleaving required: \$NEEDS_BBTOOLS" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    if [[ "\$NEEDS_BBTOOLS" != "true" ]]; then
        echo "BBTOOLS_REQUIRED=false" >> "\$STATUS_FILE"
        echo "BBMAP_ENV=SYSTEM" >> "\$STATUS_FILE"
        echo "REFORMAT_SH=" >> "\$STATUS_FILE"
        echo "BBMap / BBTools setup skipped because no paired reads require interleaving." \
            >> "\$STATUS_FILE"
        echo "BBMap / BBTools setup finished: \$(date)" \
            >> "\$STATUS_FILE"
        exit 0
    fi

    validate_reformat() {
        local executable="\$1"
        local python_executable="\$2"

        [[ -x "\$executable" ]] || return 1
        [[ -x "\$python_executable" ]] || return 1

        "\$executable" -h >> "\$STATUS_FILE" 2>&1 || return 1
        "\$python_executable" --version >> "\$STATUS_FILE" 2>&1 || return 1

        return 0
    }

    existing_reformat="\$TOOL_ENV/bin/reformat.sh"
    existing_python="\$TOOL_ENV/bin/python"
    if [[ -d "\$TOOL_ENV" ]]; then
        echo "Existing BBMap Conda environment detected." \
            >> "\$STATUS_FILE"
        if validate_reformat "\$existing_reformat" "\$existing_python"; then
            echo "Existing BBMap environment passed validation." \
                >> "\$STATUS_FILE"
            echo "BBTOOLS_REQUIRED=true" \
                >> "\$STATUS_FILE"
            echo "BBMAP_ENV=\$TOOL_ENV" \
                >> "\$STATUS_FILE"
            echo "REFORMAT_SH=\$existing_reformat" \
                >> "\$STATUS_FILE"
            echo "PYTHON_EXECUTABLE=\$existing_python" \
                >> "\$STATUS_FILE"
            echo "BBMap / BBTools setup finished: \$(date)" \
                >> "\$STATUS_FILE"
            exit 0
        fi
        echo "Existing BBMap environment is incomplete or invalid." \
            >> "\$STATUS_FILE"
        echo "Removing invalid environment: \$TOOL_ENV" \
            >> "\$STATUS_FILE"
        rm -rf "\$TOOL_ENV"
    fi
    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: BBMap environment does not exist and auto-install is disabled." \
            >> "\$STATUS_FILE"
        echo "Expected executable: \$existing_reformat" \
            >> "\$STATUS_FILE"
        exit 1
    fi
    INSTALLER=""
    if command -v mamba >/dev/null 2>&1; then
        INSTALLER="mamba"
        echo "Using mamba: \$(command -v mamba)" \
            >> "\$STATUS_FILE"
    elif command -v conda >/dev/null 2>&1; then
        INSTALLER="conda"
        echo "Using conda: \$(command -v conda)" \
            >> "\$STATUS_FILE"
    else
        echo "ERROR: Neither mamba nor conda was found in PATH." \
            >> "\$STATUS_FILE"
        echo "Install mamba/conda, or create the BBMap environment manually." \
            >> "\$STATUS_FILE"
        exit 1
    fi
    mkdir -p "\$(dirname "\$TOOL_ENV")"
    echo "Creating BBMap Conda environment: \$TOOL_ENV" \
        >> "\$STATUS_FILE"
    "\$INSTALLER" create -y \
        -p "\$TOOL_ENV" \
        -c conda-forge \
        -c bioconda \
        "python=3.11" \
        "bbmap=${params.bbmap_version}" \
        >> "\$STATUS_FILE" 2>&1
    if ! validate_reformat "\$existing_reformat" "\$existing_python"; then
        echo "ERROR: BBMap installation completed but reformat.sh failed validation." \
            >> "\$STATUS_FILE"
        echo "Expected executable: \$existing_reformat" \
            >> "\$STATUS_FILE"
        exit 1
    fi
    echo "BBTOOLS_REQUIRED=true" \
        >> "\$STATUS_FILE"
    echo "BBMAP_ENV=\$TOOL_ENV" \
        >> "\$STATUS_FILE"
    echo "REFORMAT_SH=\$existing_reformat" \
        >> "\$STATUS_FILE"
    echo "PYTHON_EXECUTABLE=\$existing_python" \
        >> "\$STATUS_FILE"
    echo "BBMap / BBTools setup finished: \$(date)" \
        >> "\$STATUS_FILE"
    """
}

/*
* Prepare all enabled read classes for MVP as interleaved FASTQ.
*
* Existing interleaved reads are used directly from the corresponding
* manifests. Paired R1/R2 reads are interleaved record-by-record into:
*
*   AuxModule_2_mvp/interleaved_reads
*/
process PREPARE_MVP_READS {
    tag "prepare_reads_for_mvp"
    publishDir "${params.outdir}/metadata", mode: "copy", pattern: "mvp_prepared_reads_manifest.tsv"
    publishDir "${params.outdir}/summary", mode: "copy", pattern: "mvp_read_preparation_stats.tsv"
    publishDir "${params.outdir}/logs", mode: "copy", pattern: "prepare_mvp_reads.log"

    input:
    val include_individual
    val include_coassemblies
    val include_subtractive
    val trimmed_manifest
    val coassembly_trimmed_manifest
    val subtractive_manifest
    val subtractive_unmapped_reads_dir
    val base_dir
    path mvp_tools_status
    path bbtools_status

    output:
    path "mvp_prepared_reads_manifest.tsv", emit: reads_manifest
    path "mvp_read_preparation_stats.tsv", emit: stats
    path "mvp_read_preparation_status.env", emit: status
    path "prepare_mvp_reads.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG="prepare_mvp_reads.log"
    READS_MANIFEST="mvp_prepared_reads_manifest.tsv"
    STATS="mvp_read_preparation_stats.tsv"
    STATUS="mvp_read_preparation_status.env"
    OUTPUT_DIR="${params.generated_interleaved_reads_dir}"

    MVP_TOOL_ENV="\$(
        grep '^TOOL_ENV=' "${mvp_tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"
    MVP_PYTHON="\$(
        grep '^PYTHON_EXECUTABLE=' "${mvp_tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"
    BBMAP_ENV="\$(
        grep '^BBMAP_ENV=' "${bbtools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"
    REFORMAT_SH="\$(
        grep '^REFORMAT_SH=' "${bbtools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"
    BBTOOLS_REQUIRED="\$(
        grep '^BBTOOLS_REQUIRED=' "${bbtools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"

    if [[ -n "\$MVP_TOOL_ENV" && "\$MVP_TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$MVP_TOOL_ENV/bin:\$PATH"
    fi
    if [[ -n "\$BBMAP_ENV" && "\$BBMAP_ENV" != "SYSTEM" ]]; then
        export PATH="\$BBMAP_ENV/bin:\$PATH"
    fi
    if [[ -z "\$MVP_PYTHON" || ! -x "\$MVP_PYTHON" ]]; then
        echo "ERROR: Python executable was not found in the MVP environment." >&2
        cat "${mvp_tools_status}" >&2 || true
        exit 1
    fi

    echo "MVP read preparation started: \$(date)" > "\$LOG"
    echo "Include individual assemblies: ${include_individual}" >> "\$LOG"
    echo "Include coassemblies: ${include_coassemblies}" >> "\$LOG"
    echo "Include subtractive assemblies: ${include_subtractive}" >> "\$LOG"
    echo "Module 1 trimmed-read manifest: ${trimmed_manifest}" >> "\$LOG"
    echo "Module 2b trimmed-read manifest: ${coassembly_trimmed_manifest}" >> "\$LOG"
    echo "Module 5 subtractive assembly manifest: ${subtractive_manifest}" >> "\$LOG"
    echo "Module 5 unmapped-read directory: ${subtractive_unmapped_reads_dir}" >> "\$LOG"
    echo "Interleaved-read output directory: \$OUTPUT_DIR" >> "\$LOG"
    echo "Path resolution base directory: ${base_dir}" >> "\$LOG"
    echo "BBTools required: \$BBTOOLS_REQUIRED" >> "\$LOG"
    echo "Task directory: \$(pwd -P)" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    mkdir -p "\$OUTPUT_DIR"

    "\$MVP_PYTHON" - \\
        "${include_individual}" \\
        "${include_coassemblies}" \\
        "${include_subtractive}" \\
        "${trimmed_manifest}" \\
        "${coassembly_trimmed_manifest}" \\
        "${subtractive_manifest}" \\
        "${subtractive_unmapped_reads_dir}" \\
        "${base_dir}" \\
        "\$OUTPUT_DIR" \\
        "\$REFORMAT_SH" \\
        "\$READS_MANIFEST" \\
        "\$STATS" \\
        "\$LOG" <<'PY'
import csv
import re
import subprocess
import sys
from pathlib import Path

(
    include_individual,
    include_coassemblies,
    include_subtractive,
    trimmed_manifest,
    coassembly_trimmed_manifest,
    subtractive_manifest,
    subtractive_unmapped_reads_dir,
    base_dir,
    output_dir,
    reformat_sh,
    reads_manifest,
    stats_file,
    log_file,
) = sys.argv[1:]

include_individual = include_individual.lower() == "true"
include_coassemblies = include_coassemblies.lower() == "true"
include_subtractive = include_subtractive.lower() == "true"

base_dir = Path(base_dir).resolve()
output_dir = Path(output_dir)
reads_manifest = Path(reads_manifest)
stats_file = Path(stats_file)
log_file = Path(log_file)

output_dir.mkdir(parents=True, exist_ok=True)

def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)

def resolve_path(value):
    value = str(value or "").strip()
    if not value:
        return ""
    path = Path(value)
    if not path.is_absolute():
        path = base_dir / path
    return str(path.resolve())

def safe_name(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("._-")
    return value or "unnamed"

def read_tsv(path_text, required_columns=None, required=True):
    required_columns = required_columns or []
    path = Path(resolve_path(path_text))
    if not path.is_file():
        if required:
            raise RuntimeError(f"Required file does not exist: {path}")
        return []
    if path.stat().st_size == 0:
        if required:
            raise RuntimeError(f"Required file is empty: {path}")
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=chr(9))
        observed_columns = set(reader.fieldnames or [])
        missing = set(required_columns) - observed_columns
        if missing:
            raise RuntimeError(
                f"File {path} is missing required columns: {sorted(missing)}"
            )
        return [
            row
            for row in reader
            if any(str(value or "").strip() for value in row.values())
        ]

def file_exists_nonempty(path_text):
    path = Path(resolve_path(path_text))
    return path.is_file() and path.stat().st_size > 0

def first_existing_column(rows, candidates):
    if not rows:
        return None
    available = set(rows[0].keys())
    for candidate in candidates:
        if candidate in available:
            return candidate
    return None

def detect_read_spec(row, label, sample_label, fallback_interleaved=""):
    layout = str(row.get("layout", "")).strip().lower()

    paired_candidates = [
        ("read1", "read2"),
        ("r1", "r2"),
        ("forward", "reverse"),
        ("fq1", "fq2"),
    ]
    interleaved_candidates = [
        "interleaved",
        "read_path",
        "reads",
        "reads_path",
        "unmapped_reads",
        "unmapped_read_path",
    ]

    def get_pair():
        for c1, c2 in paired_candidates:
            v1 = str(row.get(c1, "")).strip()
            v2 = str(row.get(c2, "")).strip()
            if v1 or v2:
                if not (v1 and v2):
                    raise RuntimeError(
                        f"{label} sample {sample_label} has incomplete paired-read fields "
                        f"({c1}={v1!r}, {c2}={v2!r})"
                    )
                return v1, v2
        return "", ""

    def get_interleaved():
        for column in interleaved_candidates:
            value = str(row.get(column, "")).strip()
            if value:
                return value
        return fallback_interleaved

    if layout:
        if layout == "paired":
            r1, r2 = get_pair()
            if not (r1 and r2):
                raise RuntimeError(
                    f"{label} sample {sample_label} declares layout=paired but no complete read pair was found"
                )
            return {
                "layout": "paired",
                "read1": r1,
                "read2": r2,
                "interleaved": "",
            }
        if layout == "interleaved":
            inter = get_interleaved()
            if not inter:
                raise RuntimeError(
                    f"{label} sample {sample_label} declares layout=interleaved but no interleaved read path was found"
                )
            return {
                "layout": "interleaved",
                "read1": "",
                "read2": "",
                "interleaved": inter,
            }
        raise RuntimeError(
            f"{label} sample {sample_label} has unsupported layout value: {layout!r}"
        )

    r1, r2 = get_pair()
    if r1 and r2:
        return {
            "layout": "paired",
            "read1": r1,
            "read2": r2,
            "interleaved": "",
        }

    inter = get_interleaved()
    if inter:
        return {
            "layout": "interleaved",
            "read1": "",
            "read2": "",
            "interleaved": inter,
        }

    raise RuntimeError(
        f"Could not determine read layout for {label} sample {sample_label}"
    )

def require_nonempty_file(path_text, label):
    resolved = Path(resolve_path(path_text))
    if not resolved.is_file():
        raise RuntimeError(f"{label} does not exist: {resolved}")
    if resolved.stat().st_size == 0:
        raise RuntimeError(f"{label} is empty: {resolved}")
    return str(resolved.resolve())

prepared_rows = []
stats_rows = []
used_output_names = set()

def prepare_record(source_class, sample_id, safe_sample_id, row, label, fallback_interleaved=""):
    sample_id = str(sample_id or "").strip()
    safe_sample_id = str(safe_sample_id or "").strip() or sample_id

    spec = detect_read_spec(
        row,
        label,
        sample_id,
        fallback_interleaved=fallback_interleaved,
    )

    action = ""
    source_read1 = ""
    source_read2 = ""
    source_interleaved = ""

    if spec["layout"] == "interleaved":
        source_interleaved = require_nonempty_file(
            spec["interleaved"],
            f"Interleaved reads for {source_class} sample {sample_id}",
        )
        read_path = source_interleaved
        action = "existing_interleaved"
        log(
            f"Using existing interleaved reads for {source_class} sample {sample_id}: "
            f"{read_path}"
        )
    else:
        source_read1 = require_nonempty_file(
            spec["read1"],
            f"Read1 for {source_class} sample {sample_id}",
        )
        source_read2 = require_nonempty_file(
            spec["read2"],
            f"Read2 for {source_class} sample {sample_id}",
        )
        if not reformat_sh or not Path(reformat_sh).is_file():
            raise RuntimeError(
                f"Paired reads were detected for {source_class} sample {sample_id}, "
                "but reformat.sh is unavailable"
            )
        output_name = (
            f"{source_class}.{safe_name(safe_sample_id)}."
            f"mvp_interleaved.fastq.gz"
        )
        if output_name in used_output_names:
            raise RuntimeError(
                f"Multiple paired read sets would produce the same output file: {output_name}"
            )
        used_output_names.add(output_name)
        destination = output_dir / output_name
        command = [
            str(Path(reformat_sh).resolve()),
            f"in1={source_read1}",
            f"in2={source_read2}",
            f"out={destination}",
            "verifypaired=t",
            "overwrite=t",
        ]
        log(
            f"Interleaving paired reads for {source_class} sample {sample_id}"
        )
        log(
            "Command: " + " ".join(map(str, command))
        )
        with log_file.open("a") as log_handle:
            subprocess.run(
                command,
                check=True,
                stdout=log_handle,
                stderr=log_handle,
            )
        if not destination.is_file() or destination.stat().st_size == 0:
            raise RuntimeError(
                f"reformat.sh did not produce a non-empty output for "
                f"{source_class} sample {sample_id}: {destination}"
            )
        read_path = str(destination.resolve())
        action = "paired_to_interleaved"
        log(
            f"Created interleaved reads for {source_class} sample {sample_id}: "
            f"{read_path}"
        )

    manifest_row = {
        "source_class": source_class,
        "sample_id": sample_id,
        "safe_sample_id": safe_sample_id,
        "read_path": read_path,
        "action": action,
    }
    stats_row = {
        "source_class": source_class,
        "sample_id": sample_id,
        "safe_sample_id": safe_sample_id,
        "source_layout": spec["layout"],
        "source_read1": source_read1,
        "source_read2": source_read2,
        "source_interleaved": source_interleaved,
        "mvp_read_path": read_path,
        "action": action,
    }
    prepared_rows.append(manifest_row)
    stats_rows.append(stats_row)

if include_individual:
    rows = read_tsv(trimmed_manifest, required_columns=["sample_id"], required=True)
    for row in rows:
        sample_id = str(row.get("sample_id", "")).strip()
        safe_sample_id = str(row.get("safe_sample_id", "")).strip() or sample_id
        if not sample_id:
            raise RuntimeError(
                "Encountered Module 1 trimmed-read manifest row with an empty sample_id"
            )
        prepare_record(
            "individual",
            sample_id,
            safe_sample_id,
            row,
            "Module 1 trimmed-read manifest",
        )

if include_coassemblies:
    if not file_exists_nonempty(coassembly_trimmed_manifest):
        log(
            "Module 2b coassembly trimmed-read manifest was not found or was empty. "
            "Proceeding without coassemblies."
        )
    else:
        rows = read_tsv(coassembly_trimmed_manifest, required_columns=[], required=True)
        id_column = first_existing_column(
            rows,
            ["sample_id", "coassembly_id", "coassembly", "sample"],
        )
        if not id_column:
            raise RuntimeError(
                "Could not identify the coassembly read sample ID column"
            )
        for row in rows:
            sample_id = str(row.get(id_column, "")).strip()
            if not sample_id:
                raise RuntimeError(
                    "Encountered Module 2b coassembly trimmed-read manifest row with an empty sample ID"
                )
            prepare_record(
                "coassembly",
                sample_id,
                sample_id,
                row,
                "Module 2b coassembly trimmed-read manifest",
            )

if include_subtractive:
    if not file_exists_nonempty(subtractive_manifest):
        log(
            "Module 5 subtractive assembly manifest was not found or was empty. "
            "Proceeding without subtractive assemblies."
        )
    else:
        rows = read_tsv(subtractive_manifest, required_columns=[], required=True)
        sample_column = first_existing_column(
            rows,
            ["sample_id", "safe_sample_id", "sample"],
        )
        if not sample_column:
            raise RuntimeError(
                "Could not identify the subtractive sample ID column"
            )
        subtractive_dir = Path(resolve_path(subtractive_unmapped_reads_dir))
        for row in rows:
            sample_id = str(row.get(sample_column, "")).strip()
            if not sample_id:
                raise RuntimeError(
                    "Encountered Module 5 subtractive manifest row with an empty sample ID"
                )
            safe_sample_id = str(row.get("safe_sample_id", "")).strip() or sample_id
            fallback_interleaved = str(
                (subtractive_dir / f"{safe_sample_id}.unmapped_interleaved.fastq.gz").resolve()
            )
            prepare_record(
                "subtractive",
                sample_id,
                safe_sample_id,
                row,
                "Module 5 subtractive assembly manifest",
                fallback_interleaved=fallback_interleaved,
            )

with reads_manifest.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter=chr(9),
        lineterminator=chr(10),
        fieldnames=[
            "source_class",
            "sample_id",
            "safe_sample_id",
            "read_path",
            "action",
        ],
    )
    writer.writeheader()
    writer.writerows(prepared_rows)

with stats_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter=chr(9),
        lineterminator=chr(10),
        fieldnames=[
            "source_class",
            "sample_id",
            "safe_sample_id",
            "source_layout",
            "source_read1",
            "source_read2",
            "source_interleaved",
            "mvp_read_path",
            "action",
        ],
    )
    writer.writeheader()
    writer.writerows(stats_rows)

log("----------------------------------------")
log(f"Total read sets prepared for MVP: {len(prepared_rows)}")
log(f"Existing interleaved read sets: {sum(row['action'] == 'existing_interleaved' for row in prepared_rows)}")
log(f"Paired read sets interleaved: {sum(row['action'] == 'paired_to_interleaved' for row in prepared_rows)}")
PY

    prepared_sample_count="\$(
        tail -n +2 "\$READS_MANIFEST" |
        awk 'NF > 0' |
        wc -l |
        tr -d ' '
    )"

    printf 'MVP_READS_MANIFEST=%s\n' \
        "\$(pwd -P)/\$READS_MANIFEST" \
        > "\$STATUS"
    printf 'MVP_READ_SAMPLE_COUNT=%s\n' \
        "\$prepared_sample_count" \
        >> "\$STATUS"
    printf 'STATUS=complete\n' \
        >> "\$STATUS"

    echo "MVP read preparation finished: \$(date)" >> "\$LOG"
    """
}

/*
* Install and validate MVP.
*/
process SETUP_AUXMODULE2_MVP {
    tag "setup_auxmodule2_mvp"
    cache false

    publishDir "${params.outdir}/setup", mode: "copy", pattern: "auxmodule2_mvp_tools_status.env"

    output:
    path "auxmodule2_mvp_tools_status.env", emit: status

    script:
    def env_dir = params.tool_env_dir
        ? absPath(params.tool_env_dir)
        : absPath("${params.outdir}/conda_envs/mvip")

    """
    set -euo pipefail

    STATUS_FILE="auxmodule2_mvp_tools_status.env"
    TOOL_ENV="${env_dir}"
    REQUESTED_PACKAGE="mvip=${params.mvp_version}"

    echo "Auxiliary Module 2 MVP setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested Conda package: \$REQUESTED_PACKAGE" >> "\$STATUS_FILE"
    echo "Requested Conda environment: \$TOOL_ENV" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    configure_mvp_r_home() {
        local r_executable="\$TOOL_ENV/bin/R"

        if [[ ! -x "\$r_executable" ]]; then
            echo "ERROR: MVP requires R, but the Conda R executable was not found: \$r_executable" \
                >> "\$STATUS_FILE"
            return 1
        fi

        export PATH="\$TOOL_ENV/bin:\$PATH"
        export R_HOME="\$("\$r_executable" RHOME)"

        if [[ -z "\$R_HOME" || ! -d "\$R_HOME" ]]; then
            echo "ERROR: Could not determine a valid R_HOME from: \$r_executable" \
                >> "\$STATUS_FILE"
            return 1
        fi

        echo "R executable: \$r_executable" >> "\$STATUS_FILE"
        echo "R_HOME=\$R_HOME" >> "\$STATUS_FILE"
    }

    validate_mvip() {
        local executable="\$1"
        local python_executable="\$2"

        [[ -x "\$executable" ]] || return 1
        [[ -x "\$python_executable" ]] || return 1
        configure_mvp_r_home || return 1

        "\$python_executable" --version >> "\$STATUS_FILE" 2>&1 || return 1
        "\$executable" --help >> "\$STATUS_FILE" 2>&1 \
            || "\$executable" -h >> "\$STATUS_FILE" 2>&1 \
            || return 1

        return 0
    }

    get_module06_source_path() {
        find "\$TOOL_ENV/lib" \
            -type f \
            -path '*/site-packages/mvip/modules/MVP_06_do_functional_annotation.py' \
            -print \
            -quit
    }

    patch_mvip_module06_evalues() {
        local python_executable="\$1"
        local module06_file="\$2"

        if [[ ! -x "\$python_executable" ]]; then
            echo "ERROR: MVP Python executable was not found: \$python_executable" \
                >> "\$STATUS_FILE"
            return 1
        fi

        if [[ ! -f "\$module06_file" ]]; then
            echo "ERROR: MVP Module 06 source file was not found: \$module06_file" \
                >> "\$STATUS_FILE"
            return 1
        fi

        echo "Patching MVP Module 06 E-value argument types." \
            >> "\$STATUS_FILE"
        echo "Module 06 source: \$module06_file" \
            >> "\$STATUS_FILE"

        "\$python_executable" - "\$module06_file" \
            >> "\$STATUS_FILE" 2>&1 <<'PYTHON_SCRIPT'
import sys
from pathlib import Path

module_file = Path(sys.argv[1]).resolve()

if not module_file.is_file():
    raise RuntimeError(
        "MVP Module 06 source file does not exist: {}".format(module_file)
    )

targets = {
    "--PHROGS_evalue",
    "--PFAM_evalue",
    "--ADS_evalue",
    "--RdRP_evalue",
}

lines = module_file.read_text().splitlines(keepends=True)

patched = []
already_float = []
found = set()
current_argument = None

for index, line in enumerate(lines):
    stripped = line.strip()

    for argument in targets:
        if argument in stripped:
            current_argument = argument
            found.add(argument)
            break

    if current_argument is None:
        continue

    if "type=float" in stripped:
        already_float.append(current_argument)
        current_argument = None
        continue

    if "type=int" in stripped:
        lines[index] = line.replace("type=int", "type=float")
        patched.append(current_argument)
        current_argument = None

missing = targets - found
if missing:
    raise RuntimeError(
        "Could not find expected Module 06 arguments: {}".format(
            ", ".join(sorted(missing))
        )
    )

module_file.write_text("".join(lines))

print("Patched MVP Module 06 source: {}".format(module_file))

if patched:
    print(
        "Changed type=int to type=float for: {}".format(
            ", ".join(patched)
        )
    )

if already_float:
    print(
        "Already type=float for: {}".format(
            ", ".join(already_float)
        )
    )
PYTHON_SCRIPT

        echo "MVP Module 06 E-value patch completed successfully." \
            >> "\$STATUS_FILE"
    }

    MVP_EXECUTABLE="\$TOOL_ENV/bin/mvip"
    MVP_PYTHON="\$TOOL_ENV/bin/python"

    if [[ -d "\$TOOL_ENV" ]]; then
        echo "Existing MVP Conda environment detected: \$TOOL_ENV" \
            >> "\$STATUS_FILE"

        if validate_mvip "\$MVP_EXECUTABLE" "\$MVP_PYTHON"; then
            echo "Existing MVP environment passed validation." \
                >> "\$STATUS_FILE"

            MVP_MODULE06_SOURCE="\$(get_module06_source_path)"

            patch_mvip_module06_evalues \
                "\$MVP_PYTHON" \
                "\$MVP_MODULE06_SOURCE"

            echo "MVP_MODULE06_EVALUES=type_float_patch_applied" \
                >> "\$STATUS_FILE"
            echo "MVP_MODULE06_SOURCE=\$MVP_MODULE06_SOURCE" \
                >> "\$STATUS_FILE"
            echo "TOOL_ENV=\$TOOL_ENV" \
                >> "\$STATUS_FILE"
            echo "MVIP_EXECUTABLE=\$MVP_EXECUTABLE" \
                >> "\$STATUS_FILE"
            echo "PYTHON_EXECUTABLE=\$MVP_PYTHON" \
                >> "\$STATUS_FILE"
            echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
                >> "\$STATUS_FILE"

            exit 0
        fi

        echo "Existing MVP environment is incomplete or failed validation." \
            >> "\$STATUS_FILE"
        echo "Removing invalid MVP environment: \$TOOL_ENV" \
            >> "\$STATUS_FILE"

        rm -rf "\$TOOL_ENV"
    fi

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: MVP is not installed and auto-install is disabled." \
            >> "\$STATUS_FILE"
        echo "Expected executable: \$MVP_EXECUTABLE" \
            >> "\$STATUS_FILE"
        exit 1
    fi

    INSTALLER=""

    if command -v mamba >/dev/null 2>&1; then
        INSTALLER="mamba"
        echo "Using mamba: \$(command -v mamba)" \
            >> "\$STATUS_FILE"
    elif command -v conda >/dev/null 2>&1; then
        INSTALLER="conda"
        echo "Using conda: \$(command -v conda)" \
            >> "\$STATUS_FILE"
    else
        echo "ERROR: Neither mamba nor conda was found in PATH." \
            >> "\$STATUS_FILE"
        exit 1
    fi

    mkdir -p "\$(dirname "\$TOOL_ENV")"

    echo "Creating MVP Conda environment: \$TOOL_ENV" \
        >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \
        -p "\$TOOL_ENV" \
        -c conda-forge \
        -c bioconda \
        "\$REQUESTED_PACKAGE" \
        "r-base" \
        >> "\$STATUS_FILE" 2>&1

    if ! validate_mvip "\$MVP_EXECUTABLE" "\$MVP_PYTHON"; then
        echo "ERROR: MVP installation completed but validation failed." \
            >> "\$STATUS_FILE"
        echo "Expected executable: \$MVP_EXECUTABLE" \
            >> "\$STATUS_FILE"
        exit 1
    fi

    MVP_MODULE06_SOURCE="\$(get_module06_source_path)"
    
    patch_mvip_module06_evalues \
        "\$MVP_PYTHON" \
        "\$MVP_MODULE06_SOURCE"

    echo "MVP_MODULE06_EVALUES=type_float_patch_applied" \
        >> "\$STATUS_FILE"
    echo "MVP_MODULE06_SOURCE=\$MVP_MODULE06_SOURCE" \
        >> "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" \
        >> "\$STATUS_FILE"
    echo "MVIP_EXECUTABLE=\$MVP_EXECUTABLE" \
        >> "\$STATUS_FILE"
    echo "PYTHON_EXECUTABLE=\$MVP_PYTHON" \
        >> "\$STATUS_FILE"
    echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
        >> "\$STATUS_FILE"
    """
}

process PREPARE_MVP_METADATA {
    tag "prepare_mvp_metadata"

    publishDir "${params.outdir}/metadata", mode: "copy", pattern: "mvp_metadata.tsv"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "mvp_input_summary.tsv"

    publishDir "${params.outdir}/logs", mode: "copy", pattern: "prepare_mvp_metadata.log"

    input:
    val include_individual
    val include_coassemblies
    val include_subtractive
    val assembly_manifest
    path mvp_reads_manifest
    val coassembly_manifest
    val subtractive_manifest
    val subtractive_assembly_dir
    val base_dir
    path mvp_tools_status

    output:
    path "mvp_metadata.tsv", emit: metadata
    path "mvp_input_summary.tsv", emit: input_summary
    path "prepare_mvp_metadata_status.env", emit: status
    path "prepare_mvp_metadata.log", emit: log

    script:
    """
    set -euo pipefail

    LOG="prepare_mvp_metadata.log"
    METADATA="mvp_metadata.tsv"
    INPUT_SUMMARY="mvp_input_summary.tsv"
    STATUS="prepare_mvp_metadata_status.env"

    TOOL_ENV="\$(
        grep '^TOOL_ENV=' "${mvp_tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"

    PYTHON_EXECUTABLE="\$(
        grep '^PYTHON_EXECUTABLE=' "${mvp_tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true
    )"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    if [[ -z "\$PYTHON_EXECUTABLE" || ! -x "\$PYTHON_EXECUTABLE" ]]; then
        echo "ERROR: Python executable was not found in the MVP environment." >&2
        cat "${mvp_tools_status}" >&2 || true
        exit 1
    fi

    MVP_READS_MANIFEST_ABS="\$(pwd -P)/${mvp_reads_manifest}"

    echo "MVP metadata preparation started: \$(date)" > "\$LOG"
    echo "Include individual assemblies: ${include_individual}" >> "\$LOG"
    echo "Include coassemblies: ${include_coassemblies}" >> "\$LOG"
    echo "Include subtractive assemblies: ${include_subtractive}" >> "\$LOG"
    echo "Module 2 assembly manifest: ${assembly_manifest}" >> "\$LOG"
    echo "Prepared MVP reads manifest: \$MVP_READS_MANIFEST_ABS" >> "\$LOG"
    echo "Coassembly manifest: ${coassembly_manifest}" >> "\$LOG"
    echo "Subtractive manifest: ${subtractive_manifest}" >> "\$LOG"
    echo "Subtractive assembly directory: ${subtractive_assembly_dir}" >> "\$LOG"
    echo "Path resolution base directory: ${base_dir}" >> "\$LOG"
    echo "Task directory: \$(pwd -P)" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    "\$PYTHON_EXECUTABLE" - \\
        "${include_individual}" \\
        "${include_coassemblies}" \\
        "${include_subtractive}" \\
        "${assembly_manifest}" \\
        "\$MVP_READS_MANIFEST_ABS" \\
        "${coassembly_manifest}" \\
        "${subtractive_manifest}" \\
        "${subtractive_assembly_dir}" \\
        "${base_dir}" \\
        "\$METADATA" \\
        "\$INPUT_SUMMARY" \\
        "\$LOG" <<'PYTHON_SCRIPT'

import csv
import re
import sys
from pathlib import Path


(
    include_individual,
    include_coassemblies,
    include_subtractive,
    assembly_manifest,
    mvp_reads_manifest,
    coassembly_manifest,
    subtractive_manifest,
    subtractive_assembly_dir,
    base_dir,
    output_metadata,
    output_summary,
    log_file,
) = sys.argv[1:]


include_individual = include_individual.lower() == "true"
include_coassemblies = include_coassemblies.lower() == "true"
include_subtractive = include_subtractive.lower() == "true"

base_dir = Path(base_dir).resolve()
output_metadata = Path(output_metadata)
output_summary = Path(output_summary)
log_file = Path(log_file)


def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)


def resolve_path(value):
    value = str(value or "").strip()

    if not value:
        return ""

    path = Path(value)

    if not path.is_absolute():
        path = base_dir / path

    return str(path.resolve())


def read_tsv(path_text, required_columns=None, required=True):
    required_columns = required_columns or []

    path = Path(path_text)

    if not path.is_absolute():
        path = Path(resolve_path(path_text))

    if not path.is_file():
        if required:
            raise RuntimeError(f"Required file does not exist: {path}")
        return []

    if path.stat().st_size == 0:
        if required:
            raise RuntimeError(f"Required file is empty: {path}")
        return []

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\\t")

        observed_columns = set(reader.fieldnames or [])
        missing_columns = set(required_columns) - observed_columns

        if missing_columns:
            raise RuntimeError(
                f"File {path} is missing required columns: "
                f"{sorted(missing_columns)}"
            )

        rows = [
            row
            for row in reader
            if any(str(value or "").strip() for value in row.values())
        ]

    return rows


def file_exists_nonempty(path_text):
    path = Path(resolve_path(path_text))
    return path.is_file() and path.stat().st_size > 0


def first_existing_column(rows, candidates):
    if not rows:
        return None

    available_columns = set(rows[0].keys())

    for candidate in candidates:
        if candidate in available_columns:
            return candidate

    return None


def infer_path_column(rows, preferred_candidates, label):
    if not rows:
        return None

    available_columns = list(rows[0].keys())

    for candidate in preferred_candidates:
        if candidate in available_columns:
            return candidate

    path_pattern = re.compile(
        r"(assembly|contig|fasta|fa|path|renamed)",
        re.IGNORECASE,
    )

    fallback_candidates = [
        column
        for column in available_columns
        if path_pattern.search(column)
    ]

    for candidate in fallback_candidates:
        if any(str(row.get(candidate, "")).strip() for row in rows):
            log(
                f"{label}: inferred path column '{candidate}' "
                f"from observed columns {available_columns}"
            )
            return candidate

    return None


def require_nonempty_file(path_text, label):
    path = Path(path_text)

    if not path.is_absolute():
        path = Path(resolve_path(path_text))

    if not path.is_file():
        raise RuntimeError(f"{label} does not exist: {path}")

    if path.stat().st_size == 0:
        raise RuntimeError(f"{label} is empty: {path}")

    return str(path.resolve())


def safe_token(value):
    value = str(value or "").strip()

    if not value:
        return ""

    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("._-")

    return value


def build_sample_name(row, primary_candidates, extra_candidates):
    primary = ""

    for column in primary_candidates:
        value = safe_token(row.get(column, ""))

        if value:
            primary = value
            break

    if not primary:
        raise RuntimeError(
            f"Could not derive a sample name from row: {row}"
        )

    extras = []

    for column in extra_candidates:
        value = safe_token(row.get(column, ""))

        if value:
            extras.append(value)

    return "_".join([primary] + extras)


prepared_read_rows = read_tsv(
    mvp_reads_manifest,
    [
        "source_class",
        "sample_id",
        "safe_sample_id",
        "read_path",
        "action",
    ],
    required=True,
)


read_lookup = {}

for row in prepared_read_rows:
    source_class = str(row.get("source_class", "")).strip()
    sample_id = str(row.get("sample_id", "")).strip()
    safe_sample_id = str(row.get("safe_sample_id", "")).strip()

    read_path = require_nonempty_file(
        row.get("read_path", ""),
        f"Prepared MVP read file for "
        f"{source_class} sample {sample_id}",
    )

    action = str(row.get("action", "")).strip()

    if not source_class or not sample_id:
        raise RuntimeError(
            "Prepared MVP reads manifest contains a row with an empty "
            "source_class or sample_id"
        )

    read_record = {
        "read_path": read_path,
        "action": action,
        "sample_id": sample_id,
        "safe_sample_id": safe_sample_id,
    }

    read_lookup[(source_class, sample_id)] = read_record

    if safe_sample_id:
        read_lookup[(source_class, safe_sample_id)] = read_record


def source_has_prepared_reads(source_class):
    return any(
        key[0] == source_class
        for key in read_lookup.keys()
    )


records = []
summary_rows = []


def add_record(
    sample_name,
    assembly_path,
    read_path,
    source,
    action="",
):
    sample_name = str(sample_name or "").strip()

    if not sample_name:
        raise RuntimeError(
            f"Cannot add {source} record with an empty sample name"
        )

    assembly_path = require_nonempty_file(
        assembly_path,
        f"Assembly for {source} sample {sample_name}",
    )

    read_path = require_nonempty_file(
        read_path,
        f"Read file for {source} sample {sample_name}",
    )

    if any(
        existing["Sample"] == sample_name
        for existing in records
    ):
        raise RuntimeError(
            f"Duplicate MVP sample name: {sample_name}"
        )

    record = {
        "Sample_number": len(records) + 1,
        "Sample": sample_name,
        "Assembly_Path": assembly_path,
        "Read_Path": read_path,
    }

    records.append(record)

    summary_rows.append(
        {
            "Sample_number": record["Sample_number"],
            "Sample": sample_name,
            "Source": source,
            "Action": action,
            "Assembly_Path": assembly_path,
            "Read_Path": read_path,
        }
    )

    log(
        f"Added MVP sample {record['Sample_number']}: "
        f"{sample_name} [{source}]"
    )


if include_individual:
    log("Preparing individual Module 2 assemblies")

    assembly_rows = read_tsv(
        assembly_manifest,
        ["sample_id"],
        required=True,
    )

    observed_columns = (
        list(assembly_rows[0].keys())
        if assembly_rows
        else []
    )

    assembly_path_column = infer_path_column(
        assembly_rows,
        [
            "assembly_path",
            "assembly",
            "assembly_file",
            "final_assembly",
            "contigs",
            "contig_path",
            "assembly_fasta",
            "fasta_path",
            "fasta",
            "contigs_fasta",
            "scaffolds",
            "scaffold_path",
            "renamed_fasta",
        ],
        "Module 2 assembly manifest",
    )

    if not assembly_path_column:
        raise RuntimeError(
            "Module 2 assembly manifest does not contain a "
            "recognizable assembly path column. Observed columns: "
            + ", ".join(observed_columns)
        )

    safe_sample_column = first_existing_column(
        assembly_rows,
        ["safe_sample_id"],
    )

    for row in assembly_rows:
        sample_id = str(row.get("sample_id", "")).strip()

        safe_sample_id = ""

        if safe_sample_column:
            safe_sample_id = str(
                row.get(safe_sample_column, "")
            ).strip()

        if not sample_id:
            raise RuntimeError(
                "Encountered Module 2 assembly row with an "
                "empty sample_id"
            )

        read_record = None

        for candidate in (sample_id, safe_sample_id):
            if not candidate:
                continue

            key = ("individual", candidate)

            if key in read_lookup:
                read_record = read_lookup[key]
                break

        if read_record is None:
            message = (
                "No prepared MVP read set could be matched to "
                f"Module 2 assembly sample '{sample_id}'. "
                "Tried source class 'individual' with "
                f"sample_id '{sample_id}'"
            )

            if safe_sample_id:
                message += (
                    f" and safe_sample_id '{safe_sample_id}'."
                )
            else:
                message += "."

            raise RuntimeError(message)

        sample_name = build_sample_name(
            row,
            [
                "assembly_sample_id",
                "safe_sample_id",
                "sample_id",
            ],
            [
                "assembler",
                "assembly_mode",
                "rarefaction_label",
                "assembly_strategy",
            ],
        )

        assembly_path = row.get(
            assembly_path_column,
            "",
        )

        add_record(
            sample_name=sample_name,
            assembly_path=assembly_path,
            read_path=read_record["read_path"],
            source="individual",
            action=read_record.get("action", ""),
        )


if include_coassemblies:
    log("Preparing coassembly entries")

    if not file_exists_nonempty(coassembly_manifest):
        log(
            "Module 2b assembly manifest was not found or was empty. "
            "Proceeding without coassemblies."
        )

    elif not source_has_prepared_reads("coassembly"):
        log(
            "No prepared coassembly read sets were found. "
            "Proceeding without coassemblies."
        )

    else:
        coassembly_rows = read_tsv(
            coassembly_manifest,
            [],
            required=True,
        )

        observed_columns = (
            list(coassembly_rows[0].keys())
            if coassembly_rows
            else []
        )

        coassembly_id_column = first_existing_column(
            coassembly_rows,
            [
                "sample_id",
                "coassembly_id",
                "coassembly",
                "sample",
            ],
        )

        coassembly_assembly_column = infer_path_column(
            coassembly_rows,
            [
                "assembly_path",
                "assembly",
                "assembly_file",
                "final_assembly",
                "contigs",
                "contig_path",
                "assembly_fasta",
                "fasta_path",
                "fasta",
                "contigs_fasta",
                "scaffolds",
                "scaffold_path",
                "renamed_fasta",
            ],
            "Module 2b coassembly manifest",
        )

        if not coassembly_id_column:
            raise RuntimeError(
                "Could not identify the coassembly sample ID column. "
                "Observed columns: "
                + ", ".join(observed_columns)
            )

        if not coassembly_assembly_column:
            raise RuntimeError(
                "Could not identify the coassembly assembly path column. "
                "Observed columns: "
                + ", ".join(observed_columns)
            )

        for row in coassembly_rows:
            raw_id = str(
                row.get(coassembly_id_column, "")
            ).strip()

            if not raw_id:
                raise RuntimeError(
                    "Encountered coassembly assembly with "
                    "an empty sample ID"
                )

            key = ("coassembly", raw_id)

            if key not in read_lookup:
                raise RuntimeError(
                    "No prepared MVP read path was found for "
                    f"coassembly {raw_id}"
                )

            sample_name = build_sample_name(
                row,
                [
                    coassembly_id_column,
                    "sample_id",
                    "coassembly_id",
                    "coassembly",
                    "sample",
                ],
                [
                    "assembler",
                    "assembly_mode",
                    "rarefaction_label",
                    "assembly_strategy",
                ],
            )

            add_record(
                sample_name=sample_name,
                assembly_path=row.get(
                    coassembly_assembly_column,
                    "",
                ),
                read_path=read_lookup[key]["read_path"],
                source="coassembly",
                action=read_lookup[key].get("action", ""),
            )


if include_subtractive:
    log("Preparing subtractive assembly entries")

    if not file_exists_nonempty(subtractive_manifest):
        log(
            "Module 5 subtractive assembly manifest was not found "
            "or was empty. Proceeding without subtractive assemblies."
        )

    elif not source_has_prepared_reads("subtractive"):
        log(
            "No prepared subtractive read sets were found. "
            "Proceeding without subtractive assemblies."
        )

    else:
        subtractive_rows = read_tsv(
            subtractive_manifest,
            [],
            required=True,
        )

        subtractive_sample_column = first_existing_column(
            subtractive_rows,
            [
                "sample_id",
                "safe_sample_id",
                "sample",
            ],
        )

        subtractive_assembly_column = infer_path_column(
            subtractive_rows,
            [
                "assembly_path",
                "assembly",
                "assembly_file",
                "final_assembly",
                "contigs",
                "contig_path",
                "assembly_fasta",
                "fasta_path",
                "fasta",
                "renamed_fasta",
            ],
            "Module 5 subtractive manifest",
        )

        if not subtractive_sample_column:
            raise RuntimeError(
                "Could not identify the subtractive sample ID column"
            )

        subtractive_assembly_root = Path(
            resolve_path(subtractive_assembly_dir)
        )

        for row in subtractive_rows:
            raw_id = str(
                row.get(subtractive_sample_column, "")
            ).strip()

            if not raw_id:
                raise RuntimeError(
                    "Encountered subtractive assembly row with "
                    "an empty sample ID"
                )

            assembly_path = ""

            if subtractive_assembly_column:
                assembly_path = row.get(
                    subtractive_assembly_column,
                    "",
                )

            if not str(assembly_path or "").strip():
                assembly_path = str(
                    (
                        subtractive_assembly_root
                        / f"{raw_id}.fa"
                    ).resolve()
                )

            key = ("subtractive", raw_id)

            if key not in read_lookup:
                raise RuntimeError(
                    "No prepared MVP read path was found for "
                    f"subtractive sample {raw_id}"
                )

            sample_name = build_sample_name(
                row,
                [
                    subtractive_sample_column,
                    "safe_sample_id",
                    "sample_id",
                    "sample",
                ],
                [
                    "assembler",
                    "assembly_mode",
                    "rarefaction_label",
                    "assembly_strategy",
                ],
            )

            try:
                resolved_assembly = require_nonempty_file(
                    assembly_path,
                    f"Subtractive assembly for sample {raw_id}",
                )
            except RuntimeError as exc:
                log(
                    "Skipping subtractive assembly with no usable "
                    f"contigs: {sample_name}; {exc}"
                )
                continue

            add_record(
                sample_name=sample_name,
                assembly_path=resolved_assembly,
                read_path=read_lookup[key]["read_path"],
                source="subtractive",
                action=read_lookup[key].get("action", ""),
            )


if not records:
    raise RuntimeError(
        "No MVP input records were generated. At least one of "
        "individual assemblies, coassemblies, or subtractive "
        "assemblies must produce valid input."
    )


metadata_fields = [
    "Sample_number",
    "Sample",
    "Assembly_Path",
    "Read_Path",
]


with output_metadata.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=metadata_fields,
    )

    writer.writeheader()
    writer.writerows(records)


summary_fields = [
    "Sample_number",
    "Sample",
    "Source",
    "Action",
    "Assembly_Path",
    "Read_Path",
]


with output_summary.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=summary_fields,
    )

    writer.writeheader()
    writer.writerows(summary_rows)


log("----------------------------------------")
log(f"Total MVP input records: {len(records)}")
log(
    "Individual records: "
    f"{sum(row['Source'] == 'individual' for row in summary_rows)}"
)
log(
    "Coassembly records: "
    f"{sum(row['Source'] == 'coassembly' for row in summary_rows)}"
)
log(
    "Subtractive records: "
    f"{sum(row['Source'] == 'subtractive' for row in summary_rows)}"
)
log("Metadata preparation completed successfully")

PYTHON_SCRIPT

    record_count="\$(
        tail -n +2 "\$METADATA" |
        awk 'NF > 0' |
        wc -l |
        tr -d ' '
    )"

    {
        echo "MVP_METADATA=\$(pwd -P)/\$METADATA"
        echo "MVP_INPUT_SUMMARY=\$(pwd -P)/\$INPUT_SUMMARY"
        echo "MVP_RECORD_COUNT=\$record_count"
        echo "STATUS=complete"
    } > "\$STATUS"

    echo "MVP metadata preparation finished: \$(date)" >> "\$LOG"
    """
}

process RUN_MVP {
    tag "run_auxmodule2_mvp"
    cache false

    stageInMode "symlink"

    publishDir "${params.outdir}/logs", mode: "copy", pattern: "mvp_commands.log"

    publishDir "${params.outdir}/logs", mode: "copy", pattern: "mvp_run.log"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "mvp_complete.txt"

    cpus {
        params.threads as int
    }

    memory {
        def gb = params.memory_gb as int
        return gb > 0 ? "${gb} GB" : null
    }

    input:
    path metadata
    path input_summary
    path tools_status

    output:
    path "mvp_commands.log", emit: commands
    path "mvp_run.log", emit: log_file
    path "mvp_complete.txt", emit: completion_marker

    script:
    def selected_modules = parseMvpModules(params.mvp_modules)
    def canonical_order = ["0", "1", "2", "3", "4", "5", "6", "7", "99", "100"]

    def ordered_modules = canonical_order.findAll { module ->
        module in selected_modules
    }

    def run_module_0 = "0" in selected_modules
    def run_module_1 = "1" in selected_modules
    def run_module_2 = "2" in selected_modules
    def run_module_3 = "3" in selected_modules
    def run_module_4 = "4" in selected_modules
    def run_module_5 = "5" in selected_modules
    def run_module_6 = "6" in selected_modules
    def run_module_7 = "7" in selected_modules
    def run_module_99 = "99" in selected_modules
    def run_module_100 = "100" in selected_modules

    def install_databases = params.install_databases
        .toString()
        .toBoolean()

    def force = params.force
        .toString()
        .toBoolean()

    def mvp_work = absPath(params.outdir)

    /*
     * MVP's default database directory.
     */
    def database_root = params.mvp_database_dir
        ? absPath(params.mvp_database_dir)
        : "${mvp_work}/00_DATABASES"

    /*
     * geNomad path remains explicitly configurable.
     */
    def genomad_database = params.genomad_db_path
        ? absPath(params.genomad_db_path)
        : "${database_root}/genomad_db"

    /*
     * An explicitly supplied CheckV path takes precedence.
     * Otherwise, the shell function below detects checkv-db-v*.
     */
    def checkv_database = params.checkv_db_path
        ? absPath(params.checkv_db_path)
        : ""

    def module00_database_flag = install_databases
        ? ""
        : "--skip_install_databases"

    def skip_check_errors_argument = params.skip_check_errors.toString().toBoolean()
        ? "--skip_check_errors"
        : ""

    def min_seq_argument = (params.min_seq_size as int) > 0
        ? "--min_seq_size ${params.min_seq_size}"
        : ""

    def genomad_filter_argument = ""

    if (params.genomad_relaxed.toString().toBoolean()) {
        genomad_filter_argument = "--genomad_relaxed"
    }
    else if (params.genomad_conservative.toString().toBoolean()) {
        genomad_filter_argument = "--genomad_conservative"
    }

    def skip_modify_headers_argument = params.skip_modify_headers.toString().toBoolean()
        ? "--skip_modify_headers"
        : ""

    def unfiltered_protein_argument = params.unfiltered_protein_file.toString().toBoolean()
        ? "--Unfiltered_protein_file"
        : ""

    def force_module01_argument = force
        ? "--force_genomad --force_checkv"
        : ""

    def interleaved_argument = "--interleaved"

    def mapping_delete_argument = params.delete_mapping_intermediates.toString().toBoolean()
        ? "--delete_files"
        : ""

    def force_mapping_argument = force
        ? "--force_read_mapping"
        : ""

    def covered_fraction_argument = params.covered_fraction != null
        ? "--covered_fraction ${params.covered_fraction}"
        : ""

    def functional_ads_argument = params.functional_ads.toString().toBoolean()
        ? "--ADS"
        : ""

    def functional_rdrp_argument = params.functional_rdrp.toString().toBoolean()
        ? "--RdRP"
        : ""

    def functional_dram_argument = params.functional_dram.toString().toBoolean()
        ? "--DRAM"
        : ""

    def functional_delete_argument = params.delete_functional_intermediates.toString().toBoolean()
        ? "--delete_files"
        : ""

    def force_functional_argument = force
        ? "--force_prodigal --force_PHROGS --force_PFAM --force_outputs"
        : ""

    if (force && params.functional_ads.toString().toBoolean()) {
        force_functional_argument += " --force_ADS"
    }

    if (force && params.functional_rdrp.toString().toBoolean()) {
        force_functional_argument += " --force_RdRP"
    }

    def keep_bam_argument = params.keep_bam.toString().toBoolean()
        ? "--keep_bam"
        : ""

    def binning_delete_argument = params.delete_binning_intermediates.toString().toBoolean()
        ? "--delete_files"
        : ""

    def binning_sample_group_argument = params.binning_sample_group != null
        ? "--binning_sample_group ${params.binning_sample_group}"
        : ""

    def read_mapping_sample_group_argument = params.read_mapping_sample_group != null
        ? "--read_mapping_sample_group ${params.read_mapping_sample_group}"
        : ""

    def force_binning_argument = force
        ? "--force_vrhyme --force_checkv --force_read_mapping --force_outputs"
        : ""

    def summary_force_argument = force
        ? "--force"
        : ""

    def miuvig_template_argument = params.miuvig_template
        ? "-t ${absPath(params.miuvig_template)}"
        : ""

    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true)"

    PYTHON_EXECUTABLE="\$(grep '^PYTHON_EXECUTABLE=' "${tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    if ! command -v mvip >/dev/null 2>&1; then
        echo "ERROR: mvip was not found after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    if [[ -z "\$PYTHON_EXECUTABLE" || ! -x "\$PYTHON_EXECUTABLE" ]]; then
        echo "ERROR: Python executable was not found after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    METADATA="\$(
        "\$PYTHON_EXECUTABLE" -c \
        'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' \
        "${metadata}"
    )"

    MVP_WORK="${mvp_work}"
    DATABASE_ROOT="${database_root}"
    GENOMAD_DB="${genomad_database}"
    CHECKV_DB="${checkv_database}"

    mkdir -p "\$MVP_WORK"
    mkdir -p "\$DATABASE_ROOT"

    : > mvp_commands.log
    : > mvp_run.log

    echo "Auxiliary Module 2 MVP started: \$(date)" >> mvp_run.log
    echo "MVP executable: \$(command -v mvip)" >> mvp_run.log
    echo "MVP working directory: \$MVP_WORK" >> mvp_run.log
    echo "MVP metadata: \$METADATA" >> mvp_run.log
    echo "Database root: \$DATABASE_ROOT" >> mvp_run.log
    echo "Selected MVP modules: ${ordered_modules.join(',')}" >> mvp_run.log
    echo "Threads: ${task.cpus}" >> mvp_run.log
    echo "----------------------------------------" >> mvp_run.log

    run_command() {
        printf '%q ' "\$@" >> mvp_commands.log
        printf '\\n' >> mvp_commands.log

        echo "==================================================" >> mvp_run.log
        echo "Started: \$(date)" >> mvp_run.log
        printf 'Command: ' >> mvp_run.log
        printf '%q ' "\$@" >> mvp_run.log
        printf '\\n' >> mvp_run.log

        set +e
        "\$@" >> mvp_run.log 2>&1
        rc=\$?
        set -e

        if [[ "\$rc" -ne 0 ]]; then
            echo "Command failed with exit code \$rc: \$*" >> mvp_run.log
            echo "Last 100 lines of mvp_run.log:" >&2
            tail -n 100 mvp_run.log >&2
            exit "\$rc"
        fi

        echo "Finished: \$(date)" >> mvp_run.log
    }

    resolve_checkv_database() {
        if [[ -n "\$CHECKV_DB" ]]; then
            echo "Using explicitly supplied CheckV database:" \
                "\$CHECKV_DB" >> mvp_run.log
        else
            echo "Searching for CheckV database under:" \
                "\$DATABASE_ROOT" >> mvp_run.log

            mapfile -t CHECKV_CANDIDATES < <(
                find "\$DATABASE_ROOT" \
                    -mindepth 1 \
                    -maxdepth 1 \
                    -type d \
                    -name 'checkv-db-v*' \
                    -print |
                sort -V
            )

            if [[ "\${#CHECKV_CANDIDATES[@]}" -eq 0 ]]; then
                echo "ERROR: No CheckV database was found under:" >&2
                echo "       \$DATABASE_ROOT" >&2
                echo "Expected a directory such as:" >&2
                echo "       \$DATABASE_ROOT/checkv-db-v1.5" >&2
                exit 1
            fi

            CHECKV_DB="\${CHECKV_CANDIDATES[\${#CHECKV_CANDIDATES[@]}-1]}"

            echo "Detected CheckV database:" \
                "\$CHECKV_DB" >> mvp_run.log
        fi

        if [[ ! -d "\$CHECKV_DB" ]]; then
            echo "ERROR: CheckV database directory does not exist:" >&2
            echo "       \$CHECKV_DB" >&2
            exit 1
        fi

        if [[ -z "\$(find "\$CHECKV_DB" -type f -print -quit)" ]]; then
            echo "ERROR: CheckV database directory contains no files:" >&2
            echo "       \$CHECKV_DB" >&2
            exit 1
        fi
    }

    run_module_00() {
        local -a module00_db_args=()

        if [[ -n "\$GENOMAD_DB" ]]; then
            module00_db_args+=(--genomad_db_path "\$GENOMAD_DB")
        fi

        if [[ -n "\$CHECKV_DB" ]]; then
            module00_db_args+=(--checkv_db_path "\$CHECKV_DB")
        fi

        run_command \
            mvip MVP_00_set_up_MVP \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${module00_database_flag} \
            ${skip_check_errors_argument} \
            "\${module00_db_args[@]}"
    }

    if [[ "${run_module_0}" == "true" &&
          "${install_databases}" == "false" &&
          -z "\$CHECKV_DB" ]]; then
        resolve_checkv_database
    fi

    if [[ "${run_module_0}" == "true" ]]; then
        run_module_00
    fi

    if [[ "${run_module_1}" == "true" ||
          "${run_module_2}" == "true" ||
          "${run_module_3}" == "true" ||
          "${run_module_4}" == "true" ||
          "${run_module_5}" == "true" ||
          "${run_module_6}" == "true" ||
          "${run_module_7}" == "true" ]]; then
        resolve_checkv_database
    fi

    echo "geNomad database: \$GENOMAD_DB" >> mvp_run.log
    echo "CheckV database: \$CHECKV_DB" >> mvp_run.log

    if [[ "${run_module_1}" == "true" ]]; then
        run_command \
            mvip MVP_01_run_genomad_checkv \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --genomad_db_path "\$GENOMAD_DB"\
            --checkv_db_path "\$CHECKV_DB" \
            ${skip_modify_headers_argument} \
            ${min_seq_argument} \
            ${genomad_filter_argument} \
            ${force_module01_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_2}" == "true" ]]; then
        run_command \
            mvip MVP_02_filter_genomad_checkv \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --viral_min_genes ${params.viral_min_genes} \
            --host_viral_genes_ratio ${params.host_viral_genes_ratio}
    fi

    if [[ "${run_module_3}" == "true" ]]; then
        run_command \
            mvip MVP_03_do_clustering \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --min_ani ${params.min_ani} \
            --min_tcov ${params.min_tcov} \
            --min_qcov ${params.min_qcov} \
            --read-type ${params.read_type} \
            ${unfiltered_protein_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_4}" == "true" ]]; then
        run_command \
            mvip MVP_04_do_read_mapping \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --read-type ${params.read_type} \
            ${interleaved_argument} \
            ${mapping_delete_argument} \
            ${force_mapping_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_5}" == "true" ]]; then
        run_command \
            mvip MVP_05_create_vOTU_table \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${covered_fraction_argument} \
            --normalization ${params.normalization} \
            --filtration ${params.filtration} \
            --viral_min_genes ${params.viral_min_genes} \
            --host_viral_genes_ratio ${params.host_viral_genes_ratio}
    fi

    if [[ "${run_module_6}" == "true" ]]; then
        run_command \
            mvip MVP_06_do_functional_annotation \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --fasta_files ${params.functional_fasta_files} \
            --PHROGS_evalue ${params.phrogs_evalue} \
            --PHROGS_score ${params.phrogs_score} \
            --PFAM_evalue ${params.pfam_evalue} \
            --PFAM_score ${params.pfam_score} \
            ${functional_ads_argument} \
            --ADS_evalue ${params.ads_evalue} \
            --ADS_score ${params.ads_score} \
            ${functional_rdrp_argument} \
            --RdRP_evalue ${params.rdrp_evalue} \
            --RdRP_score ${params.rdrp_score} \
            ${functional_dram_argument} \
            ${functional_delete_argument} \
            ${force_functional_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_7}" == "true" ]]; then
        run_command \
            mvip MVP_07_do_binning \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${binning_sample_group_argument} \
            ${read_mapping_sample_group_argument} \
            ${keep_bam_argument} \
            ${binning_delete_argument} \
            ${interleaved_argument} \
            ${force_binning_argument} \
            --read-type ${params.read_type} \
            --filtration ${params.filtration} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_99}" == "true" ]]; then
        run_command \
            mvip MVP_99_prep_MIUViG_submission \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            -g "${params.miuvig_identifier}" \
            -s "${params.miuvig_step}" \
            ${miuvig_template_argument}
    fi

    if [[ "${run_module_100}" == "true" ]]; then
        run_command \
            mvip MVP_100_summarize_outputs \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${summary_force_argument}
    fi

    echo "==================================================" >> mvp_run.log
    echo "Auxiliary Module 2 MVP finished: \$(date)" >> mvp_run.log

    cat > mvp_complete.txt <<EOF
status=complete
finished=\$(date --iso-8601=seconds 2>/dev/null || date)
mvp_work=\$MVP_WORK
metadata=\$METADATA
checkv_database=\$CHECKV_DB
selected_modules=${ordered_modules.join(',')}
EOF
    """
}
