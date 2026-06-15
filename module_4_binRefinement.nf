#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Module 4: Bin refinement / MAGScoT preparation and execution.
*/

/*
 * Parameters
 */

params.working_dir = null
params.output_dir = null
params.input_binning_manifest = null
params.dependencies_dir = "${projectDir}/dependencies"
params.tigrfam_hmm = null
params.pfam_hmm = null
params.magscot_script = null
params.magscot_profiles_dir = null
params.auto_install = true
params.tool_env_dir = null
params.threads = null
params.hmm_threads = 8
params.r_base_version = null
params.hmmer_version = null
params.prodigal_version = null
params.parallel_version = null
params.magscot_extra_args = ""
params.magscot_threshold = 0
params.publish_gathered_bins_mode = "copy"
params.publish_refined_bins_mode = "copy"
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.module3_outdir = "${params.results_dir}/module_3_binning"
params.outdir = "${params.results_dir}/module_4_binrefinement"

def absOrEmpty(value) {
    def s = value == null ? "" : value.toString().trim()
    if (!s || s == "null" || s == "NA") {
        return ""
    }
    return java.nio.file.Paths.get(s).toAbsolutePath().normalize().toString()
}

def firstExistingPath(List paths) {
        def found = paths.find { p ->
        java.nio.file.Files.exists(java.nio.file.Paths.get(p.toString()))
    }

    def selected = found ?: paths[0]

    return java.nio.file.Paths
        .get(selected.toString())
        .toAbsolutePath()
        .normalize()
        .toString()
}

workflow {

    def binning_manifest_file = params.input_binning_manifest ?: "${params.module3_outdir}/summary/binning_manifest.tsv"

    def tigrfam_hmm_file = params.tigrfam_hmm ?: firstExistingPath(
        ["${params.dependencies_dir}/hmm/gtdbtk_rel207_tigrfam.hmm", "${params.dependencies_dir}/gtdbtk_rel207_tigrfam.hmm", "${projectDir}/dependencies/hmm/gtdbtk_rel207_tigrfam.hmm", "${projectDir}/dependencies/gtdbtk_rel207_tigrfam.hmm"]
    )

    def pfam_hmm_file = params.pfam_hmm ?: firstExistingPath(
        ["${params.dependencies_dir}/hmm/gtdbtk_rel207_Pfam-A.hmm", "${params.dependencies_dir}/gtdbtk_rel207_Pfam-A.hmm", "${projectDir}/dependencies/hmm/gtdbtk_rel207_Pfam-A.hmm", "${projectDir}/dependencies/gtdbtk_rel207_Pfam-A.hmm", "${projectDir}/gtdbtk_rel207_Pfam-A.hmm"]
    )

    def magscot_script_file = params.magscot_script ?: firstExistingPath(
        ["${params.dependencies_dir}/MAGScoT.py", "${projectDir}/dependencies/MAGScoT.py", "${projectDir}/MAGScoT.py"]
    )

    def magscot_profiles_dir = params.magscot_profiles_dir ?: firstExistingPath(
        ["${params.dependencies_dir}/", "${projectDir}/dependencies/", "${params.dependencies_dir}/MAGScoT_profiles", "${projectDir}/dependencies/MAGScoT_profiles"]
    )

    if (params.magscot_threshold != null) {
        def threshold_value = params.magscot_threshold as double

        if (threshold_value < 0 || threshold_value > 1) {
            error(
                """
        Invalid --magscot_threshold value: ${params.magscot_threshold}

        Expected a value between 0 and 1.
        Example:
          --magscot_threshold 0.5
        """.stripIndent()
            )
        }
    }

    log.info("Module 4 results directory: ${params.results_dir}")
    log.info("Using Module 3 binning manifest: ${binning_manifest_file}")
    log.info("Writing Module 4 outputs to: ${params.outdir}")
    log.info("Dependencies directory: ${params.dependencies_dir}")
    log.info("TIGRFAM HMM: ${tigrfam_hmm_file}")
    log.info("Pfam HMM: ${pfam_hmm_file}")
    log.info("MAGScoT script: ${magscot_script_file}")
    log.info("HMMER threads: ${params.threads ?: params.hmm_threads}")
    log.info("MAGScoT profiles directory: ${magscot_profiles_dir}")


    def binning_manifest_ch = channel.fromPath(
        binning_manifest_file,
        type: 'file',
        checkIfExists: true,
    )

    def tigrfam_hmm_ch = channel.fromPath(
        tigrfam_hmm_file,
        type: 'file',
        checkIfExists: true,
    )

    def pfam_hmm_ch = channel.fromPath(
        pfam_hmm_file,
        type: 'file',
        checkIfExists: true,
    )

    def magscot_script_ch = channel.fromPath(
        magscot_script_file,
        type: 'file',
        checkIfExists: true,
    )

    def magscot_profiles_ch = channel.fromPath(
        magscot_profiles_dir,
        type: 'dir',
        checkIfExists: true,
    )

    SETUP_MODULE4_TOOLS()

    PREPARE_MAG_COLLECTION(
        binning_manifest_ch,
        SETUP_MODULE4_TOOLS.out.status,
    )

    RUN_PRODIGAL_ON_MAG_CONTIGS(
        PREPARE_MAG_COLLECTION.out.concatenated_fasta,
        SETUP_MODULE4_TOOLS.out.status,
    )

    RUN_HMMSEARCH_MAG_CONTIGS(
        RUN_PRODIGAL_ON_MAG_CONTIGS.out.proteins,
        tigrfam_hmm_ch,
        pfam_hmm_ch,
        SETUP_MODULE4_TOOLS.out.status,
    )

    RUN_MAGSCOT(
        PREPARE_MAG_COLLECTION.out.contigs_to_bin,
        RUN_HMMSEARCH_MAG_CONTIGS.out.hmm_table,
        magscot_script_ch,
        magscot_profiles_ch,
        SETUP_MODULE4_TOOLS.out.status,
    )

    BUILD_REFINED_MAGS(
        PREPARE_MAG_COLLECTION.out.concatenated_fasta,
        RUN_MAGSCOT.out.magscot_outputs,
    )

    WRITE_MODULE4_SUMMARY(
        PREPARE_MAG_COLLECTION.out.collection_stats,
        RUN_PRODIGAL_ON_MAG_CONTIGS.out.prodigal_status,
        RUN_HMMSEARCH_MAG_CONTIGS.out.hmm_status,
        RUN_MAGSCOT.out.magscot_status,
        BUILD_REFINED_MAGS.out.refined_stats,
    )
}

process SETUP_MODULE4_TOOLS {
    tag "setup_module4_tools"

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "module4_tools_status.env"

    output:
    path "module4_tools_status.env", emit: status

    script:
    def env_dir = params.tool_env_dir ?: "${params.outdir}/conda_envs/module4_tools"

    def packages = []
    packages << "python"
    packages << "pandas"
    packages << "r-base"
    packages << "r-optparse"
    packages << "r-dplyr"
    packages << "r-readr"
    packages << "r-funr"
    packages << "r-digest"
    packages << "hmmer"
    packages << "prodigal"
    packages << "parallel"

    if (params.r_base_version) {
        packages = packages.collect { pkg ->
            pkg == "r-base" ? "r-base=${params.r_base_version}" : pkg
        }
    }

    if (params.hmmer_version) {
        packages = packages.collect { pkg ->
            pkg == "hmmer" ? "hmmer=${params.hmmer_version}" : pkg
        }
    }

    if (params.prodigal_version) {
        packages = packages.collect { pkg ->
            pkg == "prodigal" ? "prodigal=${params.prodigal_version}" : pkg
        }
    }

    if (params.parallel_version) {
        packages = packages.collect { pkg ->
            pkg == "parallel" ? "parallel=${params.parallel_version}" : pkg
        }
    }

    def package_string = packages
        .collect { pkg ->
            "\"${pkg}\""
        }
        .join(" \\\n        ")

    """
    set -euo pipefail

    STATUS_FILE="module4_tools_status.env"
    TOOL_ENV="${env_dir}"

    echo "Module 4 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Requested packages:" >> "\$STATUS_FILE"
    echo "${packages.join(' ')}" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    find_installer() {
        if command -v mamba >/dev/null 2>&1; then
            echo "mamba"
        elif command -v conda >/dev/null 2>&1; then
            echo "conda"
        else
            echo ""
        fi
    }

    check_env() {
        local prefix="\$1"
        local ok="true"

        echo "Checking Module 4 environment: \$prefix" >> "\$STATUS_FILE"

        if [[ ! -x "\$prefix/bin/Rscript" ]]; then
            echo "Missing Rscript" >> "\$STATUS_FILE"
            ok="false"
        else
            "\$prefix/bin/Rscript" --version >> "\$STATUS_FILE" 2>&1 || true
        fi

        if [[ ! -x "\$prefix/bin/hmmsearch" ]]; then
            echo "Missing hmmsearch" >> "\$STATUS_FILE"
            ok="false"
        else
            "\$prefix/bin/hmmsearch" -h >> "\$STATUS_FILE" 2>&1 || true
        fi

        if [[ ! -x "\$prefix/bin/prodigal" ]]; then
            echo "Missing prodigal" >> "\$STATUS_FILE"
            ok="false"
        else
            "\$prefix/bin/prodigal" -v >> "\$STATUS_FILE" 2>&1 || true
        fi

        if [[ ! -x "\$prefix/bin/parallel" ]]; then
            echo "Missing GNU parallel" >> "\$STATUS_FILE"
            ok="false"
        else
            "\$prefix/bin/parallel" --version >> "\$STATUS_FILE" 2>&1 || true
        fi

        if [[ ! -x "\$prefix/bin/python" ]]; then
            echo "Missing python" >> "\$STATUS_FILE"
            ok="false"
        else
            "\$prefix/bin/python" --version >> "\$STATUS_FILE" 2>&1 || true
        fi

        echo "Checking required R packages..." >> "\$STATUS_FILE"
        "\$prefix/bin/Rscript" -e 'library(optparse); library(dplyr); library(readr); library(funr); library(digest)' >> "\$STATUS_FILE" 2>&1 || ok="false"

        echo "Checking required Python packages..." >> "\$STATUS_FILE"
        "\$prefix/bin/python" -c 'import pandas; print("pandas", pandas.__version__)' >> "\$STATUS_FILE" 2>&1 || ok="false"

        [[ "\$ok" == "true" ]]
    }

    if [[ -d "\$TOOL_ENV" ]]; then
        if check_env "\$TOOL_ENV"; then
            echo "Existing Module 4 environment passed checks." >> "\$STATUS_FILE"
            echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
            echo "Module 4 tool setup finished: \$(date)" >> "\$STATUS_FILE"
            exit 0
        else
            echo "Existing Module 4 environment failed checks. Removing." >> "\$STATUS_FILE"
            rm -rf "\$TOOL_ENV"
        fi
    fi

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: Module 4 environment missing/broken and --auto_install false." >> "\$STATUS_FILE"
        echo "Required packages include: ${packages.join(' ')}" >> "\$STATUS_FILE"
        exit 1
    fi

    INSTALLER="\$(find_installer)"

    if [[ -z "\$INSTALLER" ]]; then
        echo "ERROR: Neither mamba nor conda found in PATH." >> "\$STATUS_FILE"
        echo "Required packages include: ${packages.join(' ')}" >> "\$STATUS_FILE"
        exit 1
    fi

    echo "Using installer: \$INSTALLER" >> "\$STATUS_FILE"

    mkdir -p "\$(dirname "\$TOOL_ENV")"

    "\$INSTALLER" create -y \\
        -p "\$TOOL_ENV" \\
        -c conda-forge \\
        -c bioconda \\
        -c r \\
        ${package_string} \\
        >> "\$STATUS_FILE" 2>&1

    if ! check_env "\$TOOL_ENV"; then
        echo "ERROR: Newly created Module 4 environment failed checks." >> "\$STATUS_FILE"
        echo "Required packages include: ${packages.join(' ')}" >> "\$STATUS_FILE"
        exit 1
    fi

    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Module 4 tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}

process PREPARE_MAG_COLLECTION {
    tag "prepare_mag_collection"

    publishDir "${params.outdir}/gathered_bins", mode: params.publish_gathered_bins_mode, pattern: "gathered_bins/*", saveAs: { filename -> filename.replaceFirst(/^gathered_bins\//, '') }

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "*.tsv"

    publishDir "${params.outdir}/concat", mode: 'copy', pattern: "mag_contigs.fa"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "prepare_mag_collection.log"

    input:
    path binning_manifest
    path tools_status

    output:
    path "gathered_bins/*", optional: true, emit: gathered_bins
    path "mag_contigs.fa", emit: concatenated_fasta
    path "mag_contigs.contigs_to_bin.tsv", emit: contigs_to_bin
    path "mag_contigs.quickbin.contigs_to_bin.tsv", emit: quickbin_map
    path "mag_contigs.metabat2.contigs_to_bin.tsv", emit: metabat2_map
    path "mag_contigs.maxbin2.contigs_to_bin.tsv", emit: maxbin2_map
    path "mag_collection_manifest.tsv", emit: collection_manifest
    path "mag_collection_stats.tsv", emit: collection_stats
    path "prepare_mag_collection.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG_FILE="prepare_mag_collection.log"

    echo "MAG collection preparation started: \$(date)" > "\$LOG_FILE"
    echo "Binning manifest: ${binning_manifest}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    mkdir -p gathered_bins

    python3 - "${binning_manifest}" <<'PY'
import csv
import gzip
import re
import shutil
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1])

log_path = Path("prepare_mag_collection.log")
gathered_dir = Path("gathered_bins")
concat_path = Path("mag_contigs.fa")

combined_map_path = Path("mag_contigs.contigs_to_bin.tsv")
quickbin_map_path = Path("mag_contigs.quickbin.contigs_to_bin.tsv")
metabat2_map_path = Path("mag_contigs.metabat2.contigs_to_bin.tsv")
maxbin2_map_path = Path("mag_contigs.maxbin2.contigs_to_bin.tsv")

collection_manifest_path = Path("mag_collection_manifest.tsv")
collection_stats_path = Path("mag_collection_stats.tsv")

def log(msg):
    with log_path.open("a") as handle:
        print(msg, file=handle)

def sanitize_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("_")
    return value or "unnamed_mag"

def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def fasta_records(path):
    header = None
    seq_parts = []

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\\n")
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if header is not None:
        yield header, "".join(seq_parts)

def first_token(header):
    return header.split()[0] if header.strip() else "unnamed_contig"

methods = ["quickbin", "metabat2", "maxbin2"]

for p in [combined_map_path, quickbin_map_path, metabat2_map_path, maxbin2_map_path, concat_path]:
    p.write_text("")

rows = []

if not manifest_path.exists():
    raise SystemExit(f"ERROR: Binning manifest does not exist: {manifest_path}")

with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\\t")
    for row in reader:
        bin_fasta = row.get("bin_fasta", "").strip()
        bin_id = row.get("bin_id", "").strip()
        binner = row.get("binner", "").strip()

        if not bin_fasta:
            continue

        rows.append(row)

total_manifest_rows = len(rows)
total_bins_copied = 0
total_contigs = 0
total_bp = 0
missing_bins = 0

with collection_manifest_path.open("w") as collection_manifest, \
     combined_map_path.open("w") as combined_map, \
     quickbin_map_path.open("w") as quickbin_map, \
     metabat2_map_path.open("w") as metabat2_map, \
     maxbin2_map_path.open("w") as maxbin2_map, \
     concat_path.open("w") as concat:

    print(
        "mag_id",
        "binner",
        "source_bin_fasta",
        "gathered_bin_fasta",
        "contig_count",
        "total_bp",
        sep="\\t",
        file=collection_manifest
    )

    method_handles = {
        "quickbin": quickbin_map,
        "metabat2": metabat2_map,
        "maxbin2": maxbin2_map
    }

    used_mag_ids = set()

    for idx, row in enumerate(rows, start=1):
        source_fasta = Path(row.get("bin_fasta", "").strip())
        binner = row.get("binner", "").strip()
        raw_bin_id = row.get("bin_id", "").strip() or f"bin_{idx:06d}"
        mag_id_base = sanitize_id(raw_bin_id)

        mag_id = mag_id_base
        suffix = 1
        while mag_id in used_mag_ids:
            suffix += 1
            mag_id = f"{mag_id_base}_{suffix}"

        used_mag_ids.add(mag_id)

        if not source_fasta.exists():
            log(f"WARNING: source bin FASTA does not exist and will be skipped: {source_fasta}")
            missing_bins += 1
            continue

        gathered_fasta = gathered_dir / f"{mag_id}.fa"

        contig_count = 0
        bp_count = 0

        with gathered_fasta.open("w") as out_fa:
            for header, seq in fasta_records(source_fasta):
                if not seq:
                    continue

                contig_id = first_token(header)
                contig_count += 1
                bp_count += len(seq)

                print(f">{header}", file=out_fa)

                for start in range(0, len(seq), 80):
                    print(seq[start:start+80], file=out_fa)

                print(f">{header}", file=concat)

                for start in range(0, len(seq), 80):
                    print(seq[start:start+80], file=concat)

                print(mag_id, contig_id, binner, sep="\\t", file=combined_map)

                if binner in method_handles:
                    print(mag_id, contig_id, sep="\\t", file=method_handles[binner])

        if contig_count == 0:
            log(f"WARNING: source bin FASTA had zero usable contigs: {source_fasta}")

        total_bins_copied += 1
        total_contigs += contig_count
        total_bp += bp_count

        print(
            mag_id,
            binner,
            str(source_fasta),
            str(gathered_fasta),
            contig_count,
            bp_count,
            sep="\\t",
            file=collection_manifest
        )

with collection_stats_path.open("w") as stats:
    print(
        "total_manifest_rows",
        "total_bins_copied",
        "missing_bins",
        "total_contigs",
        "total_bp",
        "combined_fasta",
        "contigs_to_bin",
        sep="\\t",
        file=stats
    )
    print(
        total_manifest_rows,
        total_bins_copied,
        missing_bins,
        total_contigs,
        total_bp,
        "mag_contigs.fa",
        "mag_contigs.contigs_to_bin.tsv",
        sep="\\t",
        file=stats
    )

log(f"Total manifest rows: {total_manifest_rows}")
log(f"Total bins copied: {total_bins_copied}")
log(f"Missing bins: {missing_bins}")
log(f"Total contigs: {total_contigs}")
log(f"Total bp: {total_bp}")
log("MAG collection preparation finished.")
PY

    echo "MAG collection preparation finished: \$(date)" >> "\$LOG_FILE"
    """
}


process RUN_PRODIGAL_ON_MAG_CONTIGS {
    tag "prodigal_mag_contigs"

    publishDir "${params.outdir}/prodigal", mode: 'copy', pattern: "mag_contigs.prodigal.*"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "prodigal.log"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "prodigal_status.tsv"

    input:
    path mag_contigs_fasta
    path tools_status

    output:
    path "mag_contigs.prodigal.faa", emit: proteins
    path "mag_contigs.prodigal.ffn", emit: nucleotides
    path "mag_contigs.prodigal.gff", emit: prodigal_gff
    path "prodigal_status.tsv", emit: prodigal_status
    path "prodigal.log", emit: log_file

    script:
    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" && "\$TOOL_ENV" != "NOT_USED" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    LOG_FILE="prodigal.log"

    echo "Prodigal started: \$(date)" > "\$LOG_FILE"
    echo "Input FASTA: ${mag_contigs_fasta}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    touch mag_contigs.prodigal.faa
    touch mag_contigs.prodigal.ffn
    touch mag_contigs.prodigal.gff

    SEQ_COUNT="\$(grep -c '^>' "${mag_contigs_fasta}" || true)"

    if [[ "\$SEQ_COUNT" -eq 0 ]]; then
        echo "WARNING: No contigs present in ${mag_contigs_fasta}; skipping Prodigal." >> "\$LOG_FILE"

        printf 'step\\tstatus\\texit_status\\tmessage\\n' > prodigal_status.tsv
        printf 'prodigal\\tskipped\\t0\\tNo contigs available for Prodigal\\n' >> prodigal_status.tsv

        exit 0
    fi

    if ! command -v prodigal >/dev/null 2>&1; then
        echo "ERROR: prodigal is not available after tool setup." >> "\$LOG_FILE"
        exit 1
    fi

    set +e
    prodigal \\
        -p meta \\
        -a mag_contigs.prodigal.faa \\
        -d mag_contigs.prodigal.ffn \\
        -o mag_contigs.prodigal.gff \\
        -i "${mag_contigs_fasta}" \\
        >> "\$LOG_FILE" 2>&1

    STATUS="\$?"
    set -e

    if [[ "\$STATUS" -ne 0 ]]; then
        echo "ERROR: Prodigal failed with status \$STATUS" >> "\$LOG_FILE"
        printf 'step\\tstatus\\texit_status\\tmessage\\n' > prodigal_status.tsv
        printf 'prodigal\\tfailed\\t%s\\tProdigal failed\\n' "\$STATUS" >> prodigal_status.tsv
        exit "\$STATUS"
    fi

    printf 'step\\tstatus\\texit_status\\tmessage\\n' > prodigal_status.tsv
    printf 'prodigal\\tcompleted\\t0\\tProdigal completed successfully\\n' >> prodigal_status.tsv

    echo "Prodigal finished: \$(date)" >> "\$LOG_FILE"
    """
}


process RUN_HMMSEARCH_MAG_CONTIGS {
    tag "hmmsearch_mag_contigs"

    publishDir "${params.outdir}/hmm", mode: 'copy', pattern: "mag_contigs.hmm*"

    publishDir "${params.outdir}/hmm", mode: 'copy', pattern: "mag_contigs.*fam*"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "hmmsearch.log"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "hmm_status.tsv"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.hmm_threads as int
    }

    input:
    path proteins
    path tigrfam_hmm
    path pfam_hmm
    path tools_status

    output:
    path "mag_contigs.hmm", emit: hmm_table
    path "mag_contigs.tigr", emit: tigr_table
    path "mag_contigs.pfam", emit: pfam_table
    path "mag_contigs.hmm.tigr.out", emit: tigr_out
    path "mag_contigs.hmm.tigr.hit.out", emit: tigr_tblout
    path "mag_contigs.hmm.pfam.out", emit: pfam_out
    path "mag_contigs.hmm.pfam.hit.out", emit: pfam_tblout
    path "hmm_status.tsv", emit: hmm_status
    path "hmmsearch.log", emit: log_file

    script:
    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" && "\$TOOL_ENV" != "NOT_USED" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    LOG_FILE="hmmsearch.log"

    echo "HMM searches started: \$(date)" > "\$LOG_FILE"
    echo "Proteins: ${proteins}" >> "\$LOG_FILE"
    echo "TIGRFAM HMM: ${tigrfam_hmm}" >> "\$LOG_FILE"
    echo "Pfam HMM: ${pfam_hmm}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    touch mag_contigs.hmm.tigr.out
    touch mag_contigs.hmm.tigr.hit.out
    touch mag_contigs.hmm.pfam.out
    touch mag_contigs.hmm.pfam.hit.out
    touch mag_contigs.tigr
    touch mag_contigs.pfam
    touch mag_contigs.hmm

    PROTEIN_COUNT="\$(grep -c '^>' "${proteins}" || true)"

    if [[ "\$PROTEIN_COUNT" -eq 0 ]]; then
        echo "WARNING: No proteins present; skipping HMM searches." >> "\$LOG_FILE"

        printf 'step\\tstatus\\texit_status\\tmessage\\n' > hmm_status.tsv
        printf 'hmmsearch\\tskipped\\t0\\tNo proteins available for HMM search\\n' >> hmm_status.tsv

        exit 0
    fi

    if ! command -v hmmsearch >/dev/null 2>&1; then
        echo "ERROR: hmmsearch is not available after tool setup." >> "\$LOG_FILE"
        exit 1
    fi

    hmmsearch \\
        -o mag_contigs.hmm.tigr.out \\
        --tblout mag_contigs.hmm.tigr.hit.out \\
        --noali \\
        --notextw \\
        --cut_nc \\
        --cpu ${task.cpus} \\
        "${tigrfam_hmm}" \\
        "${proteins}" \\
        >> "\$LOG_FILE" 2>&1

    hmmsearch \\
        -o mag_contigs.hmm.pfam.out \\
        --tblout mag_contigs.hmm.pfam.hit.out \\
        --noali \\
        --notextw \\
        --cut_nc \\
        --cpu ${task.cpus} \\
        "${pfam_hmm}" \\
        "${proteins}" \\
        >> "\$LOG_FILE" 2>&1

    grep -v '^#' mag_contigs.hmm.tigr.hit.out \\
        | awk 'NF >= 5 {print \$1"\\t"\$3"\\t"\$5}' \\
        > mag_contigs.tigr || true

    grep -v '^#' mag_contigs.hmm.pfam.hit.out \\
        | awk 'NF >= 5 {print \$1"\\t"\$4"\\t"\$5}' \\
        > mag_contigs.pfam || true

    cat mag_contigs.pfam mag_contigs.tigr > mag_contigs.hmm

    TIGR_HITS="\$(wc -l < mag_contigs.tigr | tr -d ' ')"
    PFAM_HITS="\$(wc -l < mag_contigs.pfam | tr -d ' ')"
    TOTAL_HITS="\$(wc -l < mag_contigs.hmm | tr -d ' ')"

    printf 'step\\tstatus\\texit_status\\tmessage\\ttigr_hits\\tpfam_hits\\ttotal_hits\\n' > hmm_status.tsv
    printf 'hmmsearch\\tcompleted\\t0\\tHMM searches completed\\t%s\\t%s\\t%s\\n' "\$TIGR_HITS" "\$PFAM_HITS" "\$TOTAL_HITS" >> hmm_status.tsv

    echo "TIGRFAM hits: \$TIGR_HITS" >> "\$LOG_FILE"
    echo "Pfam hits: \$PFAM_HITS" >> "\$LOG_FILE"
    echo "Total HMM hits: \$TOTAL_HITS" >> "\$LOG_FILE"
    echo "HMM searches finished: \$(date)" >> "\$LOG_FILE"
    """
}

process RUN_MAGSCOT {
    tag "run_magscot"

    publishDir "${params.outdir}/magscot", mode: 'copy', pattern: "magscot_outputs/**", saveAs: { filename -> filename.replaceFirst(/^magscot_outputs\//, '') }

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "magscot.log"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "magscot_status.tsv"

    input:
    path contigs_to_bin
    path hmm_table
    path magscot_script
    path magscot_profiles_dir
    path tools_status

    output:
    path "magscot_outputs", emit: magscot_outputs
    path "magscot_status.tsv", emit: magscot_status
    path "magscot.log", emit: log_file

    script:

    def magscot_threshold_arg = params.magscot_threshold != null && params.magscot_threshold.toString().trim()
        ? "--threshold ${params.magscot_threshold}"
        : ""
    """

set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" && "\$TOOL_ENV" != "NOT_USED" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    LOG_FILE="magscot.log"

    mkdir -p magscot_outputs

    echo "MAGScoT started: \$(date)" > "\$LOG_FILE"
    echo "Contigs-to-bin table: ${contigs_to_bin}" >> "\$LOG_FILE"
    echo "HMM table: ${hmm_table}" >> "\$LOG_FILE"
    echo "MAGScoT script: ${magscot_script}" >> "\$LOG_FILE"
    echo "MAGScoT threshold: ${params.magscot_threshold}" >> "\$LOG_FILE"
    echo "MAGScoT extra args: ${params.magscot_extra_args}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    CONTIG_MAP_ROWS="\$(wc -l < "${contigs_to_bin}" | tr -d ' ')"
    HMM_ROWS="\$(wc -l < "${hmm_table}" | tr -d ' ')"

    if [[ "\$CONTIG_MAP_ROWS" -eq 0 ]]; then
        echo "WARNING: contigs-to-bin table is empty; skipping MAGScoT." >> "\$LOG_FILE"
        printf 'step\\tstatus\\texit_status\\tmessage\\tcontigs_to_bin_rows\\thmm_rows\\n' > magscot_status.tsv
        printf 'magscot\\tskipped\\t0\\tNo MAG/bin contig mappings available\\t%s\\t%s\\n' "\$CONTIG_MAP_ROWS" "\$HMM_ROWS" >> magscot_status.tsv
        exit 0
    fi

    if [[ "\$HMM_ROWS" -eq 0 ]]; then
        echo "WARNING: HMM table is empty; skipping MAGScoT." >> "\$LOG_FILE"
        printf 'step\\tstatus\\texit_status\\tmessage\\tcontigs_to_bin_rows\\thmm_rows\\n' > magscot_status.tsv
        printf 'magscot\\tskipped\\t0\\tNo HMM hits available\\t%s\\t%s\\n' "\$CONTIG_MAP_ROWS" "\$HMM_ROWS" >> magscot_status.tsv
        exit 0
    fi

    if [[ ! -s "${magscot_script}" ]]; then
        echo "ERROR: MAGScoT.py script is missing or empty: ${magscot_script}" >> "\$LOG_FILE"
        exit 1
    fi

    if [[ ! -d "${magscot_profiles_dir}" ]]; then
        echo "ERROR: MAGScoT profiles directory is missing: ${magscot_profiles_dir}" >> "\$LOG_FILE"
        exit 1
    fi

    if [[ ! -s "${magscot_profiles_dir}/gtdb_rel207_default_markers.tsv" ]]; then
        echo "ERROR: Required MAGScoT profile file is missing:" >> "\$LOG_FILE"
        echo "  ${magscot_profiles_dir}/gtdb_rel207_default_markers.tsv" >> "\$LOG_FILE"
        echo "" >> "\$LOG_FILE"
        echo "Contents of provided profiles directory:" >> "\$LOG_FILE"
        ls -lah "${magscot_profiles_dir}" >> "\$LOG_FILE" 2>&1 || true
        exit 1
    fi

    if command -v python >/dev/null 2>&1; then
        PYTHON_EXE="python"
    elif command -v python3 >/dev/null 2>&1; then
        PYTHON_EXE="python3"
    else
        echo "ERROR: neither python nor python3 is available after tool setup." >> "\$LOG_FILE"
        exit 1
    fi

    cp -L "${contigs_to_bin}" magscot_outputs/mag_contigs.contigs_to_bin.tsv
    cp -L "${hmm_table}" magscot_outputs/mag_contigs.hmm
    cp -L "${magscot_script}" magscot_outputs/MAGScoT.py

    # MAGScoT expects this path relative to its working directory:
    #   profiles/gtdb_rel207_default_markers.tsv
    cp -RL "${magscot_profiles_dir}/." magscot_outputs/profiles/

    echo "Copied MAGScoT profiles into magscot_outputs/profiles" >> "\$LOG_FILE"
    ls -lah magscot_outputs/profiles >> "\$LOG_FILE" 2>&1 || true

    cd magscot_outputs

    set +e
    "\$PYTHON_EXE" MAGScoT.py \\
    -i mag_contigs.contigs_to_bin.tsv \\
    --hmm mag_contigs.hmm \\
    ${magscot_threshold_arg} \\
    ${params.magscot_extra_args} \\
    >> "../\$LOG_FILE" 2>&1
STATUS="\$?"
set -e

    cd ..

    if [[ "\$STATUS" -ne 0 ]]; then
        echo "WARNING: MAGScoT exited non-zero; treating as non-fatal for workflow continuation. Exit status: \$STATUS" >> "\$LOG_FILE"
        printf 'step\\tstatus\\texit_status\\tmessage\\tcontigs_to_bin_rows\\thmm_rows\\n' > magscot_status.tsv
        printf 'magscot\\tfailed_nonfatal\\t%s\\tMAGScoT exited non-zero\\t%s\\t%s\\n' "\$STATUS" "\$CONTIG_MAP_ROWS" "\$HMM_ROWS" >> magscot_status.tsv
        exit 0
    fi

    printf 'step\\tstatus\\texit_status\\tmessage\\tcontigs_to_bin_rows\\thmm_rows\\n' > magscot_status.tsv
    printf 'magscot\\tcompleted\\t0\\tMAGScoT completed successfully\\t%s\\t%s\\n' "\$CONTIG_MAP_ROWS" "\$HMM_ROWS" >> magscot_status.tsv

    echo "MAGScoT finished: \$(date)" >> "\$LOG_FILE"
    """
}

process BUILD_REFINED_MAGS {
    tag "build_refined_mags"

    publishDir "${params.outdir}/refined_bins", mode: params.publish_refined_bins_mode, pattern: "refined_bins/*.fa", saveAs: { filename -> filename.replaceFirst(/^refined_bins\//, '') }

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "magscot_refined_bins_*.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "build_refined_mags.log"

    input:
    path mag_contigs_fasta
    path magscot_outputs_dir

    output:
    path "refined_bins/*.fa", optional: true, emit: refined_bins
    path "magscot_refined_bins_manifest.tsv", emit: refined_manifest
    path "magscot_refined_bins_stats.tsv", emit: refined_stats
    path "build_refined_mags.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG_FILE="build_refined_mags.log"
    REFINED_MAP="${magscot_outputs_dir}/MAGScoT.refined.contig_to_bin.out"

    echo "Refined MAG reconstruction started: \$(date)" > "\$LOG_FILE"
    echo "Input combined MAG contigs FASTA: ${mag_contigs_fasta}" >> "\$LOG_FILE"
    echo "MAGScoT outputs directory: ${magscot_outputs_dir}" >> "\$LOG_FILE"
    echo "Expected refined map: \$REFINED_MAP" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    mkdir -p refined_bins

    python3 - \\
        "${mag_contigs_fasta}" \\
        "\$REFINED_MAP" \\
        "refined_bins" \\
        "magscot_refined_bins_manifest.tsv" \\
        "magscot_refined_bins_stats.tsv" \\
        "\$LOG_FILE" \\
        "${params.outdir}/refined_bins" <<'PY'
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

(
    fasta_path,
    refined_map_path,
    refined_bins_dir,
    manifest_path,
    stats_path,
    log_path,
    published_refined_bins_dir
) = sys.argv[1:]

fasta_path = Path(fasta_path)
refined_map_path = Path(refined_map_path)
refined_bins_dir = Path(refined_bins_dir)
manifest_path = Path(manifest_path)
stats_path = Path(stats_path)
log_path = Path(log_path)
published_refined_bins_dir = Path(published_refined_bins_dir)

refined_bins_dir.mkdir(parents=True, exist_ok=True)

def log(message):
    with log_path.open("a") as handle:
        print(message, file=handle)

def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def first_token(header):
    header = header.strip()
    if not header:
        return ""
    return header.split()[0]

def wrap_seq(seq, width=80):
    for start in range(0, len(seq), width):
        yield seq[start:start + width]

def safe_filename(value):
    value = str(value).strip()
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    safe = safe.strip("_")
    return safe or "unnamed_refined_bin"

def iter_fasta(path):
    header = None
    seq_parts = []

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\\n")
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())

    if header is not None:
        yield header, "".join(seq_parts)

# ---------------------------------------------------------------------
# Read MAGScoT refined assignment table.
#
# Expected format:
#
#   binnew contig
#   MAGScoT_cleanbin_000001 example_contig_13943
#   MAGScoT_cleanbin_000001 example_contig_30913
#
# Whitespace-delimited is accepted, so tabs or spaces both work.
# ---------------------------------------------------------------------

assignments_by_contig = defaultdict(list)
contigs_by_bin = defaultdict(list)

map_rows_seen = 0
map_rows_used = 0
duplicate_assignment_rows = 0

if not refined_map_path.exists():
    log(f"WARNING: Refined MAGScoT map does not exist: {refined_map_path}")
else:
    with refined_map_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line:
                continue

            fields = line.split()

            if len(fields) < 2:
                log(f"WARNING: Skipping malformed refined map line {line_number}: {line}")
                continue

            binnew = fields[0].strip()
            contig = fields[1].strip()

            # Skip header.
            if binnew.lower() == "binnew" and contig.lower() == "contig":
                continue

            map_rows_seen += 1

            if not binnew or not contig:
                continue

            if binnew in assignments_by_contig[contig]:
                duplicate_assignment_rows += 1
                continue

            assignments_by_contig[contig].append(binnew)
            contigs_by_bin[binnew].append(contig)
            map_rows_used += 1

log(f"Refined map rows seen: {map_rows_seen}")
log(f"Refined map rows used: {map_rows_used}")
log(f"Duplicate refined assignment rows skipped: {duplicate_assignment_rows}")
log(f"Refined bins in map: {len(contigs_by_bin)}")

# ---------------------------------------------------------------------
# If there are no refined assignments, write empty outputs and finish.
# ---------------------------------------------------------------------

if map_rows_used == 0:
    with manifest_path.open("w") as manifest:
        print(
            "refined_bin_id",
            "refined_bin_fasta",
            "contig_count",
            "total_bp",
            sep="\\t",
            file=manifest
        )

    with stats_path.open("w") as stats:
        print(
            "refined_bin_count",
            "refined_assignment_rows_seen",
            "refined_assignment_rows_used",
            "input_fasta_contigs_seen",
            "unique_input_contigs_seen",
            "duplicate_input_contig_headers",
            "assigned_contigs_requested",
            "assigned_contigs_recovered",
            "assigned_contigs_missing",
            "total_refined_bin_bp",
            "status",
            "message",
            sep="\\t",
            file=stats
        )
        print(
            0,
            map_rows_seen,
            map_rows_used,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            "skipped",
            "No refined MAGScoT assignments were available",
            sep="\\t",
            file=stats
        )

    log("No refined assignments available; no refined MAG FASTA files were created.")
    raise SystemExit(0)

# ---------------------------------------------------------------------
# Rebuild refined bins from combined FASTA.
#
# Important:
# The combined FASTA may contain duplicated contig headers if the same
# contig was present in bins from multiple original binners. For refined
# MAG reconstruction, each contig ID is only used once. Subsequent
# duplicate headers are skipped and reported.
# ---------------------------------------------------------------------

bin_file_handles = {}
bin_file_paths = {}
bin_contig_counts = defaultdict(int)
bin_bp_counts = defaultdict(int)

input_fasta_contigs_seen = 0
unique_input_contigs_seen = 0
duplicate_input_contig_headers = 0

seen_input_contig_ids = set()
recovered_assignment_pairs = set()
recovered_contigs = set()

def get_bin_handle(bin_id):
    if bin_id in bin_file_handles:
        return bin_file_handles[bin_id]

    filename = safe_filename(bin_id) + ".fa"
    out_path = refined_bins_dir / filename

    handle = out_path.open("w")
    bin_file_handles[bin_id] = handle
    bin_file_paths[bin_id] = out_path

    return handle

if not fasta_path.exists():
    raise SystemExit(f"ERROR: Combined FASTA does not exist: {fasta_path}")

for header, seq in iter_fasta(fasta_path):
    input_fasta_contigs_seen += 1

    contig_id = first_token(header)

    if not contig_id:
        log("WARNING: Encountered FASTA record with empty contig ID; skipping.")
        continue

    if contig_id in seen_input_contig_ids:
        duplicate_input_contig_headers += 1
        continue

    seen_input_contig_ids.add(contig_id)
    unique_input_contigs_seen += 1

    if contig_id not in assignments_by_contig:
        continue

    recovered_contigs.add(contig_id)

    for bin_id in assignments_by_contig[contig_id]:
        out_handle = get_bin_handle(bin_id)

        print(f">{header}", file=out_handle)
        for chunk in wrap_seq(seq):
            print(chunk, file=out_handle)

        bin_contig_counts[bin_id] += 1
        bin_bp_counts[bin_id] += len(seq)
        recovered_assignment_pairs.add((bin_id, contig_id))

for handle in bin_file_handles.values():
    handle.close()

assigned_contigs_requested = len(assignments_by_contig)
assigned_contigs_recovered = len(recovered_contigs)
assigned_contigs_missing = assigned_contigs_requested - assigned_contigs_recovered
total_refined_bin_bp = sum(bin_bp_counts.values())

missing_contigs = sorted(set(assignments_by_contig.keys()) - recovered_contigs)

if missing_contigs:
    log(f"WARNING: {len(missing_contigs)} contig IDs from MAGScoT refined map were not found in combined FASTA.")
    preview = missing_contigs[:25]
    for contig in preview:
        log(f"  missing_contig: {contig}")
    if len(missing_contigs) > len(preview):
        log(f"  ... {len(missing_contigs) - len(preview)} more missing contigs not shown")

# ---------------------------------------------------------------------
# Write refined-bin manifest.
# ---------------------------------------------------------------------

with manifest_path.open("w") as manifest:
    print(
        "refined_bin_id",
        "refined_bin_fasta",
        "contig_count",
        "total_bp",
        sep="\\t",
        file=manifest
    )

    for bin_id in sorted(bin_file_paths):
        fasta_path_out = bin_file_paths[bin_id]
        published_path = published_refined_bins_dir / fasta_path_out.name

        print(
            bin_id,
            str(published_path),
            bin_contig_counts[bin_id],
            bin_bp_counts[bin_id],
            sep="\\t",
            file=manifest
        )

# ---------------------------------------------------------------------
# Write stats.
# ---------------------------------------------------------------------

with stats_path.open("w") as stats:
    print(
        "refined_bin_count",
        "refined_assignment_rows_seen",
        "refined_assignment_rows_used",
        "input_fasta_contigs_seen",
        "unique_input_contigs_seen",
        "duplicate_input_contig_headers",
        "assigned_contigs_requested",
        "assigned_contigs_recovered",
        "assigned_contigs_missing",
        "total_refined_bin_bp",
        "status",
        "message",
        sep="\\t",
        file=stats
    )

    status = "completed"
    message = "Refined MAG FASTA files were reconstructed"

    if len(bin_file_paths) == 0:
        status = "no_refined_bins"
        message = "Refined assignments existed, but no refined bins were reconstructed"

    print(
        len(bin_file_paths),
        map_rows_seen,
        map_rows_used,
        input_fasta_contigs_seen,
        unique_input_contigs_seen,
        duplicate_input_contig_headers,
        assigned_contigs_requested,
        assigned_contigs_recovered,
        assigned_contigs_missing,
        total_refined_bin_bp,
        status,
        message,
        sep="\\t",
        file=stats
    )

log(f"Input FASTA contigs seen: {input_fasta_contigs_seen}")
log(f"Unique input FASTA contigs seen: {unique_input_contigs_seen}")
log(f"Duplicate input FASTA contig headers skipped: {duplicate_input_contig_headers}")
log(f"Assigned contigs requested: {assigned_contigs_requested}")
log(f"Assigned contigs recovered: {assigned_contigs_recovered}")
log(f"Assigned contigs missing: {assigned_contigs_missing}")
log(f"Refined bins written: {len(bin_file_paths)}")
log(f"Total refined bin bp: {total_refined_bin_bp}")
log("Refined MAG reconstruction finished.")
PY

    echo "Refined MAG reconstruction finished: \$(date)" >> "\$LOG_FILE"
    """
}

process WRITE_MODULE4_SUMMARY {
    tag "write_module4_summary"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "module4_run_summary.tsv"

    input:
    path collection_stats
    path prodigal_status
    path hmm_status
    path magscot_status
    path refined_stats

    output:
    path "module4_run_summary.tsv", emit: summary

    script:
    """
    set -euo pipefail

    printf 'section\\tsource_file\\n' > module4_run_summary.tsv
    printf 'mag_collection\\t%s\\n' "${collection_stats}" >> module4_run_summary.tsv
    printf 'prodigal\\t%s\\n' "${prodigal_status}" >> module4_run_summary.tsv
    printf 'hmmsearch\\t%s\\n' "${hmm_status}" >> module4_run_summary.tsv
    printf 'magscot\\t%s\\n' "${magscot_status}" >> module4_run_summary.tsv
    printf 'refined_mags\\t%s\\n' "${refined_stats}" >> module4_run_summary.tsv

    echo "" >> module4_run_summary.tsv
    echo "# mag_collection_stats.tsv" >> module4_run_summary.tsv
    cat "${collection_stats}" >> module4_run_summary.tsv

    echo "" >> module4_run_summary.tsv
    echo "# prodigal_status.tsv" >> module4_run_summary.tsv
    cat "${prodigal_status}" >> module4_run_summary.tsv

    echo "" >> module4_run_summary.tsv
    echo "# hmm_status.tsv" >> module4_run_summary.tsv
    cat "${hmm_status}" >> module4_run_summary.tsv

    echo "" >> module4_run_summary.tsv
    echo "# magscot_status.tsv" >> module4_run_summary.tsv
    cat "${magscot_status}" >> module4_run_summary.tsv

    echo "" >> module4_run_summary.tsv
    echo "# magscot_refined_bins_stats.tsv" >> module4_run_summary.tsv
    cat "${refined_stats}" >> module4_run_summary.tsv
    """
}