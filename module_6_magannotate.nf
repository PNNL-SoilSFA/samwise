#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Module 6: MAG harmonization and processing.
 */


/*
 * Core parameters
 */
params.working_dir = null
params.output_dir  = null

params.input_mag_dir      = null
params.input_mag_manifest = null
params.mag_extension      = "fa"

params.threads = null
params.tool_env_dir = null
params.auto_install = true

/*
 * Which tools to run
 */
params.run_checkm2    = true
params.run_gtdbtk     = true
params.run_dram2      = true
params.run_microtrait = true

/*
 * Tool versions
 */
params.checkm2_version = null
params.gtdbtk_version  = "2.7.2"

/*
 * CheckM2 database.
 *
 * If --checkm2_db_path is supplied, that file is used directly.
 *
 * If not supplied, this module attempts to download the database from
 * Zenodo record 14897628 into:
 *
 *   <outdir>/databases/checkm2/
 */
params.checkm2_db_path = null
params.checkm2_db_dir  = null
params.checkm2_zenodo_record = "14897628"
params.checkm2_auto_download_db = true
params.checkm2_extension = "fa"

/*
 * GTDB-Tk database.
 *
 * If --gtdbtk_data_path is supplied, that directory is used.
 *
 * Otherwise the workflow attempts to run download-db.sh into:
 *
 *   <outdir>/databases/gtdbtk/
 */
params.gtdbtk_data_path = null
params.gtdbtk_db_dir = null
params.gtdbtk_auto_download_db = true
params.gtdbtk_extension = "fa"

/*
 * DRAM2.
 *
 * DRAM2 is itself a Nextflow workflow.
 */
params.dram2_work_dir = null
params.dram2_profile = "apptainer"
params.dram2_repo = "WrightonLabCSU/DRAM"
params.dram2_revision = "dev"
params.dram2_nextflow_config_url = "https://raw.githubusercontent.com/WrightonLabCSU/DRAM/refs/heads/dev/nextflow.config"

params.dram2_tiny_cpus_limit   = 1
params.dram2_small_cpus_limit  = 12
params.dram2_medium_cpus_limit = 24
params.dram2_big_cpus_limit    = 36
params.dram2_huge_cpus_limit   = 36

params.dram2_tiny_gb_mem_limit   = 1
params.dram2_small_gb_mem_limit  = 100
params.dram2_medium_gb_mem_limit = 200
params.dram2_big_gb_mem_limit    = 300
params.dram2_huge_gb_mem_limit   = 360

params.dram2_tiny_hr_time_limit   = 12
params.dram2_small_hr_time_limit  = 120
params.dram2_medium_hr_time_limit = 120
params.dram2_big_hr_time_limit    = 168
params.dram2_huge_hr_time_limit   = 168

params.nextflow_exe = "nextflow"

/*
 * microTrait.
 *
 * If --microtrait_env_dir is supplied, that environment is used.
 * Otherwise the module attempts to create an environment.
 *
 * If --microtrait_runner is supplied, that Rscript is used.
 * Otherwise the module writes a bundled runner script based on the code
 * you provided.
 *
 * If your microTrait helper R files live somewhere specific, use:
 *   --microtrait_source_dir /path/to/microtrait-022624/R
 */
params.microtrait_env_dir = null
params.microtrait_runner = null
params.microtrait_source_dir = null
params.microtrait_out_name = "samwise_final_mags"
params.microtrait_run_type = "genomic"
params.microtrait_cores = null

/*
 * Publish modes
 */
params.publish_mags_mode = "copy"
params.publish_tool_outputs_mode = "copy"

/*
 * Derived directories
 */
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.outdir      = "${params.results_dir}/module_6_magHarmonization"

params.module5_final_mag_dir = "${params.results_dir}/module_5_subtractiveAssembly/final_mag_database"
params.module4_refined_mag_dir = "${params.results_dir}/module_4_binRefinement/refined_bins"

params.checkm2_db_outdir = params.checkm2_db_dir ?: "${params.outdir}/databases/checkm2"
params.gtdbtk_db_outdir  = params.gtdbtk_db_dir  ?: "${params.outdir}/databases/gtdbtk"
params.dram2_outdir      = params.dram2_work_dir ?: "${params.outdir}/DRAM"


def firstExistingPath(List candidates) {
    def found = candidates.find { candidate ->
        java.nio.file.Files.exists(java.nio.file.Paths.get(candidate.toString()))
    }

    def selected = found ?: candidates[0]

    return java.nio.file.Paths
        .get(selected.toString())
        .toAbsolutePath()
        .normalize()
        .toString()
}


workflow {

    def run_checkm2    = params.run_checkm2.toString().toBoolean()
    def run_gtdbtk     = params.run_gtdbtk.toString().toBoolean()
    def run_dram2      = params.run_dram2.toString().toBoolean()
    def run_microtrait = params.run_microtrait.toString().toBoolean()

    if( !run_checkm2 && !run_gtdbtk && !run_dram2 && !run_microtrait ) {
        error """
        No Module 6 tools selected.

        Enable at least one of:
          --run_checkm2 true
          --run_gtdbtk true
          --run_dram2 true
          --run_microtrait true
        """.stripIndent()
    }

    def selected_mag_dir = params.input_mag_dir ?: firstExistingPath([
        params.module5_final_mag_dir,
        params.module4_refined_mag_dir
    ])

    def selected_mag_manifest = params.input_mag_manifest ?: ""

    log.info "Module 6 results directory: ${params.results_dir}"
    log.info "Writing Module 6 outputs to: ${params.outdir}"
    log.info "Selected MAG directory: ${selected_mag_dir}"
    log.info "Selected MAG manifest: ${selected_mag_manifest ?: 'not supplied'}"
    log.info "MAG extension: ${params.mag_extension}"
    log.info "Threads: ${params.threads ?: 'tool-specific defaults'}"
    log.info "Run CheckM2: ${run_checkm2}"
    log.info "Run GTDB-Tk: ${run_gtdbtk}"
    log.info "Run DRAM2: ${run_dram2}"
    log.info "Run microTrait: ${run_microtrait}"

    PREPARE_MAG_INPUTS(
        channel.value(selected_mag_dir),
        channel.value(selected_mag_manifest)
    )

    def status_ch = channel.empty()

    if( run_checkm2 ) {
        SETUP_CHECKM2()
        RUN_CHECKM2(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_CHECKM2.out.status
        )
        status_ch = status_ch.mix(RUN_CHECKM2.out.status)
    }

    if( run_gtdbtk ) {
        SETUP_GTDBTK()
        RUN_GTDBTK(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_GTDBTK.out.status
        )
        status_ch = status_ch.mix(RUN_GTDBTK.out.status)
    }

    if( run_dram2 ) {
        RUN_DRAM2(
            PREPARE_MAG_INPUTS.out.mags_dir
        )
        status_ch = status_ch.mix(RUN_DRAM2.out.status)
    }

    if( run_microtrait ) {
        SETUP_MICROTRAIT()
        RUN_MICROTRAIT(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_MICROTRAIT.out.status
        )
        status_ch = status_ch.mix(RUN_MICROTRAIT.out.status)
    }

    WRITE_MODULE6_SUMMARY(
        PREPARE_MAG_INPUTS.out.input_stats,
        status_ch.collect()
    )
}


process PREPARE_MAG_INPUTS {
    tag "prepare_final_mags"

    publishDir "${params.outdir}/refined_genomes",
        mode: params.publish_mags_mode,
        pattern: "refined_genomes/*.fa",
        saveAs: { filename -> filename.replaceFirst(/^refined_genomes\//, '') }

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "module6_mag_input_*.tsv"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "prepare_mag_inputs.log"

    input:
    val input_mag_dir
    val input_mag_manifest

    output:
    path "refined_genomes", emit: mags_dir
    path "module6_mag_input_manifest.tsv", emit: input_manifest
    path "module6_mag_input_stats.tsv", emit: input_stats
    path "prepare_mag_inputs.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG_FILE="prepare_mag_inputs.log"

    echo "Preparing MAG inputs: \$(date)" > "\$LOG_FILE"
    echo "Input MAG directory: ${input_mag_dir}" >> "\$LOG_FILE"
    echo "Input MAG manifest: ${input_mag_manifest}" >> "\$LOG_FILE"
    echo "MAG extension: ${params.mag_extension}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    mkdir -p refined_genomes

    python3 - \\
        "${input_mag_dir}" \\
        "${input_mag_manifest}" \\
        "${params.mag_extension}" \\
        "refined_genomes" \\
        "module6_mag_input_manifest.tsv" \\
        "module6_mag_input_stats.tsv" \\
        "\$LOG_FILE" \\
        "${params.outdir}/refined_genomes" <<'PY'
import csv
import gzip
import re
import shutil
import sys
from pathlib import Path

(
    input_mag_dir,
    input_mag_manifest,
    mag_extension,
    out_dir,
    out_manifest,
    out_stats,
    log_file,
    published_dir
) = sys.argv[1:]

input_mag_dir = Path(input_mag_dir)
input_mag_manifest = input_mag_manifest.strip()
mag_extension = mag_extension.lstrip(".")
out_dir = Path(out_dir)
out_manifest = Path(out_manifest)
out_stats = Path(out_stats)
log_file = Path(log_file)
published_dir = Path(published_dir)

out_dir.mkdir(parents=True, exist_ok=True)

def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)

def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("_")
    return value or "unnamed_mag"

def fasta_stats(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs = 0
    bp = 0
    cur = 0
    seen = False

    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if seen:
                    bp += cur
                cur = 0
                seen = True
                contigs += 1
            else:
                cur += len(line.strip())

    if seen:
        bp += cur

    return contigs, bp

def candidate_fastas_from_manifest(path):
    rows = []

    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\\t")
        for row in reader:
            fasta = (
                row.get("final_mag_fasta")
                or row.get("refined_bin_fasta")
                or row.get("bin_fasta")
                or ""
            ).strip()

            mag_id = (
                row.get("final_mag_id")
                or row.get("refined_bin_id")
                or row.get("bin_id")
                or Path(fasta).stem
            )

            if fasta:
                rows.append((mag_id, Path(fasta)))

    return rows

def candidate_fastas_from_dir(path):
    rows = []

    if not path.exists():
        return rows

    patterns = [
        f"*.{mag_extension}",
        f"*.{mag_extension}.gz"
    ]

    for pattern in patterns:
        for p in sorted(path.glob(pattern)):
            if p.is_file():
                rows.append((p.stem.replace(".fa", "").replace(".fasta", ""), p))

    return rows

if input_mag_manifest:
    candidates = candidate_fastas_from_manifest(input_mag_manifest)
else:
    candidates = candidate_fastas_from_dir(input_mag_dir)

if not candidates:
    raise SystemExit(
        f"ERROR: No MAG FASTA files found. input_mag_dir={input_mag_dir}, input_mag_manifest={input_mag_manifest}"
    )

seen_names = set()
copied = 0
missing = 0
total_contigs = 0
total_bp = 0

with out_manifest.open("w") as manifest:
    print(
        "mag_id",
        "source_fasta",
        "prepared_fasta",
        "contig_count",
        "total_bp",
        sep="\\t",
        file=manifest
    )

    for idx, (mag_id, src) in enumerate(candidates, start=1):
        if not src.exists():
            log(f"WARNING: missing MAG FASTA skipped: {src}")
            missing += 1
            continue

        base = safe_id(mag_id or src.stem)
        name = f"{base}.fa"
        suffix = 1

        while name in seen_names:
            suffix += 1
            name = f"{base}_{suffix}.fa"

        seen_names.add(name)

        dest = out_dir / name

        if str(src).endswith(".gz"):
            with gzip.open(src, "rt", errors="replace") as inp, dest.open("w") as out:
                shutil.copyfileobj(inp, out)
        else:
            shutil.copyfile(src, dest)

        contigs, bp = fasta_stats(dest)

        copied += 1
        total_contigs += contigs
        total_bp += bp

        print(
            dest.stem,
            str(src),
            str(published_dir / name),
            contigs,
            bp,
            sep="\\t",
            file=manifest
        )

with out_stats.open("w") as stats:
    print(
        "input_mag_dir",
        "input_mag_manifest",
        "mags_prepared",
        "missing_mags",
        "total_contigs",
        "total_bp",
        "prepared_mag_dir",
        sep="\\t",
        file=stats
    )
    print(
        str(input_mag_dir),
        input_mag_manifest,
        copied,
        missing,
        total_contigs,
        total_bp,
        str(published_dir),
        sep="\\t",
        file=stats
    )

log(f"MAGs prepared: {copied}")
log(f"Missing MAGs: {missing}")
log(f"Total contigs: {total_contigs}")
log(f"Total bp: {total_bp}")
log("MAG input preparation finished.")
PY
    """
}


process SETUP_CHECKM2 {
    tag "setup_checkm2"

    publishDir "${params.outdir}/setup",
        mode: 'copy',
        pattern: "checkm2_setup_status.env"

    output:
    path "checkm2_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir ?: "${params.outdir}/conda_envs"
    def env_dir = "${base_env}/checkm2"

    def checkm2_pkg = params.checkm2_version ? "checkm2=${params.checkm2_version}" : "checkm2"

    """
    set -euo pipefail

    STATUS="checkm2_setup_status.env"
    CHECKM2_ENV="${env_dir}"
    CHECKM2_DB_PATH="${params.checkm2_db_path ?: ''}"
    CHECKM2_DB_DIR="${params.checkm2_db_outdir}"

    echo "CheckM2 setup started: \$(date)" > "\$STATUS"
    echo "CHECKM2_ENV=\$CHECKM2_ENV" >> "\$STATUS"
    echo "Requested CheckM2 package: ${checkm2_pkg}" >> "\$STATUS"
    echo "Initial CHECKM2_DB_PATH=\$CHECKM2_DB_PATH" >> "\$STATUS"
    echo "CHECKM2_DB_DIR=\$CHECKM2_DB_DIR" >> "\$STATUS"
    echo "----------------------------------------" >> "\$STATUS"

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

        if [[ ! -x "\$prefix/bin/checkm2" ]]; then
            return 1
        fi

        "\$prefix/bin/checkm2" --help >> "\$STATUS" 2>&1 || true
        return 0
    }

    if [[ -d "\$CHECKM2_ENV" ]]; then
        if check_env "\$CHECKM2_ENV"; then
            echo "Existing CheckM2 environment passed checks." >> "\$STATUS"
        else
            echo "Existing CheckM2 environment failed checks. Removing." >> "\$STATUS"
            rm -rf "\$CHECKM2_ENV"
        fi
    fi

    if [[ ! -d "\$CHECKM2_ENV" ]]; then
        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: CheckM2 environment missing and --auto_install false." >> "\$STATUS"
            exit 1
        fi

        INSTALLER="\$(find_installer)"

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda found." >> "\$STATUS"
            exit 1
        fi

        mkdir -p "\$(dirname "\$CHECKM2_ENV")"

        "\$INSTALLER" create -y \\
            -p "\$CHECKM2_ENV" \\
            -c bioconda \\
            -c conda-forge \\
            "${checkm2_pkg}" \\
            >> "\$STATUS" 2>&1

        if ! check_env "\$CHECKM2_ENV"; then
            echo "ERROR: Newly created CheckM2 environment failed checks." >> "\$STATUS"
            exit 1
        fi
    fi

    if [[ -z "\$CHECKM2_DB_PATH" ]]; then
        mkdir -p "\$CHECKM2_DB_DIR"

        EXISTING_DB="\$(find "\$CHECKM2_DB_DIR" -type f -name 'uniref100.KO.1.dmnd' 2>/dev/null | head -n 1 || true)"

        if [[ -n "\$EXISTING_DB" ]]; then
            CHECKM2_DB_PATH="\$EXISTING_DB"
        elif [[ "${params.checkm2_auto_download_db}" == "true" ]]; then
            echo "Downloading CheckM2 database from Zenodo record ${params.checkm2_zenodo_record}" >> "\$STATUS"

            python3 - "\$CHECKM2_DB_DIR" "${params.checkm2_zenodo_record}" "\$STATUS" <<'PY'
import json
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

db_dir = Path(sys.argv[1])
record = sys.argv[2]
status = Path(sys.argv[3])

db_dir.mkdir(parents=True, exist_ok=True)

def log(msg):
    with status.open("a") as handle:
        print(msg, file=handle)

api = f"https://zenodo.org/api/records/{record}"
log(f"Fetching Zenodo API record: {api}")

with urllib.request.urlopen(api) as response:
    data = json.load(response)

files = data.get("files", [])

candidates = []
for f in files:
    key = f.get("key", "")
    links = f.get("links", {})
    url = links.get("self") or links.get("download")
    if not url:
        continue
    if key.endswith(".tar.gz") or key.endswith(".tgz") or "checkm2" in key.lower() or "uniref100" in key.lower():
        candidates.append((key, url))

if not candidates:
    raise SystemExit("ERROR: Could not find a suitable CheckM2 database tarball in Zenodo record.")

key, url = candidates[0]
archive = db_dir / key

log(f"Downloading {key}")
log(f"URL: {url}")

with urllib.request.urlopen(url) as response, archive.open("wb") as out:
    shutil.copyfileobj(response, out)

log(f"Downloaded archive: {archive}")

if tarfile.is_tarfile(archive):
    log("Extracting archive.")
    with tarfile.open(archive) as tar:
        tar.extractall(db_dir)
else:
    log("Downloaded file was not a tar archive; leaving as-is.")

matches = list(db_dir.rglob("uniref100.KO.1.dmnd"))
if not matches:
    raise SystemExit("ERROR: Could not find uniref100.KO.1.dmnd after download/extraction.")

log(f"Found CheckM2 database: {matches[0]}")
PY

            CHECKM2_DB_PATH="\$(find "\$CHECKM2_DB_DIR" -type f -name 'uniref100.KO.1.dmnd' 2>/dev/null | head -n 1 || true)"
        else
            echo "ERROR: --checkm2_db_path not supplied and auto-download disabled." >> "\$STATUS"
            exit 1
        fi
    fi

    if [[ ! -s "\$CHECKM2_DB_PATH" ]]; then
        echo "ERROR: CheckM2 database path does not exist or is empty: \$CHECKM2_DB_PATH" >> "\$STATUS"
        exit 1
    fi

    echo "CHECKM2_ENV=\$CHECKM2_ENV" >> "\$STATUS"
    echo "CHECKM2_DB_PATH=\$CHECKM2_DB_PATH" >> "\$STATUS"
    echo "CheckM2 setup finished: \$(date)" >> "\$STATUS"
    """
}


process RUN_CHECKM2 {
    tag "checkm2"

    publishDir "${params.outdir}/checkm2",
        mode: params.publish_tool_outputs_mode,
        pattern: "checkm2_out/**",
        saveAs: { filename -> filename.replaceFirst(/^checkm2_out\//, '') }

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "checkm2.log"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "checkm2_status.tsv"

    cpus {
        params.threads != null ? params.threads as int : 8
    }

    input:
    path mags_dir
    path setup_status

    output:
    path "checkm2_status.tsv", emit: status
    path "checkm2.log", emit: log_file
    path "checkm2_out", emit: checkm2_out

    script:
    """
    set -euo pipefail

    LOG="checkm2.log"

    CHECKM2_ENV="\$(grep '^CHECKM2_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    CHECKM2_DB_PATH="\$(grep '^CHECKM2_DB_PATH=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"

    export PATH="\$CHECKM2_ENV/bin:\$PATH"

    echo "CheckM2 started: \$(date)" > "\$LOG"
    echo "MAG directory: ${mags_dir}" >> "\$LOG"
    echo "CheckM2 DB: \$CHECKM2_DB_PATH" >> "\$LOG"
    echo "Threads: ${task.cpus}" >> "\$LOG"

    mkdir -p checkm2_out

    set +e
    checkm2 predict \\
        --threads ${task.cpus} \\
        --input "${mags_dir}" \\
        -x ${params.checkm2_extension} \\
        --output-directory checkm2_out \\
        --database_path "\$CHECKM2_DB_PATH" \\
        >> "\$LOG" 2>&1
    STATUS="\$?"
    set -e

    printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > checkm2_status.tsv

    if [[ "\$STATUS" -ne 0 ]]; then
        printf 'checkm2\\tfailed\\t%s\\t%s\\tCheckM2 failed\\n' "\$STATUS" "${params.outdir}/checkm2" >> checkm2_status.tsv
        exit "\$STATUS"
    fi

    printf 'checkm2\\tcompleted\\t0\\t%s\\tCheckM2 completed\\n' "${params.outdir}/checkm2" >> checkm2_status.tsv

    echo "CheckM2 finished: \$(date)" >> "\$LOG"
    """
}


process SETUP_GTDBTK {
    tag "setup_gtdbtk"

    publishDir "${params.outdir}/setup",
        mode: 'copy',
        pattern: "gtdbtk_setup_status.env"

    output:
    path "gtdbtk_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir ?: "${params.outdir}/conda_envs"
    def env_dir = "${base_env}/gtdbtk-${params.gtdbtk_version}"

    """
    set -euo pipefail

    STATUS="gtdbtk_setup_status.env"
    GTDBTK_ENV="${env_dir}"
    GTDBTK_DATA_PATH="${params.gtdbtk_data_path ?: ''}"
    GTDBTK_DB_DIR="${params.gtdbtk_db_outdir}"

    echo "GTDB-Tk setup started: \$(date)" > "\$STATUS"
    echo "GTDBTK_ENV=\$GTDBTK_ENV" >> "\$STATUS"
    echo "Initial GTDBTK_DATA_PATH=\$GTDBTK_DATA_PATH" >> "\$STATUS"
    echo "GTDBTK_DB_DIR=\$GTDBTK_DB_DIR" >> "\$STATUS"

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
        [[ -x "\$prefix/bin/gtdbtk" ]]
    }

    if [[ -d "\$GTDBTK_ENV" ]]; then
        if check_env "\$GTDBTK_ENV"; then
            echo "Existing GTDB-Tk environment passed checks." >> "\$STATUS"
        else
            echo "Existing GTDB-Tk environment failed checks. Removing." >> "\$STATUS"
            rm -rf "\$GTDBTK_ENV"
        fi
    fi

    if [[ ! -d "\$GTDBTK_ENV" ]]; then
        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: GTDB-Tk environment missing and auto_install false." >> "\$STATUS"
            exit 1
        fi

        INSTALLER="\$(find_installer)"

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda found." >> "\$STATUS"
            exit 1
        fi

        mkdir -p "\$(dirname "\$GTDBTK_ENV")"

        "\$INSTALLER" create -y \\
            -p "\$GTDBTK_ENV" \\
            -c conda-forge \\
            -c bioconda \\
            "gtdbtk=${params.gtdbtk_version}" \\
            >> "\$STATUS" 2>&1
    fi

    if ! check_env "\$GTDBTK_ENV"; then
        echo "ERROR: GTDB-Tk environment failed checks." >> "\$STATUS"
        exit 1
    fi

    export PATH="\$GTDBTK_ENV/bin:\$PATH"

    if [[ -z "\$GTDBTK_DATA_PATH" ]]; then
        mkdir -p "\$GTDBTK_DB_DIR"

        EXISTING_DATA="\$(find "\$GTDBTK_DB_DIR" -type f -name 'metadata.txt' -o -type f -name 'VERSION' 2>/dev/null | head -n 1 || true)"

        if [[ -n "\$EXISTING_DATA" ]]; then
            GTDBTK_DATA_PATH="\$(dirname "\$EXISTING_DATA")"
        elif [[ "${params.gtdbtk_auto_download_db}" == "true" ]]; then
            echo "Attempting GTDB-Tk database download with download-db.sh" >> "\$STATUS"

            cd "\$GTDBTK_DB_DIR"

            set +e
            download-db.sh >> "\$OLDPWD/\$STATUS" 2>&1
            DL_STATUS="\$?"
            set -e

            cd "\$OLDPWD"

            if [[ "\$DL_STATUS" -ne 0 ]]; then
                echo "ERROR: download-db.sh failed with status \$DL_STATUS" >> "\$STATUS"
                echo "You may need to manually download GTDB-Tk data and rerun with --gtdbtk_data_path." >> "\$STATUS"
                exit "\$DL_STATUS"
            fi

            GTDBTK_DATA_PATH="\$(find "\$GTDBTK_DB_DIR" -type f -name 'metadata.txt' -o -type f -name 'VERSION' 2>/dev/null | head -n 1 | xargs dirname || true)"
        else
            echo "ERROR: --gtdbtk_data_path not supplied and auto-download disabled." >> "\$STATUS"
            exit 1
        fi
    fi

    if [[ ! -d "\$GTDBTK_DATA_PATH" ]]; then
        echo "ERROR: GTDBTK_DATA_PATH does not exist: \$GTDBTK_DATA_PATH" >> "\$STATUS"
        exit 1
    fi

    echo "GTDBTK_ENV=\$GTDBTK_ENV" >> "\$STATUS"
    echo "GTDBTK_DATA_PATH=\$GTDBTK_DATA_PATH" >> "\$STATUS"
    echo "GTDB-Tk setup finished: \$(date)" >> "\$STATUS"
    """
}


process RUN_GTDBTK {
    tag "gtdbtk"

    publishDir "${params.outdir}/gtdbtk",
        mode: params.publish_tool_outputs_mode,
        pattern: "gtdbtk_out/**",
        saveAs: { filename -> filename.replaceFirst(/^gtdbtk_out\//, '') }

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "gtdbtk.log"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "gtdbtk_status.tsv"

    cpus {
        params.threads != null ? params.threads as int : 32
    }

    input:
    path mags_dir
    path setup_status

    output:
    path "gtdbtk_status.tsv", emit: status
    path "gtdbtk.log", emit: log_file
    path "gtdbtk_out", emit: gtdbtk_out

    script:
    """
    set -euo pipefail

    LOG="gtdbtk.log"

    GTDBTK_ENV="\$(grep '^GTDBTK_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    GTDBTK_DATA_PATH="\$(grep '^GTDBTK_DATA_PATH=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"

    export PATH="\$GTDBTK_ENV/bin:\$PATH"
    export GTDBTK_DATA_PATH="\$GTDBTK_DATA_PATH"

    echo "GTDB-Tk started: \$(date)" > "\$LOG"
    echo "MAG directory: ${mags_dir}" >> "\$LOG"
    echo "GTDBTK_DATA_PATH: \$GTDBTK_DATA_PATH" >> "\$LOG"
    echo "Threads: ${task.cpus}" >> "\$LOG"

    mkdir -p gtdbtk_out gtdbtk_mash

    set +e
    gtdbtk classify_wf \\
        --genome_dir "${mags_dir}" \\
        --out_dir gtdbtk_out \\
        --mash_db gtdbtk_mash \\
        --extension ${params.gtdbtk_extension} \\
        --cpus ${task.cpus} \\
        >> "\$LOG" 2>&1
    STATUS="\$?"
    set -e

    printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > gtdbtk_status.tsv

    if [[ "\$STATUS" -ne 0 ]]; then
        printf 'gtdbtk\\tfailed\\t%s\\t%s\\tGTDB-Tk failed\\n' "\$STATUS" "${params.outdir}/gtdbtk" >> gtdbtk_status.tsv
        exit "\$STATUS"
    fi

    printf 'gtdbtk\\tcompleted\\t0\\t%s\\tGTDB-Tk completed\\n' "${params.outdir}/gtdbtk" >> gtdbtk_status.tsv

    echo "GTDB-Tk finished: \$(date)" >> "\$LOG"
    """
}


process RUN_DRAM2 {
    tag "dram2"

    publishDir "${params.outdir}/dram2",
        mode: params.publish_tool_outputs_mode,
        pattern: "DRAM2/**",
        saveAs: { filename -> filename.replaceFirst(/^DRAM2\//, '') }

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "dram2.log"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "dram2_status.tsv"

    cpus {
        params.threads != null ? params.threads as int : 36
    }

    input:
    path mags_dir

    output:
    path "dram2_status.tsv", emit: status
    path "dram2.log", emit: log_file
    path "DRAM2", emit: dram2_out

    script:
    """
    set -euo pipefail

    LOG="dram2.log"

    echo "DRAM2 started: \$(date)" > "\$LOG"
    echo "MAG directory: ${mags_dir}" >> "\$LOG"
    echo "DRAM repo: ${params.dram2_repo}" >> "\$LOG"
    echo "DRAM revision: ${params.dram2_revision}" >> "\$LOG"
    echo "Profile: ${params.dram2_profile}" >> "\$LOG"

    mkdir -p DRAM_setup DRAM2

    cd DRAM_setup

    curl -L -o nextflow.config "${params.dram2_nextflow_config_url}" >> "../\$LOG" 2>&1

    ${params.nextflow_exe} pull "${params.dram2_repo}" -r "${params.dram2_revision}" >> "../\$LOG" 2>&1 || true

    cd ..

    set +e
    ${params.nextflow_exe} run "${params.dram2_repo}" \\
        -r "${params.dram2_revision}" \\
        -c DRAM_setup/nextflow.config \\
        --input_fasta "${mags_dir}" \\
        --outdir DRAM2 \\
        --threads ${task.cpus} \\
        --rename --annotate --qc --summarize --visualize --traits \\
        --use_camper --use_canthyd --use_dbcan --use_fegenie --use_kofam --use_dram_db --use_antismash --use_tcdb --use_vog --use_merops --use_methyl --use_sulfur --use_metals \\
        --tiny_cpus_limit ${params.dram2_tiny_cpus_limit} \\
        --small_cpus_limit ${params.dram2_small_cpus_limit} \\
        --medium_cpus_limit ${params.dram2_medium_cpus_limit} \\
        --big_cpus_limit ${params.dram2_big_cpus_limit} \\
        --huge_cpus_limit ${params.dram2_huge_cpus_limit} \\
        --tiny_gb_mem_limit ${params.dram2_tiny_gb_mem_limit} \\
        --small_gb_mem_limit ${params.dram2_small_gb_mem_limit} \\
        --medium_gb_mem_limit ${params.dram2_medium_gb_mem_limit} \\
        --big_gb_mem_limit ${params.dram2_big_gb_mem_limit} \\
        --huge_gb_mem_limit ${params.dram2_huge_gb_mem_limit} \\
        --tiny_hr_time_limit ${params.dram2_tiny_hr_time_limit} \\
        --small_hr_time_limit ${params.dram2_small_hr_time_limit} \\
        --medium_hr_time_limit ${params.dram2_medium_hr_time_limit} \\
        --big_hr_time_limit ${params.dram2_big_hr_time_limit} \\
        --huge_hr_time_limit ${params.dram2_huge_hr_time_limit} \\
        -profile "${params.dram2_profile}" \\
        >> "\$LOG" 2>&1
    STATUS="\$?"
    set -e

    printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > dram2_status.tsv

    if [[ "\$STATUS" -ne 0 ]]; then
        printf 'dram2\\tfailed\\t%s\\t%s\\tDRAM2 failed\\n' "\$STATUS" "${params.outdir}/dram2" >> dram2_status.tsv
        exit "\$STATUS"
    fi

    printf 'dram2\\tcompleted\\t0\\t%s\\tDRAM2 completed\\n' "${params.outdir}/dram2" >> dram2_status.tsv

    echo "DRAM2 finished: \$(date)" >> "\$LOG"
    """
}


process SETUP_MICROTRAIT {
    tag "setup_microtrait"

    publishDir "${params.outdir}/setup",
        mode: 'copy',
        pattern: "microtrait_setup_status.env"

    output:
    path "microtrait_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir ?: "${params.outdir}/conda_envs"
    def env_dir = params.microtrait_env_dir ?: "${base_env}/microtrait"

    """
    set -euo pipefail

    STATUS="microtrait_setup_status.env"
    MICROTRAIT_ENV="${env_dir}"

    echo "microTrait setup started: \$(date)" > "\$STATUS"
    echo "MICROTRAIT_ENV=\$MICROTRAIT_ENV" >> "\$STATUS"

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

        if [[ ! -x "\$prefix/bin/Rscript" ]]; then
            return 1
        fi

        if [[ ! -x "\$prefix/bin/hmmsearch" ]]; then
            return 1
        fi

        if [[ ! -x "\$prefix/bin/prodigal" ]]; then
            return 1
        fi

        "\$prefix/bin/Rscript" -e 'library(microtrait)' >> "\$STATUS" 2>&1
    }

    if [[ -d "\$MICROTRAIT_ENV" ]]; then
        if check_env "\$MICROTRAIT_ENV"; then
            echo "Existing microTrait environment passed checks." >> "\$STATUS"
            echo "MICROTRAIT_ENV=\$MICROTRAIT_ENV" >> "\$STATUS"
            echo "microTrait setup finished: \$(date)" >> "\$STATUS"
            exit 0
        else
            echo "Existing microTrait environment failed checks." >> "\$STATUS"

            if [[ -z "${params.microtrait_env_dir ?: ''}" ]]; then
                echo "Removing workflow-managed microTrait environment." >> "\$STATUS"
                rm -rf "\$MICROTRAIT_ENV"
            else
                echo "ERROR: User-supplied microTrait env failed checks." >> "\$STATUS"
                exit 1
            fi
        fi
    fi

    if [[ ! -d "\$MICROTRAIT_ENV" ]]; then
        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: microTrait environment missing and auto_install false." >> "\$STATUS"
            exit 1
        fi

        INSTALLER="\$(find_installer)"

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda found." >> "\$STATUS"
            exit 1
        fi

        mkdir -p "\$(dirname "\$MICROTRAIT_ENV")"

        "\$INSTALLER" create -y \\
            -p "\$MICROTRAIT_ENV" \\
            -c conda-forge \\
            -c bioconda \\
            "r-base" \\
            "r-devtools" \\
            "r-biocmanager" \\
            "r-dplyr" \\
            "r-readr" \\
            "r-stringr" \\
            "r-tidyr" \\
            "r-tibble" \\
            "r-tictoc" \\
            "r-ape" \\
            "r-assertthat" \\
            "r-checkmate" \\
            "r-r.utils" \\
            "r-rcolorbrewer" \\
            "r-corrplot" \\
            "r-doparallel" \\
            "r-futile.logger" \\
            "r-gtools" \\
            "r-lazyeval" \\
            "r-magrittr" \\
            "r-pheatmap" \\
            "hmmer" \\
            "prodigal" \\
            "infernal" \\
            "trnascan-se" \\
            "bedtools" \\
            >> "\$STATUS" 2>&1

        export PATH="\$MICROTRAIT_ENV/bin:\$PATH"

        Rscript - <<'RSCRIPT' >> "\$STATUS" 2>&1
options(repos = c(CRAN = "https://cloud.r-project.org"))

cran_pkgs <- c(
  "R.utils", "RColorBrewer", "ape", "assertthat", "checkmate",
  "corrplot", "doParallel", "dplyr", "futile.logger", "gtools",
  "lazyeval", "magrittr", "parallel", "pheatmap", "readr",
  "stringr", "tibble", "tictoc", "tidyr", "devtools"
)

missing <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
if(length(missing) > 0) {
  install.packages(missing)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("Biostrings", "coRdon", "ComplexHeatmap")
for(pkg in bioc_pkgs) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

if(!requireNamespace("gRodon", quietly = TRUE)) {
  devtools::install_github("jlw-ecoevo/gRodon", upgrade = "never")
}

if(!requireNamespace("microtrait", quietly = TRUE)) {
  devtools::install_github("ukaraoz/microtrait", upgrade = "never")
}

library(microtrait)
try(microtrait::prep.hmmmodels(), silent = TRUE)
RSCRIPT
    fi

    if ! check_env "\$MICROTRAIT_ENV"; then
        echo "ERROR: microTrait environment failed final checks." >> "\$STATUS"
        exit 1
    fi

    echo "MICROTRAIT_ENV=\$MICROTRAIT_ENV" >> "\$STATUS"
    echo "microTrait setup finished: \$(date)" >> "\$STATUS"
    """
}


process RUN_MICROTRAIT {
    tag "microtrait"

    publishDir "${params.outdir}/microtrait",
        mode: params.publish_tool_outputs_mode,
        pattern: "microtrait_out/**",
        saveAs: { filename -> filename.replaceFirst(/^microtrait_out\//, '') }

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "microtrait.log"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "microtrait_status.tsv"

    cpus {
        params.microtrait_cores != null
            ? params.microtrait_cores as int
            : (params.threads != null ? params.threads as int : 10)
    }

    input:
    path mags_dir
    path setup_status

    output:
    path "microtrait_status.tsv", emit: status
    path "microtrait.log", emit: log_file
    path "microtrait_out", emit: microtrait_out

    script:
    """
    set -euo pipefail

    LOG="microtrait.log"

    MICROTRAIT_ENV="\$(grep '^MICROTRAIT_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"

    export PATH="\$MICROTRAIT_ENV/bin:\$PATH"

    mkdir -p microtrait_out microtrait_work

    echo "microTrait started: \$(date)" > "\$LOG"
    echo "MAG directory: ${mags_dir}" >> "\$LOG"
    echo "Run type: ${params.microtrait_run_type}" >> "\$LOG"
    echo "Output name: ${params.microtrait_out_name}" >> "\$LOG"
    echo "Threads: ${task.cpus}" >> "\$LOG"

    if [[ -n "${params.microtrait_runner ?: ''}" ]]; then
        cp -L "${params.microtrait_runner}" microtrait_runner.R
    else
        cat > microtrait_runner.R <<'RSCRIPT'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

fasta_dir <- args[1]
out_name <- args[2]
type <- args[3]
source_dir <- args[4]
ncores <- as.integer(args[5])

library(microtrait)
library(dplyr)

if(!is.na(source_dir) && source_dir != "" && dir.exists(source_dir)) {
  source(file.path(source_dir, "ogt.R"))
  source(file.path(source_dir, "protein.R"))
  source(file.path(source_dir, "nucleotide.R"))
  source(file.path(source_dir, "extern.R"))
  source(file.path(source_dir, "mingentime.R"))
  source(file.path(source_dir, "utils.R"))
}

if(type == "genomic") {
  fasta.list <- list.files(path = fasta_dir, pattern = ".fa$", full.names = TRUE)

  microtrait_results <- extract.traits.parallel(fa_files = fasta.list)

  rds_files <- unlist(parallel::mclapply(
    microtrait_results,
    "[[",
    "rds_file",
    mc.cores = min(10, ncores)
  ))

  genomeset_results <- make.genomeset.results(
    rds_files = rds_files,
    ids = sub(".microtrait.rds", "", basename(rds_files)),
    ncores = 1
  )

  saveRDS(
    genomeset_results,
    file.path(fasta_dir, paste0(out_name, "_GenomeSet_Traits.rds"))
  )

  cat(paste("Finished extracting traits from", length(fasta.list), "genomes.\\n"))

} else if(type == "protein") {
  fasta.list <- list.files(path = fasta_dir, pattern = ".fa$", full.names = TRUE)
  protein.files <- list.files(path = fasta_dir, pattern = ".faa$", full.names = TRUE)

  microtrait_results <- parallel::mclapply(
    seq_along(protein.files),
    function(i) {
      extract.traits(
        protein.files[i],
        fasta_dir,
        type = "protein",
        growthrate_predict = FALSE,
        optimalT_predict = FALSE
      )
    },
    mc.cores = max(1, floor(ncores * 0.7))
  )

  rds_files <- unlist(parallel::mclapply(
    microtrait_results,
    "[[",
    "rds_file",
    mc.cores = min(10, ncores)
  ))

  genomeset_results <- make.genomeset.results(
    rds_files = rds_files,
    ids = sub(".microtrait.rds", "", basename(rds_files)),
    growthrate = FALSE,
    optimumT = FALSE,
    ncores = 1
  )

  saveRDS(
    genomeset_results,
    file.path(fasta_dir, paste0(out_name, "_GenomeSet_Traits.rds"))
  )

} else {
  stop("type must be genomic or protein")
}
RSCRIPT
    fi

    set +e
    Rscript microtrait_runner.R \\
        "${mags_dir}" \\
        "${params.microtrait_out_name}" \\
        "${params.microtrait_run_type}" \\
        "${params.microtrait_source_dir ?: ''}" \\
        "${task.cpus}" \\
        >> "\$LOG" 2>&1
    STATUS="\$?"
    set -e

    find "${mags_dir}" -maxdepth 1 -type f \\( -name '*.rds' -o -name '*.csv' \\) -exec cp -L {} microtrait_out/ \\; 2>> "\$LOG" || true

    printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > microtrait_status.tsv

    if [[ "\$STATUS" -ne 0 ]]; then
        printf 'microtrait\\tfailed\\t%s\\t%s\\tmicroTrait failed\\n' "\$STATUS" "${params.outdir}/microtrait" >> microtrait_status.tsv
        exit "\$STATUS"
    fi

    printf 'microtrait\\tcompleted\\t0\\t%s\\tmicroTrait completed\\n' "${params.outdir}/microtrait" >> microtrait_status.tsv

    echo "microTrait finished: \$(date)" >> "\$LOG"
    """
}


process WRITE_MODULE6_SUMMARY {
    tag "write_module6_summary"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "module6_run_summary.tsv"

    input:
    path input_stats
    path tool_status_files

    output:
    path "module6_run_summary.tsv", emit: summary

    script:
    def status_files = tool_status_files.collect { status_file ->
        status_file.name
    }.join(' ')

    """
    set -euo pipefail

    printf 'section\\tsource_file\\n' > module6_run_summary.tsv
    printf 'mag_inputs\\t%s\\n' "${input_stats}" >> module6_run_summary.tsv

    if [[ -n "${status_files}" ]]; then
        for f in ${status_files}; do
            printf 'tool_status\\t%s\\n' "\$f" >> module6_run_summary.tsv
        done
    fi

    echo "" >> module6_run_summary.tsv
    echo "# module6_mag_input_stats.tsv" >> module6_run_summary.tsv
    cat "${input_stats}" >> module6_run_summary.tsv

    if [[ -n "${status_files}" ]]; then
        for f in ${status_files}; do
            echo "" >> module6_run_summary.tsv
            echo "# \$f" >> module6_run_summary.tsv
            cat "\$f" >> module6_run_summary.tsv
        done
    fi
    """
}