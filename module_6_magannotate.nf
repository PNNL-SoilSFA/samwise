#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Module 6: MAG annotation / quality / taxonomy.
 *
 * Tools:
 *   - CheckM2
 *   - GTDB-Tk
 *   - EggNOG-mapper
 *
 * DRAM2 and microTrait are intentionally removed for now.
 */

/*
 * Core parameters
 */
params.working_dir = null
params.output_dir = null
params.input_mag_dir = null
params.input_mag_manifest = null
params.mag_extension = "fa"
params.threads = null
params.tool_env_dir = null
params.auto_install = true

/*
 * Workflow-local conda/mamba package cache.
 */
params.conda_pkgs_dir = null

/*
 * Which tools to run
 */
params.run_checkm2 = true
params.run_gtdbtk = true
params.run_eggnog = true

/*
 * Tool versions
 */
params.checkm2_version = null
params.gtdbtk_version = "2.7.2"

/*
 * CheckM2 database.
 *
 * If --checkm2_db_path is supplied, use that file directly.
 * Else check/download into:
 *   <outdir>/databases/checkm2
 */
params.checkm2_db_path = null
params.checkm2_db_dir = null
params.checkm2_zenodo_record = "14897628"
params.checkm2_auto_download_db = true
params.checkm2_extension = "fa"

/*
 * GTDB-Tk database.
 *
 * If --gtdbtk_data_path is supplied, use that directory directly.
 * Else check/download into:
 *   <outdir>/databases/gtdbtk
 */
params.gtdbtk_data_path = null
params.gtdbtk_db_dir = null
params.gtdbtk_auto_download_db = true
params.gtdbtk_extension = "fa"
params.gtdbtk_download_url = "https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"

/*
 * EggNOG-mapper.
 *
 * If --eggnog_data_dir is supplied, use/check/download into that directory.
 * Else check/download into:
 *   <outdir>/databases/eggnog
 */
params.eggnog_env_dir = null
params.eggnog_mapper_version = "2.1.13"
params.eggnog_data_path = null
params.eggnog_data_dir = null
params.eggnog_auto_download_db = true
params.eggnog_download_args = "-y"
params.eggnog_fixurl = true
params.eggnog_fixurl_package = "eggnog-mapper-fixurl"
params.eggnog_method = "diamond"
params.eggnog_itype = "metagenome"
params.eggnog_genepred = "prodigal"
params.eggnog_trans_table = 11
params.eggnog_output_prefix = "samwise_eggnog"
params.eggnog_extra_args = ""
params.eggnog_mmseqs_db = null
params.eggnog_fail_nonfatal = false

/*
 * Publish modes
 */
params.publish_mags_mode = "copy"
params.publish_tool_outputs_mode = "copy"

/*
 * Derived directories
 */
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.outdir = "${params.results_dir}/module_6_magannotation"

/*
 * Candidate MAG locations
 */
params.module5_final_mag_dir = "${params.results_dir}/module_5_subtractiveassembly/final_mag_database"
params.module4_refined_mag_dir = "${params.results_dir}/module_4_binrefinement/refined_bins"

params.checkm2_db_outdir = params.checkm2_db_dir ?: "${params.outdir}/databases/checkm2"
params.gtdbtk_db_outdir = params.gtdbtk_db_dir ?: "${params.outdir}/databases/gtdbtk"
params.eggnog_db_outdir = params.eggnog_data_path ?: (params.eggnog_data_dir ?: "${params.outdir}/databases/eggnog")


def absPath(value) {
    def s = value == null ? "" : value.toString().trim()
    if (!s || s == "null" || s == "NA") {
        return ""
    }
    return java.nio.file.Paths.get(s).toAbsolutePath().normalize().toString()
}


def firstExistingMagDir(candidates) {
    def suffixes = [".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"]

    def found_with_fastas = candidates.find { candidate ->
        def d = new File(candidate.toString())
        if (!d.exists() || !d.isDirectory()) {
            return false
        }
        def files = d.listFiles()
        if (files == null) {
            return false
        }
        return files.any { f ->
            f.isFile() && suffixes.any { suffix -> f.getName().endsWith(suffix) }
        }
    }

    def found_dir = candidates.find { candidate ->
        def d = new File(candidate.toString())
        d.exists() && d.isDirectory()
    }

    def selected = found_with_fastas ?: found_dir ?: candidates[0]
    return java.nio.file.Paths.get(selected.toString()).toAbsolutePath().normalize().toString()
}


workflow {
    def run_checkm2 = params.run_checkm2.toString().toBoolean()
    def run_gtdbtk = params.run_gtdbtk.toString().toBoolean()
    def run_eggnog = params.run_eggnog.toString().toBoolean()

    if (!run_checkm2 && !run_gtdbtk && !run_eggnog) {
        error(
            """
            No Module 6 tools selected.

            Enable at least one of:
              --run_checkm2 true
              --run_gtdbtk true
              --run_eggnog true
            """.stripIndent()
        )
    }

    def selected_mag_dir = params.input_mag_dir
        ? absPath(params.input_mag_dir)
        : firstExistingMagDir(
            [params.module5_final_mag_dir, params.module4_refined_mag_dir]
        )

    def selected_mag_manifest = params.input_mag_manifest
        ? absPath(params.input_mag_manifest)
        : ""

    log.info("Module 6 results directory: ${params.results_dir}")
    log.info("Writing Module 6 outputs to: ${params.outdir}")
    log.info("Selected MAG directory: ${selected_mag_dir}")
    log.info("Selected MAG manifest: ${selected_mag_manifest ?: 'not supplied'}")
    log.info("MAG extension: ${params.mag_extension}")
    log.info("Threads: ${params.threads ?: 'tool-specific defaults'}")
    log.info("Run CheckM2: ${run_checkm2}")
    log.info("Run GTDB-Tk: ${run_gtdbtk}")
    log.info("Run EggNOG-mapper: ${run_eggnog}")
    log.info("Conda package cache root: ${params.conda_pkgs_dir ?: params.outdir + '/conda_pkgs'}")
    log.info("Per-tool conda package caches will be used under the cache root.")

    PREPARE_MAG_INPUTS(
        channel.value(selected_mag_dir),
        channel.value(selected_mag_manifest),
    )

    if (run_checkm2) {
        SETUP_CHECKM2()
        RUN_CHECKM2(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_CHECKM2.out.status,
        )
    }

    if (run_gtdbtk) {
        SETUP_GTDBTK()
        RUN_GTDBTK(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_GTDBTK.out.status,
        )
    }

    if (run_eggnog) {
        SETUP_EGGNOG()
        RUN_EGGNOG(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_EGGNOG.out.status,
        )
    }

    def status_ch

    if (run_checkm2 && run_gtdbtk && run_eggnog) {
        status_ch = RUN_CHECKM2.out.status.mix(RUN_GTDBTK.out.status).mix(RUN_EGGNOG.out.status)
    }
    else if (run_checkm2 && run_gtdbtk && !run_eggnog) {
        status_ch = RUN_CHECKM2.out.status.mix(RUN_GTDBTK.out.status)
    }
    else if (run_checkm2 && !run_gtdbtk && run_eggnog) {
        status_ch = RUN_CHECKM2.out.status.mix(RUN_EGGNOG.out.status)
    }
    else if (!run_checkm2 && run_gtdbtk && run_eggnog) {
        status_ch = RUN_GTDBTK.out.status.mix(RUN_EGGNOG.out.status)
    }
    else if (run_checkm2 && !run_gtdbtk && !run_eggnog) {
        status_ch = RUN_CHECKM2.out.status
    }
    else if (!run_checkm2 && run_gtdbtk && !run_eggnog) {
        status_ch = RUN_GTDBTK.out.status
    }
    else if (!run_checkm2 && !run_gtdbtk && run_eggnog) {
        status_ch = RUN_EGGNOG.out.status
    }

    WRITE_MODULE6_SUMMARY(
        PREPARE_MAG_INPUTS.out.input_stats,
        status_ch.collect(),
    )
}


process PREPARE_MAG_INPUTS {
    tag "prepare_final_mags"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "module6_mag_input_*.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "prepare_mag_inputs.log"

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
    echo "Task working directory: \$(pwd -P)" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    mkdir -p refined_genomes

    python3 - \
    "${input_mag_dir}" \
    "${input_mag_manifest}" \
    "${params.mag_extension}" \
    "refined_genomes" \
    "module6_mag_input_manifest.tsv" \
    "module6_mag_input_stats.tsv" \
    "\$LOG_FILE" \
    "\$(pwd -P)/refined_genomes" <<'PY'
    
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
    manifest = Path(path)

    if not manifest.exists():
        raise SystemExit(f"ERROR: input MAG manifest does not exist: {manifest}")

    with manifest.open() as handle:
        reader = csv.DictReader(handle, delimiter="\\t")
        for row in reader:
            fasta = (
                row.get("final_mag_fasta")
                or row.get("refined_bin_fasta")
                or row.get("bin_fasta")
                or row.get("prepared_fasta")
                or ""
            ).strip()

            mag_id = (
                row.get("final_mag_id")
                or row.get("refined_bin_id")
                or row.get("bin_id")
                or row.get("mag_id")
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
        f"*.{mag_extension}.gz",
        "*.fa",
        "*.fa.gz",
        "*.fna",
        "*.fna.gz",
        "*.fasta",
        "*.fasta.gz",
    ]

    seen = set()

    for pattern in patterns:
        for p in sorted(path.glob(pattern)):
            if p.is_file() and p not in seen:
                seen.add(p)
                stem = p.name

                for suffix in [".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"]:
                    if stem.endswith(suffix):
                        stem = stem[:-len(suffix)]
                        break

                rows.append((stem, p))

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
zero_contig = 0
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

        if contigs == 0:
            log(f"WARNING: prepared MAG has zero contigs: {dest}")
            zero_contig += 1

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
        "zero_contig_mags",
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
        zero_contig,
        total_contigs,
        total_bp,
        str(published_dir),
        sep="\\t",
        file=stats
    )

log(f"MAGs prepared: {copied}")
log(f"Missing MAGs: {missing}")
log(f"Zero-contig MAGs: {zero_contig}")
log(f"Total contigs: {total_contigs}")
log(f"Total bp: {total_bp}")
log("MAG input preparation finished.")
PY

    echo "Preparing MAG inputs finished: \$(date)" >> "\$LOG_FILE"
    """
}


process SETUP_CHECKM2 {
    tag "setup_checkm2"

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "checkm2_setup_status.env"

    output:
    path "checkm2_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir ? absPath(params.tool_env_dir) : "${absPath(params.outdir)}/conda_envs"
    def env_dir = "${base_env}/checkm2"
    def checkm2_pkg = params.checkm2_version ? "checkm2=${params.checkm2_version}" : "checkm2"
    def db_path = params.checkm2_db_path ? absPath(params.checkm2_db_path) : ""
    def db_dir = absPath(params.checkm2_db_outdir)
    def conda_pkgs_dir = params.conda_pkgs_dir ? "${absPath(params.conda_pkgs_dir)}/checkm2" : "${absPath(params.outdir)}/conda_pkgs/checkm2"

    """
    set -euo pipefail

    STATUS="checkm2_setup_status.env"
    CHECKM2_ENV="${env_dir}"
    CHECKM2_DB_PATH="${db_path}"
    CHECKM2_DB_DIR="${db_dir}"

    CONDA_PKGS_DIRS="${conda_pkgs_dir}"
    export CONDA_PKGS_DIRS
    mkdir -p "\$CONDA_PKGS_DIRS"

    echo "CheckM2 setup started: \$(date)" > "\$STATUS"
    echo "CHECKM2_ENV=\$CHECKM2_ENV" >> "\$STATUS"
    echo "CHECKM2_DB_PATH=\$CHECKM2_DB_PATH" >> "\$STATUS"
    echo "CHECKM2_DB_DIR=\$CHECKM2_DB_DIR" >> "\$STATUS"
    echo "CONDA_PKGS_DIRS=\$CONDA_PKGS_DIRS" >> "\$STATUS"
    echo "Requested CheckM2 package: ${checkm2_pkg}" >> "\$STATUS"
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
        [[ -x "\$prefix/bin/checkm2" ]] || return 1
        [[ -x "\$prefix/bin/python" ]] || return 1
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
            --override-channels \\
            -c bioconda \\
            -c conda-forge \\
            "${checkm2_pkg}" \\
            "wget" \\
            "curl" \\
            "tar" \\
            "gzip" \\
            "pigz" \\
            >> "\$STATUS" 2>&1
    fi

    if ! check_env "\$CHECKM2_ENV"; then
        echo "ERROR: CheckM2 environment failed checks." >> "\$STATUS"
        exit 1
    fi

    export PATH="\$CHECKM2_ENV/bin:\$PATH"
    mkdir -p "\$CHECKM2_DB_DIR"

    if [[ -n "\$CHECKM2_DB_PATH" ]]; then
        echo "Using user-supplied CheckM2 database path: \$CHECKM2_DB_PATH" >> "\$STATUS"

        if [[ ! -s "\$CHECKM2_DB_PATH" ]]; then
            echo "ERROR: User-supplied CheckM2 DB path missing/empty: \$CHECKM2_DB_PATH" >> "\$STATUS"
            exit 1
        fi

    else
        EXISTING_DB="\$(find "\$CHECKM2_DB_DIR" -type f -name 'uniref100.KO.1.dmnd' 2>/dev/null | head -n 1 || true)"

        if [[ -n "\$EXISTING_DB" ]]; then
            CHECKM2_DB_PATH="\$EXISTING_DB"
            echo "Existing CheckM2 database found: \$CHECKM2_DB_PATH" >> "\$STATUS"

        elif [[ "${params.checkm2_auto_download_db}" == "true" ]]; then
            echo "No CheckM2 database found. Downloading into: \$CHECKM2_DB_DIR" >> "\$STATUS"

            find "\$CHECKM2_DB_DIR" -maxdepth 1 -type f \\( \\
                -name '*.part' -o \\
                -name '*.tmp' -o \\
                -name '*.tar' -o \\
                -name '*.tar.gz' -o \\
                -name '*.tgz' -o \\
                -name '*.zip' \\
            \\) -delete 2>/dev/null || true

            echo "Trying native CheckM2 downloader:" >> "\$STATUS"
            echo "  checkm2 database --download --path \$CHECKM2_DB_DIR" >> "\$STATUS"

            set +e
            checkm2 database \\
                --download \\
                --path "\$CHECKM2_DB_DIR" \\
                >> "\$STATUS" 2>&1
            CHECKM2_NATIVE_DOWNLOAD_STATUS="\$?"
            set -e

            echo "Native CheckM2 downloader exit status: \$CHECKM2_NATIVE_DOWNLOAD_STATUS" >> "\$STATUS"

            CHECKM2_DB_PATH="\$(find "\$CHECKM2_DB_DIR" -type f -name 'uniref100.KO.1.dmnd' 2>/dev/null | head -n 1 || true)"

            if [[ -z "\$CHECKM2_DB_PATH" ]]; then
                echo "Native downloader did not produce uniref100.KO.1.dmnd." >> "\$STATUS"
                echo "Falling back to Zenodo record ${params.checkm2_zenodo_record}" >> "\$STATUS"

                "\$CHECKM2_ENV/bin/python" - \\
                    "\$CHECKM2_DB_DIR" \\
                    "${params.checkm2_zenodo_record}" \\
                    "\$STATUS" <<'PY'
import json
import shutil
import sys
import tarfile
import zipfile
import urllib.request
from pathlib import Path

db_dir = Path(sys.argv[1])
record = sys.argv[2]
status = Path(sys.argv[3])
db_dir.mkdir(parents=True, exist_ok=True)

def log(msg):
    with status.open("a") as handle:
        print(msg, file=handle)

def download(url, dest):
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()
    with urllib.request.urlopen(url) as response, tmp.open("wb") as out:
        shutil.copyfileobj(response, out)
    if not tmp.exists() or tmp.stat().st_size == 0:
        raise RuntimeError(f"Downloaded file missing or empty: {tmp}")
    tmp.rename(dest)

api = f"https://zenodo.org/api/records/{record}"
log(f"Fetching Zenodo API record: {api}")

with urllib.request.urlopen(api) as response:
    data = json.load(response)

files = data.get("files", [])
candidates = []

for f in files:
    key = f.get("key") or f.get("filename") or ""
    links = f.get("links", {})
    url = links.get("content") or links.get("download") or links.get("self") or ""
    size = f.get("size") or f.get("filesize") or 0
    key_lower = key.lower()

    if not key or not url:
        continue

    if key_lower.endswith((".md5", ".sha256", ".txt", ".tsv", ".json")):
        continue

    if (
        key_lower.endswith((".tar.gz", ".tgz", ".tar", ".zip", ".dmnd"))
        or "checkm2" in key_lower
        or "uniref100" in key_lower
        or "database" in key_lower
    ):
        candidates.append((int(size or 0), key, url))

if not candidates:
    log("Available Zenodo files:")
    for f in files:
        log(str(f.get("key") or f.get("filename") or f))
    raise SystemExit("ERROR: Could not find suitable CheckM2 DB candidate in Zenodo record.")

candidates.sort(reverse=True)
size, key, url = candidates[0]
archive = db_dir / key

log(f"Selected candidate: {key}")
log(f"Candidate size: {size}")
log(f"URL: {url}")

download(url, archive)

log(f"Downloaded: {archive}")
log(f"Downloaded size bytes: {archive.stat().st_size}")

if archive.name.endswith((".tar.gz", ".tgz", ".tar")):
    if not tarfile.is_tarfile(archive):
        raise SystemExit(f"ERROR: Downloaded file is not a valid tar archive: {archive}")
    with tarfile.open(archive) as tar:
        tar.extractall(db_dir)

elif archive.name.endswith(".zip"):
    with zipfile.ZipFile(archive) as z:
        z.extractall(db_dir)

matches = list(db_dir.rglob("uniref100.KO.1.dmnd"))

if not matches:
    log("Files found after download/extraction:")
    for p in list(db_dir.rglob("*"))[:200]:
        log(str(p))
    raise SystemExit("ERROR: Could not find uniref100.KO.1.dmnd after download/extraction.")

log(f"Found CheckM2 DB: {matches[0]}")
PY

                CHECKM2_DB_PATH="\$(find "\$CHECKM2_DB_DIR" -type f -name 'uniref100.KO.1.dmnd' 2>/dev/null | head -n 1 || true)"
            fi

        else
            echo "ERROR: --checkm2_db_path not supplied and auto-download disabled." >> "\$STATUS"
            exit 1
        fi
    fi

    if [[ ! -s "\$CHECKM2_DB_PATH" ]]; then
        echo "ERROR: CheckM2 database path does not exist or is empty: \$CHECKM2_DB_PATH" >> "\$STATUS"
        echo "Contents of CHECKM2_DB_DIR:" >> "\$STATUS"
        find "\$CHECKM2_DB_DIR" -maxdepth 5 -type f | head -n 200 >> "\$STATUS" 2>&1 || true
        exit 1
    fi

    echo "CHECKM2_ENV=\$CHECKM2_ENV" >> "\$STATUS"
    echo "CHECKM2_DB_PATH=\$CHECKM2_DB_PATH" >> "\$STATUS"
    echo "CheckM2 setup finished: \$(date)" >> "\$STATUS"
    """
}


process RUN_CHECKM2 {
    tag "checkm2"

    publishDir "${params.outdir}/checkm2", mode: params.publish_tool_outputs_mode, pattern: "checkm2_out/**", saveAs: { filename -> filename.replaceFirst(/^checkm2_out\//, '') }
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "checkm2.log"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "checkm2_status.tsv"

    cpus {
        params.threads != null ? params.threads as int : 8
    }

    input:
    path mags_dir
    path setup_status

    output:
    path "checkm2_status.tsv", emit: status
    path "checkm2.log", emit: log_file
    path "checkm2_out/**", emit: checkm2_out

    script:
    """
    set -euo pipefail

    LOG="checkm2.log"

    CHECKM2_ENV="\$(grep '^CHECKM2_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    CHECKM2_DB_PATH="\$(grep '^CHECKM2_DB_PATH=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"

    export PATH="\$CHECKM2_ENV/bin:\$PATH"

    echo "CheckM2 started: \$(date)" > "\$LOG"
    echo "MAG directory: ${mags_dir}" >> "\$LOG"
    echo "CheckM2 env: \$CHECKM2_ENV" >> "\$LOG"
    echo "CheckM2 DB: \$CHECKM2_DB_PATH" >> "\$LOG"
    echo "Threads: ${task.cpus}" >> "\$LOG"

    mkdir -p checkm2_out

    MAG_COUNT="\$(find -L "${mags_dir}" -maxdepth 1 -type f -name '*.${params.checkm2_extension}' | wc -l | tr -d ' ')"
    echo "MAG files matching extension .${params.checkm2_extension}: \$MAG_COUNT" >> "\$LOG"

    if [[ "\$MAG_COUNT" -eq 0 ]]; then
        echo "ERROR: No MAG files found for CheckM2." >> "\$LOG"
        printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > checkm2_status.tsv
        printf 'checkm2\\tfailed\\t1\\t%s\\tNo MAG files found\\n' "${params.outdir}/checkm2" >> checkm2_status.tsv
        exit 1
    fi

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

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "gtdbtk_setup_status.env"

    output:
    path "gtdbtk_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir ? absPath(params.tool_env_dir) : "${absPath(params.outdir)}/conda_envs"
    def env_dir = "${base_env}/gtdbtk-${params.gtdbtk_version}"
    def gtdbtk_data_path = params.gtdbtk_data_path ? absPath(params.gtdbtk_data_path) : ""
    def db_dir = absPath(params.gtdbtk_db_outdir)
    def download_threads = params.threads != null ? params.threads as int : 8
    def conda_pkgs_dir = params.conda_pkgs_dir ? "${absPath(params.conda_pkgs_dir)}/gtdbtk" : "${absPath(params.outdir)}/conda_pkgs/gtdbtk"

    """
    set -euo pipefail

    STATUS="gtdbtk_setup_status.env"
    GTDBTK_ENV="${env_dir}"
    GTDBTK_DATA_PATH="${gtdbtk_data_path}"
    GTDBTK_DB_DIR="${db_dir}"
    GTDBTK_DOWNLOAD_URL="${params.gtdbtk_download_url}"
    DOWNLOAD_THREADS="${download_threads}"

    CONDA_PKGS_DIRS="${conda_pkgs_dir}"
    export CONDA_PKGS_DIRS
    mkdir -p "\$CONDA_PKGS_DIRS"

    echo "GTDB-Tk setup started: \$(date)" > "\$STATUS"
    echo "GTDBTK_ENV=\$GTDBTK_ENV" >> "\$STATUS"
    echo "GTDBTK_DATA_PATH=\$GTDBTK_DATA_PATH" >> "\$STATUS"
    echo "GTDBTK_DB_DIR=\$GTDBTK_DB_DIR" >> "\$STATUS"
    echo "GTDB-Tk version: ${params.gtdbtk_version}" >> "\$STATUS"
    echo "GTDB-Tk download URL: \$GTDBTK_DOWNLOAD_URL" >> "\$STATUS"
    echo "CONDA_PKGS_DIRS=\$CONDA_PKGS_DIRS" >> "\$STATUS"
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
            --override-channels \\
            -c conda-forge \\
            -c bioconda \\
            "gtdbtk=${params.gtdbtk_version}" \\
            "wget" \\
            "curl" \\
            "tar" \\
            "gzip" \\
            "pigz" \\
            >> "\$STATUS" 2>&1
    fi

    if ! check_env "\$GTDBTK_ENV"; then
        echo "ERROR: GTDB-Tk environment failed checks." >> "\$STATUS"
        exit 1
    fi

    export PATH="\$GTDBTK_ENV/bin:\$PATH"

    if [[ -n "\$GTDBTK_DATA_PATH" ]]; then
        echo "Using user-supplied GTDBTK_DATA_PATH: \$GTDBTK_DATA_PATH" >> "\$STATUS"

        if [[ ! -d "\$GTDBTK_DATA_PATH" ]]; then
            echo "ERROR: User-supplied GTDBTK_DATA_PATH does not exist: \$GTDBTK_DATA_PATH" >> "\$STATUS"
            exit 1
        fi

    else
        mkdir -p "\$GTDBTK_DB_DIR"

        EXISTING_DATA="\$(find "\$GTDBTK_DB_DIR" -type f \\( -name 'metadata.txt' -o -name 'VERSION' \\) 2>/dev/null | head -n 1 || true)"

        if [[ -n "\$EXISTING_DATA" ]]; then
            GTDBTK_DATA_PATH="\$(dirname "\$EXISTING_DATA")"
            echo "Existing GTDB-Tk database detected: \$GTDBTK_DATA_PATH" >> "\$STATUS"

        elif [[ "${params.gtdbtk_auto_download_db}" == "true" ]]; then
            echo "GTDB-Tk database not found. Downloading manually into: \$GTDBTK_DB_DIR" >> "\$STATUS"
            echo "Download URL: \$GTDBTK_DOWNLOAD_URL" >> "\$STATUS"

            ARCHIVE="\$GTDBTK_DB_DIR/gtdbtk_r232_data.tar.gz"
            TMP_ARCHIVE="\${ARCHIVE}.part"

            if [[ -s "\$ARCHIVE" ]]; then
                echo "Existing GTDB-Tk archive found. Validating: \$ARCHIVE" >> "\$STATUS"

                set +e
                gzip -t "\$ARCHIVE" >> "\$STATUS" 2>&1
                ARCHIVE_TEST_STATUS="\$?"
                set -e

                if [[ "\$ARCHIVE_TEST_STATUS" -ne 0 ]]; then
                    echo "Existing archive is corrupt/incomplete. Removing: \$ARCHIVE" >> "\$STATUS"
                    rm -f "\$ARCHIVE"
                else
                    echo "Existing archive passed gzip validation." >> "\$STATUS"
                fi
            fi

            rm -f "\$TMP_ARCHIVE"

            if [[ ! -s "\$ARCHIVE" ]]; then
                if command -v wget >/dev/null 2>&1; then
                    echo "Downloading GTDB-Tk DB with wget." >> "\$STATUS"
                    wget \\
                        --tries=3 \\
                        --timeout=120 \\
                        --waitretry=30 \\
                        -O "\$TMP_ARCHIVE" \\
                        "\$GTDBTK_DOWNLOAD_URL" \\
                        >> "\$STATUS" 2>&1
                elif command -v curl >/dev/null 2>&1; then
                    echo "Downloading GTDB-Tk DB with curl." >> "\$STATUS"
                    curl \\
                        -L \\
                        --retry 3 \\
                        --retry-delay 30 \\
                        --connect-timeout 120 \\
                        -o "\$TMP_ARCHIVE" \\
                        "\$GTDBTK_DOWNLOAD_URL" \\
                        >> "\$STATUS" 2>&1
                else
                    echo "ERROR: Neither wget nor curl is available for GTDB-Tk DB download." >> "\$STATUS"
                    exit 1
                fi

                if [[ ! -s "\$TMP_ARCHIVE" ]]; then
                    echo "ERROR: GTDB-Tk temporary archive missing/empty after download: \$TMP_ARCHIVE" >> "\$STATUS"
                    exit 1
                fi

                echo "Validating downloaded archive: \$TMP_ARCHIVE" >> "\$STATUS"

                set +e
                gzip -t "\$TMP_ARCHIVE" >> "\$STATUS" 2>&1
                TMP_ARCHIVE_TEST_STATUS="\$?"
                set -e

                if [[ "\$TMP_ARCHIVE_TEST_STATUS" -ne 0 ]]; then
                    echo "ERROR: Downloaded GTDB-Tk archive failed gzip validation." >> "\$STATUS"
                    rm -f "\$TMP_ARCHIVE"
                    exit 1
                fi

                mv "\$TMP_ARCHIVE" "\$ARCHIVE"
            fi

            find "\$GTDBTK_DB_DIR" \\
                -mindepth 1 \\
                -maxdepth 1 \\
                ! -name "\$(basename "\$ARCHIVE")" \\
                -exec rm -rf {} +

            echo "Extracting GTDB-Tk database archive." >> "\$STATUS"

            if command -v pigz >/dev/null 2>&1; then
                {
                    pigz -dc -p "\$DOWNLOAD_THREADS" "\$ARCHIVE" \\
                        | tar -x -C "\$GTDBTK_DB_DIR" --strip-components=1
                } >> "\$STATUS" 2>&1
            else
                tar \\
                    -xzf "\$ARCHIVE" \\
                    -C "\$GTDBTK_DB_DIR" \\
                    --strip-components=1 \\
                    >> "\$STATUS" 2>&1
            fi

            rm -f "\$ARCHIVE"

            EXISTING_DATA="\$(find "\$GTDBTK_DB_DIR" -type f \\( -name 'metadata.txt' -o -name 'VERSION' \\) 2>/dev/null | head -n 1 || true)"

            if [[ -n "\$EXISTING_DATA" ]]; then
                GTDBTK_DATA_PATH="\$(dirname "\$EXISTING_DATA")"
                echo "GTDB-Tk database extracted to: \$GTDBTK_DATA_PATH" >> "\$STATUS"
            else
                echo "ERROR: GTDB-Tk database extracted, but metadata.txt or VERSION was not found." >> "\$STATUS"
                find "\$GTDBTK_DB_DIR" -maxdepth 4 -type f | head -n 100 >> "\$STATUS" 2>&1 || true
                exit 1
            fi

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

    publishDir "${params.outdir}/gtdbtk", mode: params.publish_tool_outputs_mode, pattern: "gtdbtk_out/**", saveAs: { filename -> filename.replaceFirst(/^gtdbtk_out\//, '') }
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "gtdbtk.log"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "gtdbtk_status.tsv"

    cpus {
        params.threads != null ? params.threads as int : 32
    }

    input:
    path mags_dir
    path setup_status

    output:
    path "gtdbtk_status.tsv", emit: status
    path "gtdbtk.log", emit: log_file
    path "gtdbtk_out/**", emit: gtdbtk_out

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


    MAG_COUNT="\$(find -L "${mags_dir}" -maxdepth 1 -type f -name '*.${params.gtdbtk_extension}' | wc -l | tr -d ' ')"
    echo "MAG files matching extension .${params.gtdbtk_extension}: \$MAG_COUNT" >> "\$LOG"

    if [[ "\$MAG_COUNT" -eq 0 ]]; then
        echo "ERROR: No MAG files found for GTDB-Tk." >> "\$LOG"
        printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > gtdbtk_status.tsv
        printf 'gtdbtk\\tfailed\\t1\\t%s\\tNo MAG files found\\n' "${params.outdir}/gtdbtk" >> gtdbtk_status.tsv
        exit 1
    fi

    set +e
    gtdbtk classify_wf \\
        --genome_dir "${mags_dir}" \\
        --out_dir gtdbtk_out \\
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


process SETUP_EGGNOG {
    tag "setup_eggnog"

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "eggnog_setup_status.env"

    output:
    path "eggnog_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir ? absPath(params.tool_env_dir) : "${absPath(params.outdir)}/conda_envs"
    def env_dir = params.eggnog_env_dir ? absPath(params.eggnog_env_dir) : "${base_env}/eggnog_mapper"
    def data_dir = params.eggnog_data_path ? absPath(params.eggnog_data_path) : (params.eggnog_data_dir ? absPath(params.eggnog_data_dir) : absPath(params.eggnog_db_outdir))
    def eggnog_mmseqs_db_value = params.eggnog_mmseqs_db ? absPath(params.eggnog_mmseqs_db) : ""
    def conda_pkgs_dir = params.conda_pkgs_dir ? "${absPath(params.conda_pkgs_dir)}/eggnog_mapper" : "${absPath(params.outdir)}/conda_pkgs/eggnog_mapper"
    def eggnog_pkg = "eggnog-mapper=${params.eggnog_mapper_version}"
    def eggnog_fixurl_package = params.eggnog_fixurl_package ?: "eggnog-mapper-fixurl"
    def eggnog_download_args = params.eggnog_download_args != null && params.eggnog_download_args.toString().trim()
        ? params.eggnog_download_args.toString().trim()
        : "-y"

    """
    set -euo pipefail

    STATUS="eggnog_setup_status.env"

    EGGNOG_ENV="${env_dir}"
    EGGNOG_DATA_DIR="${data_dir}"
    USER_EGGNOG_MMSEQS_DB="${eggnog_mmseqs_db_value}"
    EGGNOG_DOWNLOAD_ARGS="${eggnog_download_args}"
    EGGNOG_FIXURL_PACKAGE="${eggnog_fixurl_package}"
    EGGNOG_INSTALL_MARKER="\$EGGNOG_ENV/.samwise_eggnog_install_mode"

    CONDA_PKGS_DIRS="${conda_pkgs_dir}"
    export CONDA_PKGS_DIRS
    mkdir -p "\$CONDA_PKGS_DIRS"

    echo "EggNOG-mapper setup started: \$(date)" > "\$STATUS"
    echo "EGGNOG_ENV=\$EGGNOG_ENV" >> "\$STATUS"
    echo "EGGNOG_DATA_DIR=\$EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "USER_EGGNOG_MMSEQS_DB=\${USER_EGGNOG_MMSEQS_DB:-not supplied}" >> "\$STATUS"
    echo "Requested EggNOG package: ${eggnog_pkg}" >> "\$STATUS"
    echo "EggNOG URL fixer enabled: ${params.eggnog_fixurl}" >> "\$STATUS"
    echo "EggNOG URL fixer package: \$EGGNOG_FIXURL_PACKAGE" >> "\$STATUS"
    echo "EGGNOG_DOWNLOAD_ARGS=\$EGGNOG_DOWNLOAD_ARGS" >> "\$STATUS"
    echo "CONDA_PKGS_DIRS=\$CONDA_PKGS_DIRS" >> "\$STATUS"
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

    check_python_version() {
        local prefix="\$1"

        PY_OK="\$("\$prefix/bin/python" - <<'PY'
import sys
major, minor = sys.version_info[:2]
print("true" if (major == 3 and minor >= 7 and minor <= 12) else "false")
PY
)"
        [[ "\$PY_OK" == "true" ]]
    }

    check_env_core() {
        local prefix="\$1"

        [[ -x "\$prefix/bin/python" ]] || return 1
        [[ -x "\$prefix/bin/pip" ]] || return 1
        [[ -x "\$prefix/bin/emapper.py" ]] || return 1
        [[ -x "\$prefix/bin/download_eggnog_data.py" ]] || return 1
        [[ -x "\$prefix/bin/diamond" ]] || return 1
        [[ -x "\$prefix/bin/prodigal" ]] || return 1

        check_python_version "\$prefix" || return 1

        return 0
    }

    check_env_full() {
        local prefix="\$1"

        check_env_core "\$prefix" || return 1

        # Require marker for the new install mode.
        if [[ ! -s "\$prefix/.samwise_eggnog_install_mode" ]]; then
            echo "EggNOG env has no SAMWISE install marker; treating as stale." >> "\$STATUS"
            return 1
        fi

        if ! grep -q '^bioconda_2.1.13_fixurl\$' "\$prefix/.samwise_eggnog_install_mode"; then
            echo "EggNOG env marker is not bioconda_2.1.13_fixurl; treating as stale." >> "\$STATUS"
            echo "Existing marker:" >> "\$STATUS"
            cat "\$prefix/.samwise_eggnog_install_mode" >> "\$STATUS" 2>&1 || true
            return 1
        fi

        if [[ "${params.eggnog_fixurl}" == "true" ]]; then
            [[ -x "\$prefix/bin/eggnog-mapper-fixurl" ]] || return 1
        fi

        return 0
    }

    if [[ -d "\$EGGNOG_ENV" ]]; then
        if check_env_full "\$EGGNOG_ENV"; then
            echo "Existing Bioconda EggNOG + fixurl environment passed checks." >> "\$STATUS"
        else
            echo "Existing EggNOG environment failed checks or is stale. Removing." >> "\$STATUS"
            rm -rf "\$EGGNOG_ENV"
        fi
    fi

    if [[ ! -d "\$EGGNOG_ENV" ]]; then
        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: EggNOG environment missing and --auto_install false." >> "\$STATUS"
            exit 1
        fi

        INSTALLER="\$(find_installer)"

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda found." >> "\$STATUS"
            exit 1
        fi

        mkdir -p "\$(dirname "\$EGGNOG_ENV")"

        echo "Creating EggNOG environment from Bioconda: \$EGGNOG_ENV" >> "\$STATUS"

        "\$INSTALLER" create -y \\
            -p "\$EGGNOG_ENV" \\
            --override-channels \\
            -c conda-forge \\
            -c bioconda \\
            "python>=3.7,<3.13" \\
            "pip" \\
            "${eggnog_pkg}" \\
            "diamond" \\
            "prodigal" \\
            "hmmer" \\
            "mmseqs2" \\
            "wget" \\
            "curl" \\
            "tar" \\
            "gzip" \\
            >> "\$STATUS" 2>&1
    fi

    if ! check_env_core "\$EGGNOG_ENV"; then
        echo "ERROR: EggNOG environment missing required tools after conda creation." >> "\$STATUS"
        echo "Environment bin preview:" >> "\$STATUS"
        ls -lah "\$EGGNOG_ENV/bin" | grep -E 'emapper|download_eggnog|diamond|prodigal|python|pip|hmmer|mmseqs' >> "\$STATUS" 2>&1 || true
        exit 1
    fi

    export PATH="\$EGGNOG_ENV/bin:\$PATH"

    echo "EggNOG executables before URL fix:" >> "\$STATUS"
    command -v python >> "\$STATUS" 2>&1 || true
    python --version >> "\$STATUS" 2>&1 || true
    command -v pip >> "\$STATUS" 2>&1 || true
    command -v emapper.py >> "\$STATUS" 2>&1 || true
    command -v download_eggnog_data.py >> "\$STATUS" 2>&1 || true
    command -v diamond >> "\$STATUS" 2>&1 || true
    command -v prodigal >> "\$STATUS" 2>&1 || true
    emapper.py --version >> "\$STATUS" 2>&1 || true
    echo "----------------------------------------" >> "\$STATUS"

    # ------------------------------------------------------------------
    # Install and run EggNOG-mapper URL fixer.
    # This patches the broken database download URL used by some versions.
    # ------------------------------------------------------------------
    if [[ "${params.eggnog_fixurl}" == "true" ]]; then
        echo "Installing EggNOG-mapper URL fixer: \$EGGNOG_FIXURL_PACKAGE" >> "\$STATUS"

        "\$EGGNOG_ENV/bin/python" -m pip install "\$EGGNOG_FIXURL_PACKAGE" >> "\$STATUS" 2>&1

        if [[ ! -x "\$EGGNOG_ENV/bin/eggnog-mapper-fixurl" ]]; then
            echo "ERROR: eggnog-mapper-fixurl was not installed into \$EGGNOG_ENV/bin" >> "\$STATUS"
            ls -lah "\$EGGNOG_ENV/bin" | grep -E 'eggnog|fixurl|download' >> "\$STATUS" 2>&1 || true
            exit 1
        fi

        echo "Running eggnog-mapper-fixurl." >> "\$STATUS"
        eggnog-mapper-fixurl >> "\$STATUS" 2>&1

        echo "EggNOG-mapper URL fixer finished." >> "\$STATUS"
    else
        echo "EggNOG-mapper URL fixer disabled by --eggnog_fixurl false" >> "\$STATUS"
    fi

    echo "bioconda_2.1.13_fixurl" > "\$EGGNOG_INSTALL_MARKER"

    if ! check_env_full "\$EGGNOG_ENV"; then
        echo "ERROR: EggNOG environment failed final checks after URL fixer setup." >> "\$STATUS"
        exit 1
    fi

    echo "EggNOG executables after URL fix:" >> "\$STATUS"
    command -v emapper.py >> "\$STATUS" 2>&1 || true
    command -v download_eggnog_data.py >> "\$STATUS" 2>&1 || true
    command -v eggnog-mapper-fixurl >> "\$STATUS" 2>&1 || true
    emapper.py --version >> "\$STATUS" 2>&1 || true
    download_eggnog_data.py --help >> "\$STATUS" 2>&1 || true
    echo "----------------------------------------" >> "\$STATUS"

    mkdir -p "\$EGGNOG_DATA_DIR"

    export EGGNOG_DATA_DIR="\$EGGNOG_DATA_DIR"
    export EGGNOG_DATA_PATH="\$EGGNOG_DATA_DIR"

    echo "Checking EggNOG data directory: \$EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "Existing data files preview:" >> "\$STATUS"
    find "\$EGGNOG_DATA_DIR" -maxdepth 3 -type f | head -n 100 >> "\$STATUS" 2>&1 || true

    EGGNOG_SQLITE_DB="\$(find "\$EGGNOG_DATA_DIR" -type f -name 'eggnog.db' 2>/dev/null | head -n 1 || true)"
    EGGNOG_DIAMOND_DB="\$(find "\$EGGNOG_DATA_DIR" -type f -name '*.dmnd' 2>/dev/null | head -n 1 || true)"

    if [[ -n "\$EGGNOG_SQLITE_DB" && -n "\$EGGNOG_DIAMOND_DB" ]]; then
        echo "Existing EggNOG database appears complete for diamond mode." >> "\$STATUS"
        echo "EggNOG SQLite DB: \$EGGNOG_SQLITE_DB" >> "\$STATUS"
        echo "EggNOG DIAMOND DB: \$EGGNOG_DIAMOND_DB" >> "\$STATUS"

    else
        echo "EggNOG database files are missing or incomplete." >> "\$STATUS"
        echo "Detected eggnog.db: \${EGGNOG_SQLITE_DB:-not found}" >> "\$STATUS"
        echo "Detected *.dmnd: \${EGGNOG_DIAMOND_DB:-not found}" >> "\$STATUS"

        if [[ "${params.eggnog_auto_download_db}" == "true" ]]; then
            echo "Downloading EggNOG database to \$EGGNOG_DATA_DIR" >> "\$STATUS"
            echo "Command: download_eggnog_data.py --data_dir \$EGGNOG_DATA_DIR \$EGGNOG_DOWNLOAD_ARGS" >> "\$STATUS"

            set +e
            download_eggnog_data.py \\
                --data_dir "\$EGGNOG_DATA_DIR" \\
                \$EGGNOG_DOWNLOAD_ARGS \\
                >> "\$STATUS" 2>&1
            DOWNLOAD_STATUS="\$?"
            set -e

            echo "EggNOG database downloader exit status: \$DOWNLOAD_STATUS" >> "\$STATUS"

            if [[ "\$DOWNLOAD_STATUS" -ne 0 ]]; then
                echo "ERROR: EggNOG database downloader failed." >> "\$STATUS"
                echo "Inspect this setup status file for the download_eggnog_data.py error above." >> "\$STATUS"
                exit "\$DOWNLOAD_STATUS"
            fi
        else
            echo "ERROR: EggNOG data directory incomplete and --eggnog_auto_download_db false." >> "\$STATUS"
            echo "Provide --eggnog_data_dir /path/to/existing/eggnog-data or enable auto-download." >> "\$STATUS"
            exit 1
        fi
    fi

    echo "Rechecking EggNOG database files after setup/download." >> "\$STATUS"

    EGGNOG_SQLITE_DB="\$(find "\$EGGNOG_DATA_DIR" -type f -name 'eggnog.db' 2>/dev/null | head -n 1 || true)"
    EGGNOG_DIAMOND_DB="\$(find "\$EGGNOG_DATA_DIR" -type f -name '*.dmnd' 2>/dev/null | head -n 1 || true)"

    if [[ -z "\$EGGNOG_SQLITE_DB" ]]; then
        echo "ERROR: EggNOG database missing eggnog.db after setup." >> "\$STATUS"
        echo "Data directory contents:" >> "\$STATUS"
        find "\$EGGNOG_DATA_DIR" -maxdepth 4 -type f | head -n 200 >> "\$STATUS" 2>&1 || true
        exit 1
    fi

    if [[ "${params.eggnog_method}" == "diamond" && -z "\$EGGNOG_DIAMOND_DB" ]]; then
        echo "ERROR: EggNOG diamond mode requested, but no .dmnd database file was found after setup." >> "\$STATUS"
        echo "Data directory contents:" >> "\$STATUS"
        find "\$EGGNOG_DATA_DIR" -maxdepth 4 -type f | head -n 200 >> "\$STATUS" 2>&1 || true
        exit 1
    fi

    MMSEQS_DB="\$USER_EGGNOG_MMSEQS_DB"

    if [[ -z "\$MMSEQS_DB" ]]; then
        if [[ -s "\$EGGNOG_DATA_DIR/mmseqs/mmseqs.db" ]]; then
            MMSEQS_DB="\$EGGNOG_DATA_DIR/mmseqs/mmseqs.db"
        elif [[ -s "\$EGGNOG_DATA_DIR/mmseqs.db" ]]; then
            MMSEQS_DB="\$EGGNOG_DATA_DIR/mmseqs.db"
        else
            MMSEQS_DB=""
        fi
    fi

    echo "Final EggNOG database files:" >> "\$STATUS"
    echo "EggNOG SQLite DB: \$EGGNOG_SQLITE_DB" >> "\$STATUS"
    echo "EggNOG DIAMOND DB: \$EGGNOG_DIAMOND_DB" >> "\$STATUS"
    echo "EggNOG MMseqs DB: \${MMSEQS_DB:-not found/not used}" >> "\$STATUS"

    echo "EGGNOG_ENV=\$EGGNOG_ENV" >> "\$STATUS"
    echo "EGGNOG_DATA_DIR=\$EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "Note: EGGNOG_DATA_DIR is resolved from --eggnog_data_path, --eggnog_data_dir, or default DB outdir." >> "\$STATUS"
    echo "EGGNOG_MMSEQS_DB=\$MMSEQS_DB" >> "\$STATUS"
    echo "EggNOG-mapper setup finished: \$(date)" >> "\$STATUS"
    """
}

process RUN_EGGNOG {
    tag "eggnog_mapper"

    publishDir "${params.outdir}/eggnog", mode: params.publish_tool_outputs_mode, pattern: "eggnog_out/**", saveAs: { filename -> filename.replaceFirst(/^eggnog_out\//, '') }
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "eggnog.log"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "eggnog_status.tsv"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "eggnog_input_manifest.tsv"

    cpus {
        params.threads != null ? params.threads as int : 30
    }

    input:
    path mags_dir
    path setup_status

    output:
    path "eggnog_status.tsv", emit: status
    path "eggnog.log", emit: log_file
    path "eggnog_input_manifest.tsv", emit: input_manifest
    path "eggnog_out/**", emit: eggnog_out

    script:
    """
    set -euo pipefail

    LOG="eggnog.log"

    EGGNOG_ENV="\$(grep '^EGGNOG_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    EGGNOG_DATA_DIR="\$(grep '^EGGNOG_DATA_DIR=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    EGGNOG_MMSEQS_DB="\$(grep '^EGGNOG_MMSEQS_DB=' "${setup_status}" | tail -n 1 | cut -d= -f2- || true)"

    export PATH="\$EGGNOG_ENV/bin:\$PATH"
    export EGGNOG_DATA_DIR="\$EGGNOG_DATA_DIR"
    export EGGNOG_DATA_PATH="\$EGGNOG_DATA_DIR"

    echo "EggNOG-mapper started: \$(date)" > "\$LOG"
    echo "MAG directory: ${mags_dir}" >> "\$LOG"
    echo "EggNOG env: \$EGGNOG_ENV" >> "\$LOG"
    echo "EggNOG data dir: \$EGGNOG_DATA_DIR" >> "\$LOG"
    echo "Threads: ${task.cpus}" >> "\$LOG"

    mkdir -p eggnog_out eggnog_inputs

    python3 - \\
        "${mags_dir}" \\
        "eggnog_inputs/module6_eggnog_input.fasta" \\
        "eggnog_input_manifest.tsv" \\
        "\$LOG" <<'PY'
import re
import sys
from pathlib import Path

mags_dir = Path(sys.argv[1])
out_fasta = Path(sys.argv[2])
manifest = Path(sys.argv[3])
log_file = Path(sys.argv[4])
suffixes = [".fa", ".fna", ".fasta"]

def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)

def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("_")
    return value or "unnamed"

def strip_suffix(name):
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem

def iter_fasta(path):
    header = None
    seq_parts = []
    with path.open() as handle:
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

fasta_files = []
for path in sorted(mags_dir.iterdir()):
    if path.is_file() and any(path.name.endswith(suffix) for suffix in suffixes):
        fasta_files.append(path)

if not fasta_files:
    raise SystemExit(f"ERROR: No MAG FASTA files found in {mags_dir}")

total_records = 0
total_bp = 0

with out_fasta.open("w") as out, manifest.open("w") as man:
    print("mag_id", "source_fasta", "original_contig_id", "eggnog_contig_id", "bp", sep="\\t", file=man)

    for fasta in fasta_files:
        mag_id = safe_id(strip_suffix(fasta.name))

        for idx, (header, seq) in enumerate(iter_fasta(fasta), start=1):
            if not seq:
                continue

            original_contig = header.split()[0] if header.strip() else f"contig_{idx}"
            new_contig = f"{mag_id}|{safe_id(original_contig)}"

            print(f">{new_contig}", file=out)
            for start in range(0, len(seq), 80):
                print(seq[start:start+80], file=out)

            print(mag_id, str(fasta), original_contig, new_contig, len(seq), sep="\\t", file=man)

            total_records += 1
            total_bp += len(seq)

log(f"EggNOG input FASTA files: {len(fasta_files)}")
log(f"EggNOG input contigs: {total_records}")
log(f"EggNOG input bp: {total_bp}")
PY

    INPUT_FASTA="eggnog_inputs/module6_eggnog_input.fasta"

    if [[ ! -s "\$INPUT_FASTA" ]]; then
        echo "ERROR: EggNOG input FASTA empty." >> "\$LOG"
        printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > eggnog_status.tsv
        printf 'eggnog\\tfailed\\t1\\t%s\\tInput FASTA empty\\n' "${params.outdir}/eggnog" >> eggnog_status.tsv
        exit 1
    fi

    CONTIG_COUNT="\$(grep -c '^>' "\$INPUT_FASTA" || true)"
    echo "EggNOG input contig count: \$CONTIG_COUNT" >> "\$LOG"

    MMSEQS_ARG=""
    if [[ -n "\${EGGNOG_MMSEQS_DB:-}" && -s "\$EGGNOG_MMSEQS_DB" ]]; then
        MMSEQS_ARG="--mmseqs_db \$EGGNOG_MMSEQS_DB"
    fi

    echo "Running EggNOG-mapper:" >> "\$LOG"
    echo "emapper.py -m ${params.eggnog_method} --cpu ${task.cpus} -i \$INPUT_FASTA --itype ${params.eggnog_itype} --genepred ${params.eggnog_genepred} --trans_table ${params.eggnog_trans_table} --data_dir \$EGGNOG_DATA_DIR \$MMSEQS_ARG --output ${params.eggnog_output_prefix} --output_dir eggnog_out --excel ${params.eggnog_extra_args}" >> "\$LOG"

    set +e
    emapper.py \\
        -m ${params.eggnog_method} \\
        --cpu ${task.cpus} \\
        -i "\$INPUT_FASTA" \\
        --itype ${params.eggnog_itype} \\
        --genepred ${params.eggnog_genepred} \\
        --trans_table ${params.eggnog_trans_table} \\
        --data_dir "\$EGGNOG_DATA_DIR" \\
        \$MMSEQS_ARG \\
        --output "${params.eggnog_output_prefix}" \\
        --output_dir eggnog_out \\
        --excel \\
        ${params.eggnog_extra_args} \\
        >> "\$LOG" 2>&1
    STATUS="\$?"
    set -e

    printf 'tool\\tstatus\\texit_status\\toutput_dir\\tinput_contigs\\tmessage\\n' > eggnog_status.tsv

    if [[ "\$STATUS" -ne 0 ]]; then
        if [[ "${params.eggnog_fail_nonfatal}" == "true" ]]; then
            printf 'eggnog\\tfailed_nonfatal\\t%s\\t%s\\t%s\\tEggNOG failed; workflow continued\\n' "\$STATUS" "${params.outdir}/eggnog" "\$CONTIG_COUNT" >> eggnog_status.tsv
            exit 0
        else
            printf 'eggnog\\tfailed\\t%s\\t%s\\t%s\\tEggNOG failed\\n' "\$STATUS" "${params.outdir}/eggnog" "\$CONTIG_COUNT" >> eggnog_status.tsv
            exit "\$STATUS"
        fi
    fi

    OUTPUT_FILES="\$(find eggnog_out -maxdepth 2 -type f | wc -l | tr -d ' ')"
    printf 'eggnog\\tcompleted\\t0\\t%s\\t%s\\tEggNOG completed; output_files=%s\\n' "${params.outdir}/eggnog" "\$CONTIG_COUNT" "\$OUTPUT_FILES" >> eggnog_status.tsv

    echo "EggNOG output files: \$OUTPUT_FILES" >> "\$LOG"
    echo "EggNOG-mapper finished: \$(date)" >> "\$LOG"
    """
}


process WRITE_MODULE6_SUMMARY {
    tag "write_module6_summary"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "module6_run_summary.tsv"

    input:
    path input_stats
    path tool_status_files

    output:
    path "module6_run_summary.tsv", emit: summary

    script:
    def status_files = tool_status_files.collect { status_file -> status_file.name }.join(' ')

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
