#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Module 6: MAG annotate / quality / taxonomy.
 */

params.working_dir = null
params.output_dir = null
params.input_mag_dir = null
params.input_mag_manifest = null
params.mag_extension = "fa"
params.run_drep = true
params.drep_version = null
params.drep_env_dir = null
params.drep_extension = "fa"
params.drep_threads = null
params.drep_extra_args = "-sa 0.99 -comp 50 -con 10"
params.threads = null
params.tool_env_dir = null
params.auto_install = true
params.conda_pkgs_dir = null
params.run_checkm2 = true
params.run_gtdbtk = true
params.run_eggnog = false
params.checkm2_version = null
params.gtdbtk_version = "2.7.2"
params.checkm2_db_path = null
params.checkm2_db_dir = null
params.checkm2_zenodo_record = "14897628"
params.checkm2_auto_download_db = true
params.checkm2_extension = "fa"
params.gtdbtk_data_path = null
params.gtdbtk_db_dir = null
params.gtdbtk_auto_download_db = true
params.gtdbtk_extension = "fa"
params.gtdbtk_download_url = "https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz"
params.gtdbtk_pplacer_cpus = 2
params.gtdbtk_extra_args = ""
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
params.publish_mags_mode = "copy"
params.publish_tool_outputs_mode = "copy"
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.outdir = "${params.results_dir}/module_6_magannotate"
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
    def run_drep = params.run_drep.toString().toBoolean()
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

    if (run_drep && !run_checkm2) {
        error(
            """
            dRep is enabled, but CheckM2 is disabled.

            This dRep implementation requires CheckM2 because it uses:
              checkm2_out/**/quality_report.tsv

            to create:
              modified_quality_report.csv

            Required:
              --run_drep true --run_checkm2 true
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
    log.info("Run dRep: ${run_drep}")
    log.info("Run CheckM2: ${run_checkm2}")
    log.info("Run GTDB-Tk: ${run_gtdbtk}")
    log.info("Run EggNOG-mapper: ${run_eggnog}")
    log.info("dRep extra args: ${params.drep_extra_args}")
    log.info("Conda package cache root: ${params.conda_pkgs_dir ?: params.outdir + '/conda_pkgs'}")
    log.info("Per-tool conda package caches will be used under the cache root.")

    PREPARE_MAG_INPUTS(
        channel.value(selected_mag_dir),
        channel.value(selected_mag_manifest),
    )

    /*
     * Default downstream MAG directory.
     * If dRep runs, this is replaced below with RUN_DREP.out.derep_mags_dir.
     */
    def mags_for_annotation = PREPARE_MAG_INPUTS.out.mags_dir

    /*
     * CheckM2 runs once.
     *
     * If dRep is enabled, CheckM2 is run on the prepared MAGs before dRep.
     * Its quality_report.tsv is converted to modified_quality_report.csv
     * and passed into dRep with --genomeInfo.
     */
    if (run_checkm2) {
        SETUP_CHECKM2()

        RUN_CHECKM2(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_CHECKM2.out.status,
        )
    }

    if (run_drep) {
        PREPARE_DREP_GENOME_INFO(
            RUN_CHECKM2.out.quality_report,
        )

        SETUP_DREP()

        RUN_DREP(
            PREPARE_MAG_INPUTS.out.mags_dir,
            SETUP_DREP.out.status,
            PREPARE_DREP_GENOME_INFO.out.genome_info,
        )

        mags_for_annotation = RUN_DREP.out.derep_mags_dir
    }

        if (run_gtdbtk) {
        SETUP_GTDBTK()
        RUN_GTDBTK(
            mags_for_annotation,
            SETUP_GTDBTK.out.status,
        )
    }

    if (run_drep && run_checkm2 && run_gtdbtk) {
        WRITE_DREP_QUALITY_GTDBTK_SUMMARY(
            RUN_DREP.out.derep_mags_dir,
            PREPARE_DREP_GENOME_INFO.out.genome_info,
            RUN_GTDBTK.out.gtdbtk_out,
        )
    }

    if (run_eggnog) {
        SETUP_EGGNOG()
        RUN_EGGNOG(
            mags_for_annotation,
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

    if (run_drep) {
        status_ch = RUN_DREP.out.status.mix(status_ch)
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

process SETUP_DREP {
    tag "setup_drep"

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "drep_setup_status.env"

    output:
    path "drep_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir ? absPath(params.tool_env_dir) : "${absPath(params.outdir)}/conda_envs"
    def env_dir = params.drep_env_dir ? absPath(params.drep_env_dir) : "${base_env}/drep"
    def drep_pkg = params.drep_version ? "drep=${params.drep_version}" : "drep"
    def conda_pkgs_dir = params.conda_pkgs_dir ? "${absPath(params.conda_pkgs_dir)}/drep" : "${absPath(params.outdir)}/conda_pkgs/drep"

    """
    set -euo pipefail

    STATUS="drep_setup_status.env"
    DREP_ENV="${env_dir}"
    CONDA_PKGS_DIRS="${conda_pkgs_dir}"
    DREP_INSTALL_MARKER="\$DREP_ENV/.samwise_drep_install_mode"

    export CONDA_PKGS_DIRS

    mkdir -p "\$CONDA_PKGS_DIRS"

    echo "dRep setup started: \$(date)" > "\$STATUS"
    echo "DREP_ENV=\$DREP_ENV" >> "\$STATUS"
    echo "Requested dRep package: ${drep_pkg}" >> "\$STATUS"
    echo "Pinned Python: python>=3.8,<3.11" >> "\$STATUS"
    echo "Pinned pandas: pandas<2.2" >> "\$STATUS"
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

    check_drep_python_compatibility() {
        local prefix="\$1"

        "\$prefix/bin/python" - <<'PY'
import sys
import io
import pandas as pd

major, minor = sys.version_info[:2]

# dRep currently behaves more safely outside Python 3.13.
if not (major == 3 and minor >= 8 and minor < 11):
    raise SystemExit(f"Bad Python version for this dRep wrapper: {sys.version}")

# This is the exact compatibility issue seen in the failed run.
# dRep calls pandas.read_csv(..., delim_whitespace=True).
text = "a b\\n1 2\\n"
pd.read_csv(io.StringIO(text), delim_whitespace=True)

print("python_ok")
print(sys.version.replace("\\n", " "))
print("pandas_ok")
print(pd.__version__)
PY
    }

    check_env() {
        local prefix="\$1"

        [[ -x "\$prefix/bin/dRep" ]] || return 1
        [[ -x "\$prefix/bin/python" ]] || return 1
        [[ -x "\$prefix/bin/fastANI" ]] || return 1
        [[ -x "\$prefix/bin/mash" ]] || return 1

        check_drep_python_compatibility "\$prefix" >> "\$STATUS" 2>&1 || return 1

        if [[ ! -s "\$prefix/.samwise_drep_install_mode" ]]; then
            echo "dRep env has no SAMWISE install marker; treating as stale." >> "\$STATUS"
            return 1
        fi

        if ! grep -q '^drep_py38_310_pandas_lt22\$' "\$prefix/.samwise_drep_install_mode"; then
            echo "dRep env marker is not drep_py38_310_pandas_lt22; treating as stale." >> "\$STATUS"
            echo "Existing marker:" >> "\$STATUS"
            cat "\$prefix/.samwise_drep_install_mode" >> "\$STATUS" 2>&1 || true
            return 1
        fi

        return 0
    }

    if [[ -d "\$DREP_ENV" ]]; then
        if check_env "\$DREP_ENV"; then
            echo "Existing dRep environment passed checks." >> "\$STATUS"
        else
            echo "Existing dRep environment failed checks or is stale. Removing." >> "\$STATUS"
            rm -rf "\$DREP_ENV"
        fi
    fi

    if [[ ! -d "\$DREP_ENV" ]]; then
        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: dRep environment missing and --auto_install false." >> "\$STATUS"
            exit 1
        fi

        INSTALLER="\$(find_installer)"

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda found." >> "\$STATUS"
            exit 1
        fi

        mkdir -p "\$(dirname "\$DREP_ENV")"

        echo "Creating dRep environment with pinned Python/pandas." >> "\$STATUS"

        "\$INSTALLER" create -y \\
            -p "\$DREP_ENV" \\
            --override-channels \\
            -c conda-forge \\
            -c bioconda \\
            "python>=3.8,<3.11" \\
            "pandas<2.2" \\
            "${drep_pkg}" \\
            "fastani" \\
            "mash" \\
            >> "\$STATUS" 2>&1

        echo "drep_py38_310_pandas_lt22" > "\$DREP_INSTALL_MARKER"
    fi

    if ! check_env "\$DREP_ENV"; then
        echo "ERROR: dRep environment failed final checks." >> "\$STATUS"
        echo "Environment bin preview:" >> "\$STATUS"
        ls -lah "\$DREP_ENV/bin" >> "\$STATUS" 2>&1 || true

        echo "Python/pandas diagnostic:" >> "\$STATUS"
        "\$DREP_ENV/bin/python" - <<'PY' >> "\$STATUS" 2>&1 || true
import sys
print(sys.version)
try:
    import pandas as pd
    print("pandas", pd.__version__)
except Exception as e:
    print("pandas import failed:", e)
PY

        exit 1
    fi

    export PATH="\$DREP_ENV/bin:\$PATH"

    echo "dRep executable:" >> "\$STATUS"
    command -v dRep >> "\$STATUS" 2>&1 || true
    dRep --version >> "\$STATUS" 2>&1 || true

    echo "fastANI executable:" >> "\$STATUS"
    command -v fastANI >> "\$STATUS" 2>&1 || true
    fastANI --version >> "\$STATUS" 2>&1 || true

    echo "mash executable:" >> "\$STATUS"
    command -v mash >> "\$STATUS" 2>&1 || true
    mash --version >> "\$STATUS" 2>&1 || true

    echo "Python/pandas final diagnostic:" >> "\$STATUS"
    "\$DREP_ENV/bin/python" - <<'PY' >> "\$STATUS" 2>&1
import sys
import pandas as pd
print(sys.version)
print("pandas", pd.__version__)
PY

    echo "DREP_ENV=\$DREP_ENV" >> "\$STATUS"
    echo "dRep setup finished: \$(date)" >> "\$STATUS"
    """
}

process RUN_DREP {
    tag "drep_dereplicate"

    publishDir "${params.outdir}/drep_out", mode: 'copy', pattern: "drep_out/**", saveAs: { filename -> filename.replaceFirst(/^drep_out\//, '') }
    publishDir "${params.outdir}/drep_out/dereplicated_genomes", mode: 'copy', pattern: "dereplicated_genomes/*"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "drep.log"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "drep_status.tsv"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "drep_genomeInfo_used.csv"

    cpus {
        params.drep_threads != null ? params.drep_threads as int : (params.threads != null ? params.threads as int : 16)
    }

    input:
    path mags_dir
    path setup_status
    path genome_info

    output:
    path "drep_status.tsv", emit: status
    path "drep.log", emit: log_file
    path "drep_genomeInfo_used.csv", emit: genome_info_used
    path "dereplicated_genomes", emit: derep_mags_dir
    path "drep_out/**", emit: drep_out

    script:
    """
    set -euo pipefail

    LOG="drep.log"

    DREP_ENV="\$(grep '^DREP_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"

    export PATH="\$DREP_ENV/bin:\$PATH"

    echo "dRep dereplication started: \$(date)" > "\$LOG"
    echo "Input MAG directory: ${mags_dir}" >> "\$LOG"
    echo "dRep env: \$DREP_ENV" >> "\$LOG"
    echo "dRep extension: .${params.drep_extension}" >> "\$LOG"
    echo "Threads: ${task.cpus}" >> "\$LOG"
    echo "dRep genomeInfo input: ${genome_info}" >> "\$LOG"
    echo "Extra dRep args: ${params.drep_extra_args}" >> "\$LOG"
    echo "Published dRep output directory: ${params.outdir}/drep_out" >> "\$LOG"
    echo "Published dereplicated genomes directory: ${params.outdir}/drep_out/dereplicated_genomes" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    cp "${genome_info}" drep_genomeInfo_used.csv

    rm -rf drep_out dereplicated_genomes
    mkdir -p dereplicated_genomes

    mapfile -d '' GENOMES < <(
        find -L "${mags_dir}" \\
            -maxdepth 1 \\
            -type f \\
            -name '*.${params.drep_extension}' \\
            -print0 \\
            | sort -z
    )

    MAG_COUNT="\${#GENOMES[@]}"

    echo "Input MAG files found for dRep: \$MAG_COUNT" >> "\$LOG"

    if [[ "\$MAG_COUNT" -eq 0 ]]; then
        echo "ERROR: No MAG files found for dRep." >> "\$LOG"

        printf 'tool\\tstatus\\texit_status\\tinput_genomes\\tdereplicated_genomes\\toutput_dir\\tmessage\\n' > drep_status.tsv
        printf 'drep\\tfailed\\t1\\t0\\t0\\t%s\\tNo MAG files found\\n' "${params.outdir}/drep_out" >> drep_status.tsv

        exit 1
    fi

    if [[ ! -s drep_genomeInfo_used.csv ]]; then
        echo "ERROR: dRep genomeInfo file is missing or empty: drep_genomeInfo_used.csv" >> "\$LOG"

        printf 'tool\\tstatus\\texit_status\\tinput_genomes\\tdereplicated_genomes\\toutput_dir\\tmessage\\n' > drep_status.tsv
        printf 'drep\\tfailed\\t1\\t%s\\t0\\t%s\\tdRep genomeInfo file missing or empty\\n' "\$MAG_COUNT" "${params.outdir}/drep_out" >> drep_status.tsv

        exit 1
    fi

    echo "dRep genomeInfo preview:" >> "\$LOG"
    head -n 10 drep_genomeInfo_used.csv >> "\$LOG" 2>&1 || true

    if [[ "\$MAG_COUNT" -eq 1 ]]; then
        echo "Only one MAG supplied. Skipping dRep clustering and copying the genome as dereplicated." >> "\$LOG"

        mkdir -p drep_out/dereplicated_genomes

        ONE_GENOME="\${GENOMES[0]}"
        ONE_BASE="\$(basename "\$ONE_GENOME")"

        cp -L "\$ONE_GENOME" "drep_out/dereplicated_genomes/\$ONE_BASE"
    else
        echo "Running dRep dereplicate." >> "\$LOG"
        echo "Command template:" >> "\$LOG"
        echo "dRep dereplicate drep_out --processors ${task.cpus} -g <MAG files> --genomeInfo drep_genomeInfo_used.csv ${params.drep_extra_args}" >> "\$LOG"

        set +e
        dRep dereplicate \\
            drep_out \\
            --processors ${task.cpus} \\
            -g "\${GENOMES[@]}" \\
            --genomeInfo drep_genomeInfo_used.csv \\
            ${params.drep_extra_args} \\
            >> "\$LOG" 2>&1

        DREP_EXIT="\$?"
        set -e

        if [[ "\$DREP_EXIT" -ne 0 ]]; then
            echo "ERROR: dRep failed with exit status \$DREP_EXIT" >> "\$LOG"
            echo "----------------------------------------" >> "\$LOG"
            echo "dRep output directory preview after failure:" >> "\$LOG"
            find drep_out -maxdepth 5 \\( -type f -o -type l -o -type d \\) | head -n 300 >> "\$LOG" 2>&1 || true
            echo "----------------------------------------" >> "\$LOG"
            echo "GenomeInfo file used by dRep:" >> "\$LOG"
            cat drep_genomeInfo_used.csv >> "\$LOG" 2>&1 || true

            printf 'tool\\tstatus\\texit_status\\tinput_genomes\\tdereplicated_genomes\\toutput_dir\\tmessage\\n' > drep_status.tsv
            printf 'drep\\tfailed\\t%s\\t%s\\t0\\t%s\\tdRep dereplicate failed\\n' "\$DREP_EXIT" "\$MAG_COUNT" "${params.outdir}/drep_out" >> drep_status.tsv

            exit "\$DREP_EXIT"
        fi
    fi

    if [[ ! -d drep_out/dereplicated_genomes ]]; then
        echo "ERROR: dRep did not create drep_out/dereplicated_genomes." >> "\$LOG"
        echo "dRep output preview:" >> "\$LOG"
        find drep_out -maxdepth 4 \\( -type f -o -type l \\) | head -n 200 >> "\$LOG" 2>&1 || true

        printf 'tool\\tstatus\\texit_status\\tinput_genomes\\tdereplicated_genomes\\toutput_dir\\tmessage\\n' > drep_status.tsv
        printf 'drep\\tfailed\\t1\\t%s\\t0\\t%s\\tdRep dereplicated_genomes directory missing\\n' "\$MAG_COUNT" "${params.outdir}/drep_out" >> drep_status.tsv

        exit 1
    fi

    echo "Replacing symlinks inside drep_out with real copied files/directories." >> "\$LOG"

    while IFS= read -r -d '' LINK_PATH; do
        TARGET_PATH="\$(readlink -f "\$LINK_PATH" || true)"

        if [[ -z "\$TARGET_PATH" || ! -e "\$TARGET_PATH" ]]; then
            echo "WARNING: broken symlink skipped: \$LINK_PATH" >> "\$LOG"
            continue
        fi

        rm -f "\$LINK_PATH"

        if [[ -d "\$TARGET_PATH" ]]; then
            cp -aL "\$TARGET_PATH" "\$LINK_PATH"
        else
            cp -L "\$TARGET_PATH" "\$LINK_PATH"
        fi
    done < <(find drep_out -type l -print0)

    rm -rf dereplicated_genomes
    mkdir -p dereplicated_genomes

    find -L drep_out/dereplicated_genomes \\
        -maxdepth 1 \\
        -type f \\
        -name '*.${params.drep_extension}' \\
        -print0 \\
        | while IFS= read -r -d '' genome; do
            cp -L "\$genome" "dereplicated_genomes/\$(basename "\$genome")"
        done

    DEREP_COUNT="\$(find dereplicated_genomes -maxdepth 1 -type f -name '*.${params.drep_extension}' | wc -l | tr -d ' ')"

    echo "Dereplicated genome count: \$DEREP_COUNT" >> "\$LOG"

    if [[ "\$DEREP_COUNT" -eq 0 ]]; then
        echo "ERROR: No dereplicated genomes found after dRep." >> "\$LOG"

        printf 'tool\\tstatus\\texit_status\\tinput_genomes\\tdereplicated_genomes\\toutput_dir\\tmessage\\n' > drep_status.tsv
        printf 'drep\\tfailed\\t1\\t%s\\t0\\t%s\\tNo dereplicated genomes produced\\n' "\$MAG_COUNT" "${params.outdir}/drep_out" >> drep_status.tsv

        exit 1
    fi

    echo "Final dereplicated genomes copied for downstream/publishing:" >> "\$LOG"
    find dereplicated_genomes -maxdepth 1 -type f -name '*.${params.drep_extension}' -printf '%f\\n' | sort >> "\$LOG" 2>&1 || true

    REMAINING_SYMLINKS="\$(find drep_out -type l | wc -l | tr -d ' ')"

    echo "Remaining symlinks in drep_out: \$REMAINING_SYMLINKS" >> "\$LOG"

    if [[ "\$REMAINING_SYMLINKS" -ne 0 ]]; then
        echo "WARNING: Some symlinks remain in drep_out:" >> "\$LOG"
        find drep_out -type l >> "\$LOG" 2>&1 || true
    fi

    printf 'tool\\tstatus\\texit_status\\tinput_genomes\\tdereplicated_genomes\\toutput_dir\\tmessage\\n' > drep_status.tsv
    printf 'drep\\tcompleted\\t0\\t%s\\t%s\\t%s\\tdRep completed\\n' "\$MAG_COUNT" "\$DEREP_COUNT" "${params.outdir}/drep_out" >> drep_status.tsv

    echo "dRep dereplication finished: \$(date)" >> "\$LOG"
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
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "checkm2_quality_report.tsv"

    cpus {
        params.threads != null ? params.threads as int : 8
    }

    input:
    path mags_dir
    path setup_status

    output:
    path "checkm2_status.tsv", emit: status
    path "checkm2.log", emit: log_file
    path "checkm2_quality_report.tsv", emit: quality_report
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

    QUALITY_REPORT="\$(find checkm2_out -maxdepth 4 -type f -name 'quality_report.tsv' | head -n 1 || true)"

    if [[ -z "\$QUALITY_REPORT" || ! -s "\$QUALITY_REPORT" ]]; then
        echo "ERROR: CheckM2 completed but quality_report.tsv was not found." >> "\$LOG"
        echo "CheckM2 output preview:" >> "\$LOG"
        find checkm2_out -maxdepth 5 -type f | head -n 200 >> "\$LOG" 2>&1 || true

        printf 'checkm2\\tfailed\\t1\\t%s\\tCheckM2 quality_report.tsv missing\\n' "${params.outdir}/checkm2" >> checkm2_status.tsv

        exit 1
    fi

    cp "\$QUALITY_REPORT" checkm2_quality_report.tsv

    echo "Copied CheckM2 quality report for dRep:" >> "\$LOG"
    echo "  source: \$QUALITY_REPORT" >> "\$LOG"
    echo "  staged: checkm2_quality_report.tsv" >> "\$LOG"

    printf 'checkm2\\tcompleted\\t0\\t%s\\tCheckM2 completed\\n' "${params.outdir}/checkm2" >> checkm2_status.tsv

    echo "CheckM2 finished: \$(date)" >> "\$LOG"
    """
}

process PREPARE_DREP_GENOME_INFO {
    tag "prepare_drep_genome_info"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "modified_quality_report.csv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "prepare_drep_genome_info.log"

    input:
    path checkm2_quality_report

    output:
    path "modified_quality_report.csv", emit: genome_info
    path "prepare_drep_genome_info.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG="prepare_drep_genome_info.log"

    echo "Preparing dRep genomeInfo file from CheckM2 quality report: \$(date)" > "\$LOG"
    echo "Input CheckM2 quality report: ${checkm2_quality_report}" >> "\$LOG"
    echo "Output dRep genomeInfo CSV: modified_quality_report.csv" >> "\$LOG"
    echo "Required dRep headers: genome,completeness,contamination" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    python3 - \\
        "${checkm2_quality_report}" \\
        "modified_quality_report.csv" \\
        "\$LOG" <<'PY'

import csv
import sys
from pathlib import Path

input_tsv = Path(sys.argv[1])
output_csv = Path(sys.argv[2])
log_file = Path(sys.argv[3])

def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)

if not input_tsv.exists() or input_tsv.stat().st_size == 0:
    raise SystemExit(f"ERROR: CheckM2 quality report missing or empty: {input_tsv}")

with input_tsv.open(newline="") as inp:
    reader = csv.DictReader(inp, delimiter="\\t")

    required = ["Name", "Completeness", "Contamination"]

    if reader.fieldnames is None:
        raise SystemExit("ERROR: CheckM2 quality report has no header row.")

    missing = [col for col in required if col not in reader.fieldnames]

    if missing:
        log(f"Input header fields were: {reader.fieldnames}")
        raise SystemExit(f"ERROR: CheckM2 quality report missing required columns: {missing}")

    written = 0

    with output_csv.open("w", newline="") as out:
        writer = csv.writer(out)

        # dRep is case-sensitive and requires exactly these headers.
        writer.writerow(["genome", "completeness", "contamination"])

        for row in reader:
            genome = str(row.get("Name", "")).strip()
            completeness = str(row.get("Completeness", "")).strip()
            contamination = str(row.get("Contamination", "")).strip()

            if not genome:
                continue

            # CheckM2 usually reports the genome name without .fa.
            # dRep wants this to match the basename passed through -g.
            if not genome.endswith(".fa"):
                genome = genome + ".fa"

            writer.writerow([genome, completeness, contamination])
            written += 1

if written == 0:
    raise SystemExit("ERROR: No rows written to modified_quality_report.csv")

log(f"Rows written to modified_quality_report.csv: {written}")
log("modified_quality_report.csv preview:")

with output_csv.open() as handle:
    for idx, line in enumerate(handle):
        if idx >= 10:
            break
        log(line.rstrip("\\n"))

PY

    echo "Finished preparing dRep genomeInfo file: \$(date)" >> "\$LOG"
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
    def pplacer_cpus = params.gtdbtk_pplacer_cpus != null ? params.gtdbtk_pplacer_cpus as int : 6
    def gtdbtk_extra_args = params.gtdbtk_extra_args != null ? params.gtdbtk_extra_args.toString().trim() : ""

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
    echo "GTDB-Tk CPUs: ${task.cpus}" >> "\$LOG"
    echo "GTDB-Tk pplacer CPUs: ${pplacer_cpus}" >> "\$LOG"
    echo "GTDB-Tk extra args: ${gtdbtk_extra_args}" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    MAG_COUNT="\$(find -L "${mags_dir}" -maxdepth 1 -type f -name '*.${params.gtdbtk_extension}' | wc -l | tr -d ' ')"

    echo "MAG files matching extension .${params.gtdbtk_extension}: \$MAG_COUNT" >> "\$LOG"

    if [[ "\$MAG_COUNT" -eq 0 ]]; then
        echo "ERROR: No MAG files found for GTDB-Tk." >> "\$LOG"

        printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > gtdbtk_status.tsv
        printf 'gtdbtk\\tfailed\\t1\\t%s\\tNo MAG files found\\n' "${params.outdir}/gtdbtk" >> gtdbtk_status.tsv

        exit 1
    fi

    rm -rf gtdbtk_out
    mkdir -p gtdbtk_out

    echo "Running GTDB-Tk:" >> "\$LOG"
    echo "gtdbtk classify_wf --genome_dir ${mags_dir} --out_dir gtdbtk_out --extension ${params.gtdbtk_extension} --cpus ${task.cpus} --pplacer_cpus ${pplacer_cpus} ${gtdbtk_extra_args}" >> "\$LOG"

    set +e
    gtdbtk classify_wf \\
        --genome_dir "${mags_dir}" \\
        --out_dir gtdbtk_out \\
        --extension ${params.gtdbtk_extension} \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${pplacer_cpus} \\
        ${gtdbtk_extra_args} \\
        >> "\$LOG" 2>&1
    STATUS="\$?"
    set -e

    printf 'tool\\tstatus\\texit_status\\toutput_dir\\tmessage\\n' > gtdbtk_status.tsv

    if [[ "\$STATUS" -ne 0 ]]; then
        echo "ERROR: GTDB-Tk failed with exit status \$STATUS" >> "\$LOG"
        echo "----------------------------------------" >> "\$LOG"
        echo "GTDB-Tk pplacer logs found after failure:" >> "\$LOG"

        find gtdbtk_out \\
            -type f \\
            \\( -name '*pplacer*.out' -o -name '*pplacer*.log' -o -name '*pplacer*.err' \\) \\
            -print \\
            >> "\$LOG" 2>&1 || true

        echo "----------------------------------------" >> "\$LOG"
        echo "Tail of pplacer output files:" >> "\$LOG"

        find gtdbtk_out \\
            -type f \\
            \\( -name '*pplacer*.out' -o -name '*pplacer*.log' -o -name '*pplacer*.err' \\) \\
            -print0 \\
            | while IFS= read -r -d '' pf; do
                echo "" >> "\$LOG"
                echo "### \$pf" >> "\$LOG"
                tail -n 120 "\$pf" >> "\$LOG" 2>&1 || true
            done

        echo "----------------------------------------" >> "\$LOG"
        echo "GTDB-Tk output file preview:" >> "\$LOG"
        find gtdbtk_out -maxdepth 6 -type f | head -n 300 >> "\$LOG" 2>&1 || true

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

process WRITE_DREP_QUALITY_GTDBTK_SUMMARY {
    tag "write_drep_quality_gtdbtk_summary"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "dereplicated_genomes_quality_taxonomy_summary.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "write_drep_quality_gtdbtk_summary.log"

    input:
    path derep_mags_dir
    path drep_genome_info
    path gtdbtk_outputs

    output:
    path "dereplicated_genomes_quality_taxonomy_summary.tsv", emit: summary
    path "write_drep_quality_gtdbtk_summary.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG="write_drep_quality_gtdbtk_summary.log"

    echo "Writing dRep + CheckM2 quality + GTDB-Tk taxonomy summary: \$(date)" > "\$LOG"
    echo "Dereplicated MAG directory: ${derep_mags_dir}" >> "\$LOG"
    echo "dRep genomeInfo CSV: ${drep_genome_info}" >> "\$LOG"
    echo "GTDB-Tk outputs staged from RUN_GTDBTK." >> "\$LOG"
    echo "Output summary: dereplicated_genomes_quality_taxonomy_summary.tsv" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    python3 - \\
        "${derep_mags_dir}" \\
        "${drep_genome_info}" \\
        "${params.drep_extension}" \\
        "dereplicated_genomes_quality_taxonomy_summary.tsv" \\
        "\$LOG" <<'PY'

import csv
import sys
from pathlib import Path

derep_mags_dir = Path(sys.argv[1])
drep_genome_info = Path(sys.argv[2])
extension = sys.argv[3].lstrip(".")
out_tsv = Path(sys.argv[4])
log_file = Path(sys.argv[5])

def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)

def strip_known_fasta_suffix(name):
    name = str(name).strip()
    for suffix in [".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem

def add_name_keys(mapping, name, value):
    name = str(name).strip()
    if not name:
        return

    stem = strip_known_fasta_suffix(name)

    keys = set()
    keys.add(name)
    keys.add(stem)
    keys.add(f"{stem}.{extension}")

    for key in keys:
        mapping[key] = value

if not derep_mags_dir.exists() or not derep_mags_dir.is_dir():
    raise SystemExit(f"ERROR: Dereplicated MAG directory missing: {derep_mags_dir}")

if not drep_genome_info.exists() or drep_genome_info.stat().st_size == 0:
    raise SystemExit(f"ERROR: dRep genomeInfo CSV missing or empty: {drep_genome_info}")

derep_fastas = sorted([
    p for p in derep_mags_dir.iterdir()
    if p.is_file() and p.name.endswith("." + extension)
])

if not derep_fastas:
    raise SystemExit(f"ERROR: No dereplicated .{extension} genomes found in {derep_mags_dir}")

derep_genomes = [p.name for p in derep_fastas]

log(f"Dereplicated genomes found: {len(derep_genomes)}")
log("First dereplicated genomes:")
for name in derep_genomes[:20]:
    log(f"  {name}")

# --------------------------------------------------------------------
# Load quality information from modified_quality_report.csv /
# drep_genomeInfo_used.csv style file.
#
# Required headers:
#   genome,completeness,contamination
# --------------------------------------------------------------------
quality_by_key = {}

with drep_genome_info.open(newline="") as handle:
    reader = csv.DictReader(handle)

    if reader.fieldnames is None:
        raise SystemExit("ERROR: dRep genomeInfo CSV has no header row.")

    required = ["genome", "completeness", "contamination"]
    missing = [col for col in required if col not in reader.fieldnames]

    if missing:
        raise SystemExit(f"ERROR: dRep genomeInfo CSV missing required columns: {missing}")

    loaded_quality_rows = 0

    for row in reader:
        genome = str(row.get("genome", "")).strip()
        completeness = str(row.get("completeness", "")).strip()
        contamination = str(row.get("contamination", "")).strip()

        if not genome:
            continue

        value = {
            "completeness": completeness,
            "contamination": contamination,
        }

        add_name_keys(quality_by_key, genome, value)
        loaded_quality_rows += 1

log(f"Quality rows loaded: {loaded_quality_rows}")

# --------------------------------------------------------------------
# Find GTDB-Tk taxonomy summary files.
#
# Current GTDB-Tk files of interest:
#   gtdbtk.bac120.summary.tsv
#   gtdbtk.ar53.summary.tsv
#
# Both bacteria and archaea are handled when present.
# --------------------------------------------------------------------
search_root = Path(".")
gtdbtk_summary_files = []

for p in sorted(search_root.rglob("gtdbtk.*.summary.tsv")):
    name = p.name

    if name == "gtdbtk.bac120.summary.tsv" or name == "gtdbtk.ar53.summary.tsv":
        gtdbtk_summary_files.append(p)

if not gtdbtk_summary_files:
    log("Could not find gtdbtk.bac120.summary.tsv or gtdbtk.ar53.summary.tsv in staged GTDB-Tk outputs.")
    log("Staged files preview:")
    for p in list(search_root.rglob("*"))[:300]:
        log(str(p))
    raise SystemExit("ERROR: No GTDB-Tk taxonomy summary files found.")

log("GTDB-Tk taxonomy summary files found:")
for p in gtdbtk_summary_files:
    log(f"  {p}")

taxonomy_by_key = {}
taxonomy_headers = []

for summary_file in gtdbtk_summary_files:
    if summary_file.name == "gtdbtk.bac120.summary.tsv":
        marker_set = "bac120"
    elif summary_file.name == "gtdbtk.ar53.summary.tsv":
        marker_set = "ar53"
    else:
        marker_set = "unknown"

    with summary_file.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\\t")

        if reader.fieldnames is None:
            log(f"WARNING: GTDB-Tk summary file has no header and will be skipped: {summary_file}")
            continue

        if "user_genome" not in reader.fieldnames:
            log(f"WARNING: GTDB-Tk summary file lacks user_genome column and will be skipped: {summary_file}")
            log(f"Headers: {reader.fieldnames}")
            continue

        local_headers = ["gtdbtk_marker_set", "gtdbtk_summary_file"] + [
            h for h in reader.fieldnames if h != "user_genome"
        ]

        for h in local_headers:
            if h not in taxonomy_headers:
                taxonomy_headers.append(h)

        row_count = 0

        for row in reader:
            user_genome = str(row.get("user_genome", "")).strip()

            if not user_genome:
                continue

            tax_value = {
                "gtdbtk_marker_set": marker_set,
                "gtdbtk_summary_file": summary_file.name,
            }

            for h in reader.fieldnames:
                if h == "user_genome":
                    continue
                tax_value[h] = row.get(h, "")

            add_name_keys(taxonomy_by_key, user_genome, tax_value)
            row_count += 1

        log(f"Loaded {row_count} taxonomy rows from {summary_file.name}")

log(f"Taxonomy lookup entries loaded: {len(taxonomy_by_key)}")
log(f"GTDB-Tk taxonomy columns retained: {taxonomy_headers}")

missing_quality = []
missing_taxonomy = []

with out_tsv.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\\t")

    header = ["genome", "completeness", "contamination"] + taxonomy_headers
    writer.writerow(header)

    for genome_file in derep_genomes:
        genome_stem = strip_known_fasta_suffix(genome_file)

        quality = (
            quality_by_key.get(genome_file)
            or quality_by_key.get(genome_stem)
            or {}
        )

        taxonomy = (
            taxonomy_by_key.get(genome_file)
            or taxonomy_by_key.get(genome_stem)
            or {}
        )

        completeness = quality.get("completeness", "")
        contamination = quality.get("contamination", "")

        if not quality:
            missing_quality.append(genome_file)

        if not taxonomy:
            missing_taxonomy.append(genome_file)

        row = [
            genome_file,
            completeness,
            contamination,
        ]

        for h in taxonomy_headers:
            row.append(taxonomy.get(h, ""))

        writer.writerow(row)

if missing_quality:
    log(f"WARNING: Missing quality information for {len(missing_quality)} dereplicated genomes.")
    for name in missing_quality[:50]:
        log(f"  missing_quality: {name}")

if missing_taxonomy:
    log(f"WARNING: Missing GTDB-Tk taxonomy for {len(missing_taxonomy)} dereplicated genomes.")
    for name in missing_taxonomy[:50]:
        log(f"  missing_taxonomy: {name}")

log(f"Final summary rows written: {len(derep_genomes)}")
log("Final summary preview:")

with out_tsv.open() as handle:
    for idx, line in enumerate(handle):
        if idx >= 10:
            break
        log(line.rstrip("\\n"))

PY

    echo "Finished writing dRep + CheckM2 quality + GTDB-Tk taxonomy summary: \$(date)" >> "\$LOG"
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
