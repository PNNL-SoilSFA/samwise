#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.working_dir = null
params.min_scaffold_bp = 1000
params.threads = null
params.auto_install = true
params.tool_env_dir = null
params.conda_pkgs_dir = null
params.publish_filtered_assemblies_mode = "copy"
params.publish_tool_outputs_mode = "copy"
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
params.eggnog_output_prefix = "samwise_assembly_eggnog"
params.eggnog_extra_args = ""
params.eggnog_mmseqs_db = null
params.eggnog_fail_nonfatal = false


def absPath(value) {
    def text = value == null ? "" : value.toString().trim()

    if (!text || text == "null" || text == "NA") {
        return ""
    }

    return java.nio.file.Paths.get(text)
        .toAbsolutePath()
        .normalize()
        .toString()
}


def requireWorkingDir(value) {
    def workingDir = absPath(value)

    if (!workingDir) {
        error(
            """
            Missing required SAMWISE root directory.

            Supply the path to the samwise-main directory:

              --working_dir /path/to/samwise-main
            """.stripIndent()
        )
    }

    return workingDir
}


params.results_dir = requireWorkingDir(params.working_dir)

params.module2_assembly_dir =
    "${params.results_dir}/module_2_readassembly/assemblies"

params.outdir =
    "${params.results_dir}/AuxModule_1_assemblyAnnotate"

params.eggnog_db_outdir = params.eggnog_data_path
    ? absPath(params.eggnog_data_path)
    : (
        params.eggnog_data_dir
            ? absPath(params.eggnog_data_dir)
            : "${params.outdir}/databases/eggnog"
    )

/*
 * Main workflow.
 */
workflow {

    def minimum_scaffold_bp = params.min_scaffold_bp as int

    if (minimum_scaffold_bp < 1) {
        error(
            """
            Invalid scaffold cutoff: ${minimum_scaffold_bp}

            --min_scaffold_bp must be at least 1.
            """.stripIndent()
        )
    }

    def samwise_root = params.results_dir
    def selected_assembly_dir = params.module2_assembly_dir
    def selected_output_dir = params.outdir

    def root_path =
        java.nio.file.Paths.get(samwise_root)

    def assembly_path =
        java.nio.file.Paths.get(selected_assembly_dir)

    if (!java.nio.file.Files.isDirectory(root_path)) {
        error(
            """
            SAMWISE root directory does not exist or is not a directory:

              ${samwise_root}
            """.stripIndent()
        )
    }

    if (!java.nio.file.Files.isDirectory(assembly_path)) {
        error(
            """
            The standardized Module 2 assembly directory was not found.

            SAMWISE root:
              ${samwise_root}

            Expected assembly directory:
              ${selected_assembly_dir}
            """.stripIndent()
        )
    }

    log.info("Auxiliary Module 1: Assembly annotation")
    log.info("SAMWISE root: ${samwise_root}")
    log.info("Input assembly directory: ${selected_assembly_dir}")
    log.info("Output directory: ${selected_output_dir}")
    log.info("Minimum scaffold size: ${minimum_scaffold_bp} bp")
    log.info("EggNOG method: ${params.eggnog_method}")
    log.info("EggNOG input type: ${params.eggnog_itype}")
    log.info("Threads: ${params.threads ?: 16}")

    /*
     * Stage the standardized assembly directory for the task.
     *
     * selected_assembly_dir is also passed separately as a value so that
     * manifests contain the original absolute input paths rather than paths
     * inside the Nextflow work directory.
     */
    def assembly_dir_ch = channel.fromPath(
        selected_assembly_dir,
        type: "dir",
        checkIfExists: true,
    )

    PREPARE_ASSEMBLIES(
        assembly_dir_ch,
        channel.value(selected_assembly_dir),
        channel.value(minimum_scaffold_bp),
    )

    SETUP_EGGNOG()

    RUN_EGGNOG(
        PREPARE_ASSEMBLIES.out.combined_fasta,
        PREPARE_ASSEMBLIES.out.scaffold_manifest,
        SETUP_EGGNOG.out.status,
    )

    WRITE_ANNOTATION_SUMMARY(
        PREPARE_ASSEMBLIES.out.filter_stats,
        PREPARE_ASSEMBLIES.out.assembly_manifest,
        PREPARE_ASSEMBLIES.out.scaffold_manifest,
        RUN_EGGNOG.out.status,
    )
}


process PREPARE_ASSEMBLIES {

    tag "prepare_assembly_scaffolds"

    publishDir(
        "${params.outdir}/filtered_assemblies",
        mode: params.publish_filtered_assemblies_mode,
        pattern: "filtered_assemblies/*.fa",
        saveAs: { filename ->
            filename.tokenize('/').last()
        },
    )

    publishDir(
        "${params.outdir}/inputs",
        mode: "copy",
        pattern: "eggnog_assembly_scaffolds.fasta",
    )

    publishDir(
        "${params.outdir}/summary",
        mode: "copy",
        pattern: "assembly_filtering_manifest.tsv",
    )

    publishDir(
        "${params.outdir}/summary",
        mode: "copy",
        pattern: "scaffold_filtering_manifest.tsv",
    )

    publishDir(
        "${params.outdir}/summary",
        mode: "copy",
        pattern: "scaffold_filtering_stats.tsv",
    )

    publishDir(
        "${params.outdir}/logs",
        mode: "copy",
        pattern: "prepare_assemblies.log",
    )

    input:
    path staged_assembly_dir
    val original_assembly_dir
    val minimum_scaffold_bp

    output:
    path "filtered_assemblies/*.fa",
        emit: filtered_assemblies

    path "eggnog_assembly_scaffolds.fasta",
        emit: combined_fasta

    path "assembly_filtering_manifest.tsv",
        emit: assembly_manifest

    path "scaffold_filtering_manifest.tsv",
        emit: scaffold_manifest

    path "scaffold_filtering_stats.tsv",
        emit: filter_stats

    path "prepare_assemblies.log",
        emit: log_file

    script:
    """
    set -euo pipefail

    LOG="prepare_assemblies.log"

    echo "Assembly scaffold preparation started: \$(date)" > "\$LOG"
    echo "Original input directory: ${original_assembly_dir}" >> "\$LOG"
    echo "Staged input directory: ${staged_assembly_dir}" >> "\$LOG"
    echo "Output directory: ${params.outdir}" >> "\$LOG"
    echo "Minimum scaffold length: ${minimum_scaffold_bp} bp" >> "\$LOG"
    echo "Task directory: \$(pwd -P)" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    mkdir -p filtered_assemblies

    python3 - \\
        "${staged_assembly_dir}" \\
        "${original_assembly_dir}" \\
        "${minimum_scaffold_bp}" \\
        "filtered_assemblies" \\
        "eggnog_assembly_scaffolds.fasta" \\
        "assembly_filtering_manifest.tsv" \\
        "scaffold_filtering_manifest.tsv" \\
        "scaffold_filtering_stats.tsv" \\
        "\$LOG" \\
        "${params.outdir}/filtered_assemblies" \\
        "${params.outdir}/inputs/eggnog_assembly_scaffolds.fasta" <<'PY'
import gzip
import re
import sys
from contextlib import ExitStack
from pathlib import Path

(
    staged_input_dir,
    original_input_dir,
    minimum_scaffold_bp,
    filtered_dir,
    combined_fasta,
    assembly_manifest,
    scaffold_manifest,
    stats_file,
    log_file,
    published_filtered_dir,
    published_combined_fasta,
) = sys.argv[1:]

staged_input_dir = Path(staged_input_dir)
original_input_dir = Path(original_input_dir)
minimum_scaffold_bp = int(minimum_scaffold_bp)
filtered_dir = Path(filtered_dir)
combined_fasta = Path(combined_fasta)
assembly_manifest = Path(assembly_manifest)
scaffold_manifest = Path(scaffold_manifest)
stats_file = Path(stats_file)
log_file = Path(log_file)
published_filtered_dir = Path(published_filtered_dir)
published_combined_fasta = Path(published_combined_fasta)

filtered_dir.mkdir(parents=True, exist_ok=True)

FASTA_SUFFIXES = (
    ".renamed.fasta.gz",
    ".renamed.fna.gz",
    ".renamed.fa.gz",
    ".renamed.fasta",
    ".renamed.fna",
    ".renamed.fa",
    ".fasta.gz",
    ".fna.gz",
    ".fa.gz",
    ".fasta",
    ".fna",
    ".fa",
)


def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)


def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("._-")
    return value or "unnamed"


def strip_fasta_suffix(name):
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            return name[:-len(suffix)]

    return Path(name).stem


def open_fasta(path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")

    return path.open("rt", errors="replace")


def read_fasta(path):
    header = None
    sequence_parts = []

    with open_fasta(path) as handle:
        for line in handle:
            # rstrip() avoids Groovy/Nextflow escape processing of newline
            # characters in this embedded Python script.
            line = line.rstrip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise RuntimeError(
                        f"Sequence data found before a FASTA header: {path}"
                    )

                sequence_parts.append(line.strip())

    if header is not None:
        yield header, "".join(sequence_parts)


def write_wrapped(handle, sequence, width=80):
    for start in range(0, len(sequence), width):
        print(sequence[start:start + width], file=handle)


if minimum_scaffold_bp < 1:
    raise SystemExit(
        "ERROR: minimum_scaffold_bp must be at least 1."
    )

if not staged_input_dir.exists():
    raise SystemExit(
        f"ERROR: Staged assembly directory does not exist: "
        f"{staged_input_dir}"
    )

if not staged_input_dir.is_dir():
    raise SystemExit(
        f"ERROR: Staged assembly input is not a directory: "
        f"{staged_input_dir}"
    )

if not original_input_dir.is_absolute():
    raise SystemExit(
        f"ERROR: Original assembly directory is not absolute: "
        f"{original_input_dir}"
    )

assemblies = sorted(
    path
    for path in staged_input_dir.iterdir()
    if path.is_file()
    and path.name.endswith(FASTA_SUFFIXES)
)

if not assemblies:
    raise SystemExit(
        f"ERROR: No assembly FASTA files found in "
        f"{original_input_dir}"
    )

log(f"Assembly FASTA files discovered: {len(assemblies)}")

input_assemblies = 0
retained_assemblies = 0
empty_assemblies = 0

input_scaffolds = 0
retained_scaffolds = 0
removed_scaffolds = 0
empty_scaffolds = 0

input_bp = 0
retained_bp = 0
removed_bp = 0

used_assembly_names = set()
used_combined_headers = set()

with ExitStack() as stack:
    assembly_out = stack.enter_context(
        assembly_manifest.open("w")
    )
    scaffold_out = stack.enter_context(
        scaffold_manifest.open("w")
    )
    combined_out = stack.enter_context(
        combined_fasta.open("w")
    )

    print(
        "assembly_id",
        "source_fasta",
        "filtered_fasta",
        "combined_eggnog_fasta",
        "minimum_scaffold_bp",
        "input_scaffolds",
        "retained_scaffolds",
        "removed_scaffolds",
        "input_bp",
        "retained_bp",
        "removed_bp",
        "status",
        sep="\\t",
        file=assembly_out,
    )

    print(
        "assembly_id",
        "source_fasta",
        "original_scaffold_id",
        "eggnog_scaffold_id",
        "scaffold_bp",
        "minimum_scaffold_bp",
        "retained",
        "filter_reason",
        sep="\\t",
        file=scaffold_out,
    )

    for source in assemblies:
        input_assemblies += 1

        original_source = original_input_dir / source.name

        base_assembly_id = safe_id(
            strip_fasta_suffix(source.name)
        )

        assembly_id = base_assembly_id
        duplicate_number = 1

        while assembly_id in used_assembly_names:
            duplicate_number += 1
            assembly_id = (
                f"{base_assembly_id}_{duplicate_number}"
            )

        used_assembly_names.add(assembly_id)

        filtered_name = f"{assembly_id}.fa"
        filtered_path = filtered_dir / filtered_name

        assembly_input_scaffolds = 0
        assembly_retained_scaffolds = 0
        assembly_removed_scaffolds = 0

        assembly_input_bp = 0
        assembly_retained_bp = 0
        assembly_removed_bp = 0

        used_filtered_headers = set()

        with filtered_path.open("w") as filtered_out:
            for scaffold_index, (header, sequence) in enumerate(
                read_fasta(source),
                start=1,
            ):
                scaffold_length = len(sequence)

                assembly_input_scaffolds += 1
                assembly_input_bp += scaffold_length

                original_scaffold_id = (
                    header.split()[0]
                    if header.strip()
                    else f"scaffold_{scaffold_index}"
                )

                filtered_scaffold_id = safe_id(
                    original_scaffold_id
                )

                if filtered_scaffold_id in used_filtered_headers:
                    filtered_scaffold_id = (
                        f"{filtered_scaffold_id}_"
                        f"duplicate_{scaffold_index}"
                    )

                used_filtered_headers.add(
                    filtered_scaffold_id
                )

                eggnog_scaffold_id = (
                    f"{assembly_id}|{filtered_scaffold_id}"
                )

                if eggnog_scaffold_id in used_combined_headers:
                    eggnog_scaffold_id = (
                        f"{eggnog_scaffold_id}|"
                        f"duplicate_{scaffold_index}"
                    )

                used_combined_headers.add(
                    eggnog_scaffold_id
                )

                if scaffold_length == 0:
                    retain = False
                    filter_reason = "empty_sequence"
                    empty_scaffolds += 1
                elif scaffold_length < minimum_scaffold_bp:
                    retain = False
                    filter_reason = "below_minimum_length"
                else:
                    retain = True
                    filter_reason = "retained"

                print(
                    assembly_id,
                    str(original_source),
                    original_scaffold_id,
                    eggnog_scaffold_id,
                    scaffold_length,
                    minimum_scaffold_bp,
                    str(retain).lower(),
                    filter_reason,
                    sep="\\t",
                    file=scaffold_out,
                )

                if retain:
                    print(
                        f">{filtered_scaffold_id}",
                        file=filtered_out,
                    )
                    write_wrapped(
                        filtered_out,
                        sequence,
                    )

                    print(
                        f">{eggnog_scaffold_id}",
                        file=combined_out,
                    )
                    write_wrapped(
                        combined_out,
                        sequence,
                    )

                    assembly_retained_scaffolds += 1
                    assembly_retained_bp += scaffold_length
                else:
                    assembly_removed_scaffolds += 1
                    assembly_removed_bp += scaffold_length

        if assembly_retained_scaffolds == 0:
            if filtered_path.exists():
                filtered_path.unlink()

            empty_assemblies += 1
            assembly_status = (
                "skipped_no_scaffolds_above_cutoff"
            )
            published_filtered_fasta = ""

            log(
                f"WARNING: No scaffolds >= "
                f"{minimum_scaffold_bp} bp were retained "
                f"for {source.name}"
            )
        else:
            retained_assemblies += 1
            assembly_status = "retained"
            published_filtered_fasta = str(
                published_filtered_dir / filtered_name
            )

        print(
            assembly_id,
            str(original_source),
            published_filtered_fasta,
            str(published_combined_fasta),
            minimum_scaffold_bp,
            assembly_input_scaffolds,
            assembly_retained_scaffolds,
            assembly_removed_scaffolds,
            assembly_input_bp,
            assembly_retained_bp,
            assembly_removed_bp,
            assembly_status,
            sep="\\t",
            file=assembly_out,
        )

        input_scaffolds += assembly_input_scaffolds
        retained_scaffolds += assembly_retained_scaffolds
        removed_scaffolds += assembly_removed_scaffolds

        input_bp += assembly_input_bp
        retained_bp += assembly_retained_bp
        removed_bp += assembly_removed_bp

        log(
            f"{source.name}: "
            f"input_scaffolds={assembly_input_scaffolds}, "
            f"retained_scaffolds={assembly_retained_scaffolds}, "
            f"removed_scaffolds={assembly_removed_scaffolds}, "
            f"input_bp={assembly_input_bp}, "
            f"retained_bp={assembly_retained_bp}, "
            f"removed_bp={assembly_removed_bp}"
        )

if retained_scaffolds == 0:
    if combined_fasta.exists():
        combined_fasta.unlink()

    raise SystemExit(
        "ERROR: No assembly scaffolds passed the minimum "
        f"length cutoff of {minimum_scaffold_bp} bp."
    )

if not combined_fasta.exists() or combined_fasta.stat().st_size == 0:
    raise SystemExit(
        "ERROR: The combined EggNOG input FASTA is missing "
        "or empty."
    )

with stats_file.open("w") as stats:
    print(
        "input_assembly_dir",
        "minimum_scaffold_bp",
        "input_assemblies",
        "retained_assemblies",
        "assemblies_without_retained_scaffolds",
        "input_scaffolds",
        "retained_scaffolds",
        "removed_scaffolds",
        "empty_scaffolds",
        "input_bp",
        "retained_bp",
        "removed_bp",
        "combined_eggnog_fasta",
        sep="\\t",
        file=stats,
    )

    print(
        str(original_input_dir),
        minimum_scaffold_bp,
        input_assemblies,
        retained_assemblies,
        empty_assemblies,
        input_scaffolds,
        retained_scaffolds,
        removed_scaffolds,
        empty_scaffolds,
        input_bp,
        retained_bp,
        removed_bp,
        str(published_combined_fasta),
        sep="\\t",
        file=stats,
    )

log("----------------------------------------")
log(f"Input assembly directory: {original_input_dir}")
log(f"Published filtered directory: {published_filtered_dir}")
log(f"Published combined FASTA: {published_combined_fasta}")
log(f"Input assemblies: {input_assemblies}")
log(f"Retained assemblies: {retained_assemblies}")
log(f"Assemblies omitted: {empty_assemblies}")
log(f"Input scaffolds: {input_scaffolds}")
log(f"Retained scaffolds: {retained_scaffolds}")
log(f"Removed scaffolds: {removed_scaffolds}")
log(f"Empty scaffolds: {empty_scaffolds}")
log(f"Input bp: {input_bp}")
log(f"Retained bp: {retained_bp}")
log(f"Removed bp: {removed_bp}")
PY

    echo "Assembly scaffold preparation finished: \$(date)" >> "\$LOG"
    """
}

/*
 * Install and configure EggNOG-mapper.
 */
process SETUP_EGGNOG {

    tag "setup_eggnog"

    publishDir(
        "${params.outdir}/setup",
        mode: "copy",
        pattern: "eggnog_setup_status.env",
    )

    output:
    path "eggnog_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir
        ? absPath(params.tool_env_dir)
        : "${absPath(params.outdir)}/conda_envs"

    def env_dir = params.eggnog_env_dir
        ? absPath(params.eggnog_env_dir)
        : "${base_env}/eggnog_mapper"

    def data_dir = params.eggnog_data_path
        ? absPath(params.eggnog_data_path)
        : (
            params.eggnog_data_dir
                ? absPath(params.eggnog_data_dir)
                : absPath(params.eggnog_db_outdir)
        )

    def conda_pkgs_dir = params.conda_pkgs_dir
        ? "${absPath(params.conda_pkgs_dir)}/eggnog_mapper"
        : "${absPath(params.outdir)}/conda_pkgs/eggnog_mapper"

    def eggnog_package =
        "eggnog-mapper=${params.eggnog_mapper_version}"

    def fixurl_package =
        params.eggnog_fixurl_package
            ?: "eggnog-mapper-fixurl"

    def download_args =
        params.eggnog_download_args != null &&
        params.eggnog_download_args.toString().trim()
            ? params.eggnog_download_args.toString().trim()
            : "-y"

    def configured_mmseqs_db = params.eggnog_mmseqs_db
        ? absPath(params.eggnog_mmseqs_db)
        : ""

    """
    set -euo pipefail

    STATUS="eggnog_setup_status.env"

    EGGNOG_ENV="${env_dir}"
    EGGNOG_DATA_DIR="${data_dir}"
    USER_EGGNOG_MMSEQS_DB="${configured_mmseqs_db}"
    CONDA_PKGS_DIRS="${conda_pkgs_dir}"

    export CONDA_PKGS_DIRS

    mkdir -p "\$CONDA_PKGS_DIRS"

    echo "EggNOG-mapper setup started: \$(date)" > "\$STATUS"
    echo "EGGNOG_ENV=\$EGGNOG_ENV" >> "\$STATUS"
    echo "EGGNOG_DATA_DIR=\$EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "CONDA_PKGS_DIRS=\$CONDA_PKGS_DIRS" >> "\$STATUS"
    echo "Requested package: ${eggnog_package}" >> "\$STATUS"
    echo "EggNOG method: ${params.eggnog_method}" >> "\$STATUS"
    echo "EggNOG URL fixer enabled: ${params.eggnog_fixurl}" >> "\$STATUS"
    echo "EggNOG download arguments: ${download_args}" >> "\$STATUS"
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

    check_environment() {
        local prefix="\$1"

        [[ -x "\$prefix/bin/python" ]] || return 1
        [[ -x "\$prefix/bin/pip" ]] || return 1
        [[ -x "\$prefix/bin/emapper.py" ]] || return 1
        [[ -x "\$prefix/bin/download_eggnog_data.py" ]] || return 1
        [[ -x "\$prefix/bin/diamond" ]] || return 1
        [[ -x "\$prefix/bin/prodigal" ]] || return 1

        return 0
    }

    if [[ -d "\$EGGNOG_ENV" ]]; then
        if check_environment "\$EGGNOG_ENV"; then
            echo "Existing EggNOG environment passed checks." >> "\$STATUS"
        else
            echo "Existing EggNOG environment failed checks. Removing it." >> "\$STATUS"
            rm -rf "\$EGGNOG_ENV"
        fi
    fi

    if [[ ! -d "\$EGGNOG_ENV" ]]; then
        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: EggNOG environment is missing and --auto_install is false." >> "\$STATUS"
            exit 1
        fi

        INSTALLER="\$(find_installer)"

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda was found in PATH." >> "\$STATUS"
            exit 1
        fi

        mkdir -p "\$(dirname "\$EGGNOG_ENV")"

        echo "Creating EggNOG environment: \$EGGNOG_ENV" >> "\$STATUS"

        "\$INSTALLER" create -y \\
            -p "\$EGGNOG_ENV" \\
            --override-channels \\
            -c conda-forge \\
            -c bioconda \\
            "python>=3.8,<3.13" \\
            "pip" \\
            "${eggnog_package}" \\
            "diamond" \\
            "prodigal" \\
            "hmmer" \\
            "mmseqs2" \\
            "openpyxl" \\
            "wget" \\
            "curl" \\
            "tar" \\
            "gzip" \\
            >> "\$STATUS" 2>&1
    fi

    if ! check_environment "\$EGGNOG_ENV"; then
        echo "ERROR: EggNOG environment failed final checks." >> "\$STATUS"
        echo "Environment bin directory:" >> "\$STATUS"
        ls -lah "\$EGGNOG_ENV/bin" >> "\$STATUS" 2>&1 || true
        exit 1
    fi

    export PATH="\$EGGNOG_ENV/bin:\$PATH"

    echo "EggNOG executable information:" >> "\$STATUS"
    command -v emapper.py >> "\$STATUS" 2>&1
    emapper.py --version >> "\$STATUS" 2>&1 || true
    command -v diamond >> "\$STATUS" 2>&1
    diamond version >> "\$STATUS" 2>&1 || true
    command -v prodigal >> "\$STATUS" 2>&1
    prodigal -v >> "\$STATUS" 2>&1 || true

    if [[ "${params.eggnog_fixurl}" == "true" ]]; then
        echo "Installing EggNOG URL fixer: ${fixurl_package}" >> "\$STATUS"

        "\$EGGNOG_ENV/bin/python" -m pip install \\
            "${fixurl_package}" \\
            >> "\$STATUS" 2>&1

        if [[ ! -x "\$EGGNOG_ENV/bin/eggnog-mapper-fixurl" ]]; then
            echo "ERROR: eggnog-mapper-fixurl was not installed." >> "\$STATUS"
            exit 1
        fi

        echo "Running eggnog-mapper-fixurl." >> "\$STATUS"
        eggnog-mapper-fixurl >> "\$STATUS" 2>&1
    fi

    mkdir -p "\$EGGNOG_DATA_DIR"

    export EGGNOG_DATA_DIR="\$EGGNOG_DATA_DIR"
    export EGGNOG_DATA_PATH="\$EGGNOG_DATA_DIR"

    find_database_file() {
        local pattern="\$1"

        find "\$EGGNOG_DATA_DIR" \\
            -type f \\
            -name "\$pattern" \\
            2>/dev/null \\
            | head -n 1 || true
    }

    EGGNOG_SQLITE_DB="\$(find_database_file 'eggnog.db')"
    EGGNOG_DIAMOND_DB="\$(find_database_file '*.dmnd')"

    echo "Existing EggNOG SQLite database: \${EGGNOG_SQLITE_DB:-not found}" >> "\$STATUS"
    echo "Existing EggNOG DIAMOND database: \${EGGNOG_DIAMOND_DB:-not found}" >> "\$STATUS"

    NEED_DOWNLOAD="false"

    if [[ -z "\$EGGNOG_SQLITE_DB" ]]; then
        NEED_DOWNLOAD="true"
    fi

    if [[ "${params.eggnog_method}" == "diamond" && -z "\$EGGNOG_DIAMOND_DB" ]]; then
        NEED_DOWNLOAD="true"
    fi

    if [[ "\$NEED_DOWNLOAD" == "true" ]]; then
        if [[ "${params.eggnog_auto_download_db}" != "true" ]]; then
            echo "ERROR: EggNOG database is incomplete and automatic download is disabled." >> "\$STATUS"
            exit 1
        fi

        echo "Downloading EggNOG data into: \$EGGNOG_DATA_DIR" >> "\$STATUS"

        set +e

        download_eggnog_data.py \\
            --data_dir "\$EGGNOG_DATA_DIR" \\
            ${download_args} \\
            >> "\$STATUS" 2>&1

        DOWNLOAD_EXIT="\$?"

        set -e

        echo "EggNOG downloader exit status: \$DOWNLOAD_EXIT" >> "\$STATUS"

        if [[ "\$DOWNLOAD_EXIT" -ne 0 ]]; then
            echo "ERROR: EggNOG database download failed." >> "\$STATUS"
            exit "\$DOWNLOAD_EXIT"
        fi
    fi

    EGGNOG_SQLITE_DB="\$(find_database_file 'eggnog.db')"
    EGGNOG_DIAMOND_DB="\$(find_database_file '*.dmnd')"

    if [[ -z "\$EGGNOG_SQLITE_DB" || ! -s "\$EGGNOG_SQLITE_DB" ]]; then
        echo "ERROR: eggnog.db is missing after setup." >> "\$STATUS"
        exit 1
    fi

    if [[ "${params.eggnog_method}" == "diamond" ]]; then
        if [[ -z "\$EGGNOG_DIAMOND_DB" || ! -s "\$EGGNOG_DIAMOND_DB" ]]; then
            echo "ERROR: A DIAMOND database is required but was not found." >> "\$STATUS"
            exit 1
        fi
    fi

    MMSEQS_DB="\$USER_EGGNOG_MMSEQS_DB"

    if [[ -z "\$MMSEQS_DB" ]]; then
        if [[ -s "\$EGGNOG_DATA_DIR/mmseqs/mmseqs.db" ]]; then
            MMSEQS_DB="\$EGGNOG_DATA_DIR/mmseqs/mmseqs.db"
        elif [[ -s "\$EGGNOG_DATA_DIR/mmseqs.db" ]]; then
            MMSEQS_DB="\$EGGNOG_DATA_DIR/mmseqs.db"
        fi
    fi

    echo "Final EggNOG SQLite database: \$EGGNOG_SQLITE_DB" >> "\$STATUS"
    echo "Final EggNOG DIAMOND database: \${EGGNOG_DIAMOND_DB:-not used}" >> "\$STATUS"
    echo "Final EggNOG MMseqs database: \${MMSEQS_DB:-not used}" >> "\$STATUS"

    echo "EGGNOG_ENV=\$EGGNOG_ENV" >> "\$STATUS"
    echo "EGGNOG_DATA_DIR=\$EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "EGGNOG_MMSEQS_DB=\$MMSEQS_DB" >> "\$STATUS"

    echo "EggNOG-mapper setup finished: \$(date)" >> "\$STATUS"
    """
}


/*
 * Run EggNOG-mapper on the combined retained assembly scaffolds.
 */
process RUN_EGGNOG {

    tag "eggnog_assembly_scaffolds"

    publishDir(
        "${params.outdir}/eggnog",
        mode: params.publish_tool_outputs_mode,
        pattern: "eggnog_out/*",
        saveAs: { filename ->
            filename.tokenize('/').last()
        },
    )

    publishDir(
        "${params.outdir}/logs",
        mode: "copy",
        pattern: "eggnog.log",
    )

    publishDir(
        "${params.outdir}/summary",
        mode: "copy",
        pattern: "eggnog_status.tsv",
    )

    publishDir(
        "${params.outdir}/summary",
        mode: "copy",
        pattern: "eggnog_scaffold_manifest.tsv",
    )

    cpus {
        params.threads != null
            ? params.threads as int
            : 16
    }

    input:
    path combined_fasta
    path input_scaffold_manifest
    path setup_status

    output:
    path "eggnog_status.tsv",
        emit: status

    path "eggnog.log",
        emit: log_file

    path "eggnog_scaffold_manifest.tsv",
        emit: scaffold_manifest

    path "eggnog_out/*",
        emit: eggnog_out

    script:
    """
    set -euo pipefail

    LOG="eggnog.log"

    EGGNOG_ENV="\$(grep '^EGGNOG_ENV=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    EGGNOG_DATA_DIR="\$(grep '^EGGNOG_DATA_DIR=' "${setup_status}" | tail -n 1 | cut -d= -f2-)"
    EGGNOG_MMSEQS_DB="\$(grep '^EGGNOG_MMSEQS_DB=' "${setup_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -z "\$EGGNOG_ENV" || ! -d "\$EGGNOG_ENV" ]]; then
        echo "ERROR: Invalid EggNOG environment from setup status: \$EGGNOG_ENV" >&2
        exit 1
    fi

    if [[ -z "\$EGGNOG_DATA_DIR" || ! -d "\$EGGNOG_DATA_DIR" ]]; then
        echo "ERROR: Invalid EggNOG data directory from setup status: \$EGGNOG_DATA_DIR" >&2
        exit 1
    fi

    export PATH="\$EGGNOG_ENV/bin:\$PATH"
    export EGGNOG_DATA_DIR="\$EGGNOG_DATA_DIR"
    export EGGNOG_DATA_PATH="\$EGGNOG_DATA_DIR"

    echo "EggNOG-mapper started: \$(date)" > "\$LOG"
    echo "Input FASTA: ${combined_fasta}" >> "\$LOG"
    echo "Input scaffold manifest: ${input_scaffold_manifest}" >> "\$LOG"
    echo "EggNOG environment: \$EGGNOG_ENV" >> "\$LOG"
    echo "EggNOG data directory: \$EGGNOG_DATA_DIR" >> "\$LOG"
    echo "Method: ${params.eggnog_method}" >> "\$LOG"
    echo "Input type: ${params.eggnog_itype}" >> "\$LOG"
    echo "Gene prediction: ${params.eggnog_genepred}" >> "\$LOG"
    echo "Translation table: ${params.eggnog_trans_table}" >> "\$LOG"
    echo "Threads: ${task.cpus}" >> "\$LOG"
    echo "Output prefix: ${params.eggnog_output_prefix}" >> "\$LOG"
    echo "Extra arguments: ${params.eggnog_extra_args}" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    mkdir -p eggnog_out

    cp "${input_scaffold_manifest}" eggnog_scaffold_manifest.tsv

    printf 'tool\\tstatus\\texit_status\\tminimum_scaffold_bp\\tinput_scaffolds\\tinput_bp\\toutput_dir\\tmessage\\n' > eggnog_status.tsv

    if [[ ! -s "${combined_fasta}" ]]; then
        echo "ERROR: Combined EggNOG input FASTA is missing or empty." >> "\$LOG"

        printf 'eggnog\\tfailed\\t1\\t%s\\t0\\t0\\t%s\\tInput FASTA missing or empty\\n' \\
            "${params.min_scaffold_bp}" \\
            "${params.outdir}/eggnog" \\
            >> eggnog_status.tsv

        exit 1
    fi

    INPUT_SCAFFOLDS="\$(grep -c '^>' "${combined_fasta}" || true)"

    INPUT_BP="\$(python3 - "${combined_fasta}" <<'PY'
import sys
from pathlib import Path

fasta = Path(sys.argv[1])
total = 0

with fasta.open() as handle:
    for line in handle:
        if not line.startswith(">"):
            total += len(line.strip())

print(total)
PY
)"

    echo "Input scaffolds: \$INPUT_SCAFFOLDS" >> "\$LOG"
    echo "Input base pairs: \$INPUT_BP" >> "\$LOG"

    if [[ "\$INPUT_SCAFFOLDS" -eq 0 ]]; then
        echo "ERROR: EggNOG input contains no FASTA records." >> "\$LOG"

        printf 'eggnog\\tfailed\\t1\\t%s\\t0\\t0\\t%s\\tInput FASTA contains no records\\n' \\
            "${params.min_scaffold_bp}" \\
            "${params.outdir}/eggnog" \\
            >> eggnog_status.tsv

        exit 1
    fi

    MMSEQS_ARG=""

    if [[ -n "\${EGGNOG_MMSEQS_DB:-}" && -s "\$EGGNOG_MMSEQS_DB" ]]; then
        MMSEQS_ARG="--mmseqs_db \$EGGNOG_MMSEQS_DB"
    fi

    echo "Running EggNOG-mapper." >> "\$LOG"
    echo "Command:" >> "\$LOG"
    echo "emapper.py -m ${params.eggnog_method} --cpu ${task.cpus} -i ${combined_fasta} --itype ${params.eggnog_itype} --genepred ${params.eggnog_genepred} --trans_table ${params.eggnog_trans_table} --data_dir \$EGGNOG_DATA_DIR \$MMSEQS_ARG --output ${params.eggnog_output_prefix} --output_dir eggnog_out --excel ${params.eggnog_extra_args}" >> "\$LOG"

    set +e

    emapper.py \\
        -m "${params.eggnog_method}" \\
        --cpu "${task.cpus}" \\
        -i "${combined_fasta}" \\
        --itype "${params.eggnog_itype}" \\
        --genepred "${params.eggnog_genepred}" \\
        --trans_table "${params.eggnog_trans_table}" \\
        --data_dir "\$EGGNOG_DATA_DIR" \\
        \$MMSEQS_ARG \\
        --output "${params.eggnog_output_prefix}" \\
        --output_dir eggnog_out \\
        --excel \\
        ${params.eggnog_extra_args} \\
        >> "\$LOG" 2>&1

    EGGNOG_EXIT="\$?"

    set -e

    if [[ "\$EGGNOG_EXIT" -ne 0 ]]; then
        echo "ERROR: EggNOG-mapper failed with exit status \$EGGNOG_EXIT." >> "\$LOG"

        if [[ "${params.eggnog_fail_nonfatal}" == "true" ]]; then
            printf 'eggnog\\tfailed_nonfatal\\t%s\\t%s\\t%s\\t%s\\t%s\\tEggNOG failed; workflow continued\\n' \\
                "\$EGGNOG_EXIT" \\
                "${params.min_scaffold_bp}" \\
                "\$INPUT_SCAFFOLDS" \\
                "\$INPUT_BP" \\
                "${params.outdir}/eggnog" \\
                >> eggnog_status.tsv

            exit 0
        fi

        printf 'eggnog\\tfailed\\t%s\\t%s\\t%s\\t%s\\t%s\\tEggNOG-mapper failed\\n' \\
            "\$EGGNOG_EXIT" \\
            "${params.min_scaffold_bp}" \\
            "\$INPUT_SCAFFOLDS" \\
            "\$INPUT_BP" \\
            "${params.outdir}/eggnog" \\
            >> eggnog_status.tsv

        exit "\$EGGNOG_EXIT"
    fi

    OUTPUT_FILES="\$(find eggnog_out -type f | wc -l | tr -d ' ')"

    if [[ "\$OUTPUT_FILES" -eq 0 ]]; then
        echo "ERROR: EggNOG completed but produced no output files." >> "\$LOG"

        printf 'eggnog\\tfailed\\t1\\t%s\\t%s\\t%s\\t%s\\tEggNOG produced no output files\\n' \\
            "${params.min_scaffold_bp}" \\
            "\$INPUT_SCAFFOLDS" \\
            "\$INPUT_BP" \\
            "${params.outdir}/eggnog" \\
            >> eggnog_status.tsv

        exit 1
    fi

    printf 'eggnog\\tcompleted\\t0\\t%s\\t%s\\t%s\\t%s\\tEggNOG completed; output_files=%s\\n' \\
        "${params.min_scaffold_bp}" \\
        "\$INPUT_SCAFFOLDS" \\
        "\$INPUT_BP" \\
        "${params.outdir}/eggnog" \\
        "\$OUTPUT_FILES" \\
        >> eggnog_status.tsv

    echo "EggNOG output files: \$OUTPUT_FILES" >> "\$LOG"
    echo "EggNOG output listing:" >> "\$LOG"
    find eggnog_out -maxdepth 2 -type f -printf '%p\\n' | sort >> "\$LOG"

    echo "EggNOG-mapper finished: \$(date)" >> "\$LOG"
    """
}


/*
 * Write a combined run summary.
 */
process WRITE_ANNOTATION_SUMMARY {

    tag "write_assembly_annotation_summary"

    publishDir(
        "${params.outdir}/summary",
        mode: "copy",
        pattern: "assembly_annotation_run_summary.tsv",
    )

    input:
    path filter_stats
    path assembly_manifest
    path scaffold_manifest
    path eggnog_status

    output:
    path "assembly_annotation_run_summary.tsv",
        emit: summary

    script:
    """
    set -euo pipefail

    printf 'section\\tsource_file\\n' > assembly_annotation_run_summary.tsv
    printf 'filtering_statistics\\t%s\\n' "${filter_stats}" >> assembly_annotation_run_summary.tsv
    printf 'assembly_manifest\\t%s\\n' "${assembly_manifest}" >> assembly_annotation_run_summary.tsv
    printf 'scaffold_manifest\\t%s\\n' "${scaffold_manifest}" >> assembly_annotation_run_summary.tsv
    printf 'eggnog_status\\t%s\\n' "${eggnog_status}" >> assembly_annotation_run_summary.tsv

    echo "" >> assembly_annotation_run_summary.tsv
    echo "# Scaffold filtering statistics" >> assembly_annotation_run_summary.tsv
    cat "${filter_stats}" >> assembly_annotation_run_summary.tsv

    echo "" >> assembly_annotation_run_summary.tsv
    echo "# EggNOG status" >> assembly_annotation_run_summary.tsv
    cat "${eggnog_status}" >> assembly_annotation_run_summary.tsv
    """
}