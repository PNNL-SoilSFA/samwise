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

    return java.nio.file.Paths
        .get(text)
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

params.module2_assembly_dir = "${params.results_dir}/module_2_readassembly/assemblies"

params.module2b_assembly_dir = "${params.results_dir}/module_2b_coassembly/assemblies"

params.module5_assembly_dir = "${params.results_dir}/module_5_subassembly/assemblies"

params.outdir = "${params.results_dir}/AuxModule_1_assemblyAnnotate"
params.module6_eggnog_data_dir = "${params.results_dir}/module_6_magannotate/databases/eggnog"

params.eggnog_db_outdir = params.eggnog_data_path
    ? absPath(params.eggnog_data_path)
    : (params.eggnog_data_dir
        ? absPath(params.eggnog_data_dir)
        : "${params.outdir}/databases/eggnog")

process CLEAN_AUXMODULE1_PUBLISHED_OUTPUTS {
    tag "clean_auxmodule1_published_outputs"
    cache false

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "auxmodule1_publication_cleanup_status.tsv"

    output:
    path "auxmodule1_publication_cleanup_status.tsv", emit: status

    script:
    """
    set -euo pipefail

    printf 'step\\tstatus\\tmessage\\n' > auxmodule1_publication_cleanup_status.tsv

    for output_dir in \\
        "${params.outdir}/filtered_assemblies" \\
        "${params.outdir}/inputs" \\
        "${params.outdir}/eggnog"; do
        if [[ -e "\$output_dir" ]]; then
            rm -rf "\$output_dir"
        fi
    done

    printf 'auxmodule1_publication_cleanup\\tcompleted\\tRemoved prior managed published-output directories\\n' \\
        >> auxmodule1_publication_cleanup_status.tsv
    """
}

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
    def root_path = java.nio.file.Paths.get(samwise_root)

    if (!java.nio.file.Files.isDirectory(root_path)) {
        error(
            """
            SAMWISE root directory does not exist or is not a directory:

              ${samwise_root}
            """.stripIndent()
        )
    }

    def assembly_sources = [
        [
            source_type: "module2",
            path: params.module2_assembly_dir.toString(),
        ],
        [
            source_type: "coassembly",
            path: params.module2b_assembly_dir.toString(),
        ],
        [
            source_type: "subtractive",
            path: params.module5_assembly_dir.toString(),
        ],
    ]

    def existing_sources = assembly_sources.findAll { source ->
        java.nio.file.Files.isDirectory(
            java.nio.file.Paths.get(source.path)
        )
    }

    if (!existing_sources) {
        error(
            """
            No assembly directories were found.

            Checked:

              ${params.module2_assembly_dir}
              ${params.module2b_assembly_dir}
              ${params.module5_assembly_dir}

            Run at least one assembly-producing module before running
            AuxModule_1_assemblyAnnotate.
            """.stripIndent()
        )
    }

    /*
     * A path input cannot be null. If an optional module was not run,
     * use an existing assembly directory as a staging placeholder and
     * pass a Boolean flag telling PREPARE_ASSEMBLIES not to read it.
     */
    def fallback_assembly_dir = existing_sources[0].path.toString()

    def module2_exists = java.nio.file.Files.isDirectory(
        java.nio.file.Paths.get(
            params.module2_assembly_dir.toString()
        )
    )

    def module2b_exists = java.nio.file.Files.isDirectory(
        java.nio.file.Paths.get(
            params.module2b_assembly_dir.toString()
        )
    )

    def module5_exists = java.nio.file.Files.isDirectory(
        java.nio.file.Paths.get(
            params.module5_assembly_dir.toString()
        )
    )

    def module2_input = module2_exists
        ? params.module2_assembly_dir.toString()
        : fallback_assembly_dir

    def module2b_input = module2b_exists
        ? params.module2b_assembly_dir.toString()
        : fallback_assembly_dir

    def module5_input = module5_exists
        ? params.module5_assembly_dir.toString()
        : fallback_assembly_dir

    log.info("Auxiliary Module 1: Assembly annotation")
    log.info("SAMWISE root: ${samwise_root}")
    log.info("Output directory: ${params.outdir}")
    log.info("Minimum scaffold size: ${minimum_scaffold_bp} bp")
    log.info("Exact full-length scaffold deduplication: enabled")
    log.info("Reverse-complement equivalence: enabled")

    log.info(
        "Module 2 assemblies: " + (module2_exists
            ? params.module2_assembly_dir.toString()
            : "not present")
    )

    log.info(
        "Module 2b coassemblies: " + (module2b_exists
            ? params.module2b_assembly_dir.toString()
            : "not present")
    )

    log.info(
        "Module 5 subtractive assemblies: " + (module5_exists
            ? params.module5_assembly_dir.toString()
            : "not present")
    )

    log.info("EggNOG method: ${params.eggnog_method}")
    log.info("EggNOG input type: ${params.eggnog_itype}")
    log.info("Threads: ${params.threads ?: 16}")

    def module2_dir_ch = channel.fromPath(
        module2_input,
        type: "dir",
        checkIfExists: true,
    )

    def module2b_dir_ch = channel.fromPath(
        module2b_input,
        type: "dir",
        checkIfExists: true,
    )

    def module5_dir_ch = channel.fromPath(
        module5_input,
        type: "dir",
        checkIfExists: true,
    )

    CLEAN_AUXMODULE1_PUBLISHED_OUTPUTS()

    PREPARE_ASSEMBLIES(
        module2_dir_ch,
        module2b_dir_ch,
        module5_dir_ch,
        channel.value(module2_exists),
        channel.value(module2b_exists),
        channel.value(module5_exists),
        channel.value(params.module2_assembly_dir.toString()),
        channel.value(params.module2b_assembly_dir.toString()),
        channel.value(params.module5_assembly_dir.toString()),
        channel.value(minimum_scaffold_bp),
        CLEAN_AUXMODULE1_PUBLISHED_OUTPUTS.out.status,
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


/*
 * Collect assembly FASTAs from all assembly-producing modules, apply the
 * minimum scaffold-length filter, collapse exact full-length duplicates,
 * and create the combined EggNOG input FASTA.
 */
process PREPARE_ASSEMBLIES {

    tag "prepare_filter_and_deduplicate_assembly_scaffolds"

    publishDir "${params.outdir}/filtered_assemblies", mode: params.publish_filtered_assemblies_mode, pattern: "filtered_assemblies/*.fa", saveAs: { filename ->
        filename.tokenize('/').last()
    }

    publishDir "${params.outdir}/inputs", mode: "copy", pattern: "eggnog_assembly_scaffolds.fasta"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "assembly_filtering_manifest.tsv"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "scaffold_filtering_manifest.tsv"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "scaffold_filtering_stats.tsv"

    publishDir "${params.outdir}/logs", mode: "copy", pattern: "prepare_assemblies.log"

    input:
    path staged_module2_dir, stageAs: "module2_assemblies"

    path staged_module2b_dir, stageAs: "module2b_assemblies"

    path staged_module5_dir, stageAs: "module5_assemblies"

    val use_module2
    val use_module2b
    val use_module5

    val original_module2_dir
    val original_module2b_dir
    val original_module5_dir

    val minimum_scaffold_bp
    path publication_cleanup_status

    output:
    path "filtered_assemblies/*.fa", emit: filtered_assemblies

    path "eggnog_assembly_scaffolds.fasta", emit: combined_fasta

    path "assembly_filtering_manifest.tsv", emit: assembly_manifest

    path "scaffold_filtering_manifest.tsv", emit: scaffold_manifest

    path "scaffold_filtering_stats.tsv", emit: filter_stats

    path "prepare_assemblies.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG="prepare_assemblies.log"

    echo "Assembly scaffold preparation started: \$(date)" > "\$LOG"
    echo "Module 2 enabled: ${use_module2}" >> "\$LOG"
    echo "Module 2 directory: ${original_module2_dir}" >> "\$LOG"
    echo "Module 2b enabled: ${use_module2b}" >> "\$LOG"
    echo "Module 2b directory: ${original_module2b_dir}" >> "\$LOG"
    echo "Module 5 enabled: ${use_module5}" >> "\$LOG"
    echo "Module 5 directory: ${original_module5_dir}" >> "\$LOG"
    echo "Minimum scaffold length: ${minimum_scaffold_bp} bp" >> "\$LOG"
    echo "Deduplication identity: 100%" >> "\$LOG"
    echo "Deduplication coverage: 100% of both scaffolds" >> "\$LOG"
    echo "Reverse-complement equivalence: enabled" >> "\$LOG"
    echo "Task directory: \$(pwd -P)" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    mkdir -p filtered_assemblies

    python3 - \\
        "${staged_module2_dir}" \\
        "${staged_module2b_dir}" \\
        "${staged_module5_dir}" \\
        "${use_module2}" \\
        "${use_module2b}" \\
        "${use_module5}" \\
        "${original_module2_dir}" \\
        "${original_module2b_dir}" \\
        "${original_module5_dir}" \\
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
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

(
    staged_module2_dir,
    staged_module2b_dir,
    staged_module5_dir,
    use_module2,
    use_module2b,
    use_module5,
    original_module2_dir,
    original_module2b_dir,
    original_module5_dir,
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

SOURCE_DEFINITIONS = [
    {
        "source_type": "module2",
        "enabled": use_module2.lower() == "true",
        "staged_dir": Path(staged_module2_dir),
        "original_dir": Path(original_module2_dir),
    },
    {
        "source_type": "coassembly",
        "enabled": use_module2b.lower() == "true",
        "staged_dir": Path(staged_module2b_dir),
        "original_dir": Path(original_module2b_dir),
    },
    {
        "source_type": "subtractive",
        "enabled": use_module5.lower() == "true",
        "staged_dir": Path(staged_module5_dir),
        "original_dir": Path(original_module5_dir),
    },
]


def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)


def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._|-]+", "_", value)
    value = value.strip("._-|")
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
            line = line.strip()

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
                        f"Sequence data found before a FASTA header: "
                        f"{path}"
                    )

                sequence_parts.append(line)

    if header is not None:
        yield header, "".join(sequence_parts)


def write_wrapped(handle, sequence, width=80):
    for start in range(0, len(sequence), width):
        print(sequence[start:start + width], file=handle)


COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def canonical_sequence(sequence):
    sequence = sequence.upper()
    reverse = reverse_complement(sequence)

    if reverse < sequence:
        return reverse

    return sequence


def sequence_hash(sequence):
    return hashlib.sha256(
        sequence.encode("ascii", errors="strict")
    ).hexdigest()


if minimum_scaffold_bp < 1:
    raise SystemExit(
        "ERROR: minimum_scaffold_bp must be at least 1."
    )

assemblies = []

for source_definition in SOURCE_DEFINITIONS:
    source_type = source_definition["source_type"]
    enabled = source_definition["enabled"]
    staged_dir = source_definition["staged_dir"]
    original_dir = source_definition["original_dir"]

    if not enabled:
        log(
            f"Assembly source not present; skipping: "
            f"{source_type} ({original_dir})"
        )
        continue

    if not staged_dir.exists():
        raise SystemExit(
            f"ERROR: staged directory does not exist for "
            f"{source_type}: {staged_dir}"
        )

    if not staged_dir.is_dir():
        raise SystemExit(
            f"ERROR: staged input is not a directory for "
            f"{source_type}: {staged_dir}"
        )

    source_fastas = sorted(
        path
        for path in staged_dir.iterdir()
        if path.is_file()
        and path.name.endswith(FASTA_SUFFIXES)
    )

    log(
        f"{source_type} assembly FASTAs discovered: "
        f"{len(source_fastas)}"
    )

    for staged_path in source_fastas:
        assemblies.append(
            {
                "source_type": source_type,
                "staged_path": staged_path,
                "original_path": original_dir / staged_path.name,
            }
        )

if not assemblies:
    raise SystemExit(
        "ERROR: No assembly FASTA files were found in any enabled "
        "assembly directory."
    )

assemblies.sort(
    key=lambda item: (
        item["source_type"],
        item["staged_path"].name,
    )
)

used_assembly_ids = set()
assembly_records = []

for item in assemblies:
    source_type = item["source_type"]
    staged_path = item["staged_path"]

    base_assembly_id = safe_id(
        f"{source_type}_{strip_fasta_suffix(staged_path.name)}"
    )

    assembly_id = base_assembly_id
    duplicate_number = 1

    while assembly_id in used_assembly_ids:
        duplicate_number += 1
        assembly_id = (
            f"{base_assembly_id}_{duplicate_number}"
        )

    used_assembly_ids.add(assembly_id)

    item["assembly_id"] = assembly_id
    assembly_records.append(item)

log(f"Total assembly FASTA files discovered: {len(assembly_records)}")

representatives_by_hash = defaultdict(list)
representatives = []
representatives_by_assembly = defaultdict(list)

scaffold_rows = []
assembly_counts = {}

total_input_scaffolds = 0
total_input_bp = 0

total_passed_length_scaffolds = 0
total_passed_length_bp = 0

total_unique_scaffolds = 0
total_unique_bp = 0

total_duplicate_scaffolds = 0
total_duplicate_bp = 0

total_below_cutoff_scaffolds = 0
total_below_cutoff_bp = 0

total_empty_scaffolds = 0

for assembly in assembly_records:
    assembly_id = assembly["assembly_id"]
    source_type = assembly["source_type"]
    staged_path = assembly["staged_path"]
    original_path = assembly["original_path"]

    counts = {
        "input_scaffolds": 0,
        "input_bp": 0,
        "passed_length_scaffolds": 0,
        "passed_length_bp": 0,
        "unique_scaffolds": 0,
        "unique_bp": 0,
        "duplicate_scaffolds": 0,
        "duplicate_bp": 0,
        "below_cutoff_scaffolds": 0,
        "below_cutoff_bp": 0,
        "empty_scaffolds": 0,
    }

    assembly_counts[assembly_id] = counts

    used_scaffold_ids = set()

    for scaffold_index, (header, sequence) in enumerate(
        read_fasta(staged_path),
        start=1,
    ):
        sequence = "".join(sequence.split()).upper()
        scaffold_bp = len(sequence)

        counts["input_scaffolds"] += 1
        counts["input_bp"] += scaffold_bp

        total_input_scaffolds += 1
        total_input_bp += scaffold_bp

        original_scaffold_id = (
            header.split()[0]
            if header.strip()
            else f"scaffold_{scaffold_index}"
        )

        original_scaffold_id = safe_id(
            original_scaffold_id
        )

        unique_original_scaffold_id = original_scaffold_id
        duplicate_header_number = 1

        while unique_original_scaffold_id in used_scaffold_ids:
            duplicate_header_number += 1
            unique_original_scaffold_id = (
                f"{original_scaffold_id}_"
                f"duplicate_{duplicate_header_number}"
            )

        used_scaffold_ids.add(unique_original_scaffold_id)

        candidate_representative_id = (
            f"{assembly_id}|{unique_original_scaffold_id}"
        )

        representative_id = ""
        sequence_sha256 = ""
        retained_as_representative = False
        duplicate_status = "not_evaluated"
        orientation_to_representative = ""
        filter_reason = ""

        if scaffold_bp == 0:
            filter_reason = "empty_sequence"

            counts["empty_scaffolds"] += 1
            total_empty_scaffolds += 1

        elif scaffold_bp < minimum_scaffold_bp:
            filter_reason = "below_minimum_length"

            counts["below_cutoff_scaffolds"] += 1
            counts["below_cutoff_bp"] += scaffold_bp

            total_below_cutoff_scaffolds += 1
            total_below_cutoff_bp += scaffold_bp

        else:
            counts["passed_length_scaffolds"] += 1
            counts["passed_length_bp"] += scaffold_bp

            total_passed_length_scaffolds += 1
            total_passed_length_bp += scaffold_bp

            canonical = canonical_sequence(sequence)
            sequence_sha256 = sequence_hash(canonical)

            matched_representative = None

            for candidate in representatives_by_hash[
                sequence_sha256
            ]:
                if candidate["canonical_sequence"] == canonical:
                    matched_representative = candidate
                    break

            if matched_representative is None:
                representative_id = (
                    candidate_representative_id
                )

                retained_as_representative = True
                duplicate_status = "representative"
                orientation_to_representative = "same"
                filter_reason = "retained_unique"

                representative = {
                    "representative_id": representative_id,
                    "assembly_id": assembly_id,
                    "source_type": source_type,
                    "source_fasta": str(original_path),
                    "source_scaffold_id": (
                        unique_original_scaffold_id
                    ),
                    "sequence": sequence,
                    "canonical_sequence": canonical,
                    "sequence_sha256": sequence_sha256,
                }

                representatives_by_hash[
                    sequence_sha256
                ].append(representative)

                representatives.append(representative)

                representatives_by_assembly[
                    assembly_id
                ].append(representative)

                counts["unique_scaffolds"] += 1
                counts["unique_bp"] += scaffold_bp

                total_unique_scaffolds += 1
                total_unique_bp += scaffold_bp

            else:
                representative_id = matched_representative[
                    "representative_id"
                ]

                duplicate_status = "exact_duplicate"
                filter_reason = "removed_exact_duplicate"

                representative_sequence = (
                    matched_representative["sequence"]
                )

                if sequence == representative_sequence:
                    orientation_to_representative = "same"
                elif (
                    reverse_complement(sequence)
                    == representative_sequence
                ):
                    orientation_to_representative = (
                        "reverse_complement"
                    )
                else:
                    orientation_to_representative = (
                        "canonical_match"
                    )

                counts["duplicate_scaffolds"] += 1
                counts["duplicate_bp"] += scaffold_bp

                total_duplicate_scaffolds += 1
                total_duplicate_bp += scaffold_bp

        scaffold_rows.append(
            {
                "source_type": source_type,
                "assembly_id": assembly_id,
                "source_fasta": str(original_path),
                "original_scaffold_id": (
                    unique_original_scaffold_id
                ),
                "representative_scaffold_id": (
                    representative_id
                ),
                "sequence_sha256": sequence_sha256,
                "scaffold_bp": scaffold_bp,
                "minimum_scaffold_bp": (
                    minimum_scaffold_bp
                ),
                "retained_as_representative": str(
                    retained_as_representative
                ).lower(),
                "duplicate_status": duplicate_status,
                "orientation_to_representative": (
                    orientation_to_representative
                ),
                "filter_reason": filter_reason,
            }
        )

if not representatives:
    raise SystemExit(
        "ERROR: No assembly scaffolds passed the minimum length "
        f"cutoff of {minimum_scaffold_bp} bp."
    )

with combined_fasta.open("w") as combined_out:
    for representative in representatives:
        print(
            f">{representative['representative_id']}",
            file=combined_out,
        )

        write_wrapped(
            combined_out,
            representative["sequence"],
        )

if not combined_fasta.exists():
    raise SystemExit(
        "ERROR: The combined EggNOG FASTA was not created."
    )

if combined_fasta.stat().st_size == 0:
    raise SystemExit(
        "ERROR: The combined EggNOG FASTA is empty."
    )

with assembly_manifest.open("w") as assembly_out:
    print(
        "source_type",
        "assembly_id",
        "source_fasta",
        "filtered_fasta",
        "combined_eggnog_fasta",
        "minimum_scaffold_bp",
        "input_scaffolds",
        "passed_length_scaffolds",
        "unique_representative_scaffolds",
        "exact_duplicate_scaffolds",
        "below_cutoff_scaffolds",
        "empty_scaffolds",
        "input_bp",
        "passed_length_bp",
        "unique_representative_bp",
        "exact_duplicate_bp",
        "below_cutoff_bp",
        "status",
        sep="\t",
        file=assembly_out,
    )

    for assembly in assembly_records:
        assembly_id = assembly["assembly_id"]
        source_type = assembly["source_type"]
        original_path = assembly["original_path"]

        counts = assembly_counts[assembly_id]
        assembly_representatives = (
            representatives_by_assembly[assembly_id]
        )

        filtered_name = f"{assembly_id}.fa"
        filtered_path = filtered_dir / filtered_name

        if assembly_representatives:
            with filtered_path.open("w") as filtered_out:
                for representative in assembly_representatives:
                    print(
                        f">{representative['representative_id']}",
                        file=filtered_out,
                    )

                    write_wrapped(
                        filtered_out,
                        representative["sequence"],
                    )

            filtered_fasta_value = str(
                published_filtered_dir / filtered_name
            )

            status = "retained_unique_scaffolds"

        elif counts["passed_length_scaffolds"] > 0:
            filtered_fasta_value = ""
            status = (
                "all_passing_scaffolds_were_duplicates"
            )

        else:
            filtered_fasta_value = ""
            status = "no_scaffolds_above_cutoff"

        print(
            source_type,
            assembly_id,
            str(original_path),
            filtered_fasta_value,
            str(published_combined_fasta),
            minimum_scaffold_bp,
            counts["input_scaffolds"],
            counts["passed_length_scaffolds"],
            counts["unique_scaffolds"],
            counts["duplicate_scaffolds"],
            counts["below_cutoff_scaffolds"],
            counts["empty_scaffolds"],
            counts["input_bp"],
            counts["passed_length_bp"],
            counts["unique_bp"],
            counts["duplicate_bp"],
            counts["below_cutoff_bp"],
            status,
            sep="\t",
            file=assembly_out,
        )

with scaffold_manifest.open("w") as scaffold_out:
    print(
        "source_type",
        "assembly_id",
        "source_fasta",
        "original_scaffold_id",
        "representative_scaffold_id",
        "sequence_sha256",
        "scaffold_bp",
        "minimum_scaffold_bp",
        "retained_as_representative",
        "duplicate_status",
        "orientation_to_representative",
        "filter_reason",
        sep="\t",
        file=scaffold_out,
    )

    for row in scaffold_rows:
        print(
            row["source_type"],
            row["assembly_id"],
            row["source_fasta"],
            row["original_scaffold_id"],
            row["representative_scaffold_id"],
            row["sequence_sha256"],
            row["scaffold_bp"],
            row["minimum_scaffold_bp"],
            row["retained_as_representative"],
            row["duplicate_status"],
            row["orientation_to_representative"],
            row["filter_reason"],
            sep="\t",
            file=scaffold_out,
        )

enabled_source_count = sum(
    1
    for source_definition in SOURCE_DEFINITIONS
    if source_definition["enabled"]
)

with stats_file.open("w") as stats:
    print(
        "minimum_scaffold_bp",
        "assembly_sources_enabled",
        "input_assemblies",
        "input_scaffolds",
        "input_bp",
        "passed_length_scaffolds",
        "passed_length_bp",
        "unique_representative_scaffolds",
        "unique_representative_bp",
        "exact_duplicate_scaffolds_removed",
        "exact_duplicate_bp_removed",
        "below_cutoff_scaffolds",
        "below_cutoff_bp",
        "empty_scaffolds",
        "deduplication_identity",
        "deduplication_coverage",
        "reverse_complement_equivalent",
        "combined_eggnog_fasta",
        sep="\t",
        file=stats,
    )

    print(
        minimum_scaffold_bp,
        enabled_source_count,
        len(assembly_records),
        total_input_scaffolds,
        total_input_bp,
        total_passed_length_scaffolds,
        total_passed_length_bp,
        total_unique_scaffolds,
        total_unique_bp,
        total_duplicate_scaffolds,
        total_duplicate_bp,
        total_below_cutoff_scaffolds,
        total_below_cutoff_bp,
        total_empty_scaffolds,
        "100%",
        "100%_reciprocal",
        "true",
        str(published_combined_fasta),
        sep="\t",
        file=stats,
    )

log("----------------------------------------")
log(f"Enabled assembly sources: {enabled_source_count}")
log(f"Input assembly FASTAs: {len(assembly_records)}")
log(f"Input scaffolds: {total_input_scaffolds}")
log(f"Input bp: {total_input_bp}")
log(
    f"Scaffolds passing length filter: "
    f"{total_passed_length_scaffolds}"
)
log(
    f"Scaffolds retained after exact deduplication: "
    f"{total_unique_scaffolds}"
)
log(
    f"Exact duplicate scaffolds removed: "
    f"{total_duplicate_scaffolds}"
)
log(
    f"Scaffolds below cutoff: "
    f"{total_below_cutoff_scaffolds}"
)
log(f"Empty scaffolds: {total_empty_scaffolds}")
log(
    f"Deduplicated EggNOG FASTA: "
    f"{published_combined_fasta}"
)
PY

    echo "Assembly scaffold preparation finished: \$(date)" >> "\$LOG"
    """
}

/*
 * Install and configure EggNOG-mapper.
 */
process SETUP_EGGNOG {

    tag "setup_eggnog"
    cache false

    publishDir "${params.outdir}/setup", mode: "copy", pattern: "eggnog_setup_status.env"

    output:
    path "eggnog_setup_status.env", emit: status

    script:
    def base_env = params.tool_env_dir
        ? absPath(params.tool_env_dir)
        : "${absPath(params.outdir)}/conda_envs"

    def env_dir = params.eggnog_env_dir
        ? absPath(params.eggnog_env_dir)
        : "${base_env}/eggnog_mapper"

    def configured_data_dir = params.eggnog_data_path
        ? absPath(params.eggnog_data_path)
        : (params.eggnog_data_dir
            ? absPath(params.eggnog_data_dir)
            : "")

    def default_data_dir = absPath(params.eggnog_db_outdir)
    def module6_data_dir = absPath(params.module6_eggnog_data_dir)

    def conda_pkgs_dir = params.conda_pkgs_dir
        ? "${absPath(params.conda_pkgs_dir)}/eggnog_mapper"
        : "${absPath(params.outdir)}/conda_pkgs/eggnog_mapper"

    def eggnog_package = "eggnog-mapper=${params.eggnog_mapper_version}"

    def fixurl_package = params.eggnog_fixurl_package ?: "eggnog-mapper-fixurl"

    def download_args = params.eggnog_download_args != null && params.eggnog_download_args.toString().trim()
        ? params.eggnog_download_args.toString().trim()
        : "-y"

    def configured_mmseqs_db = params.eggnog_mmseqs_db
        ? absPath(params.eggnog_mmseqs_db)
        : ""

    """
    set -euo pipefail

    STATUS="eggnog_setup_status.env"

    EGGNOG_ENV="${env_dir}"
    CONFIGURED_EGGNOG_DATA_DIR="${configured_data_dir}"
    DEFAULT_EGGNOG_DATA_DIR="${default_data_dir}"
    MODULE6_EGGNOG_DATA_DIR="${module6_data_dir}"
    USER_EGGNOG_MMSEQS_DB="${configured_mmseqs_db}"
    EGGNOG_INSTALL_MARKER="\$EGGNOG_ENV/.samwise_eggnog_install_mode"
    CONDA_PKGS_DIRS="${conda_pkgs_dir}"

    export CONDA_PKGS_DIRS

    mkdir -p "\$CONDA_PKGS_DIRS"

    echo "EggNOG-mapper setup started: \$(date)" > "\$STATUS"
    echo "EGGNOG_ENV=\$EGGNOG_ENV" >> "\$STATUS"
    echo "Configured EggNOG data directory: \${CONFIGURED_EGGNOG_DATA_DIR:-not supplied}" >> "\$STATUS"
    echo "Module 6 EggNOG data directory candidate: \$MODULE6_EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "Default AuxModule 1 EggNOG data directory: \$DEFAULT_EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "CONDA_PKGS_DIRS=\$CONDA_PKGS_DIRS" >> "\$STATUS"
    echo "Requested package: ${eggnog_package}" >> "\$STATUS"
    echo "EggNOG method: ${params.eggnog_method}" >> "\$STATUS"
    echo "EggNOG URL fixer enabled: ${params.eggnog_fixurl}" >> "\$STATUS"
    echo "EggNOG URL fixer package: ${fixurl_package}" >> "\$STATUS"
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

        if [[ ! -s "\$prefix/.samwise_eggnog_install_mode" ]]; then
            echo "EggNOG environment has no SAMWISE install marker; treating as stale." >> "\$STATUS"
            return 1
        fi

        if ! grep -q '^bioconda_2.1.13_fixurl\$' "\$prefix/.samwise_eggnog_install_mode"; then
            echo "EggNOG environment marker is not bioconda_2.1.13_fixurl; treating as stale." >> "\$STATUS"
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
            echo "Existing EggNOG environment failed checks or is stale. Removing it." >> "\$STATUS"
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

    if ! check_env_core "\$EGGNOG_ENV"; then
        echo "ERROR: EggNOG environment failed core checks." >> "\$STATUS"
        ls -lah "\$EGGNOG_ENV/bin" >> "\$STATUS" 2>&1 || true
        exit 1
    fi

    export PATH="\$EGGNOG_ENV/bin:\$PATH"

    {
        echo "EggNOG executable information:"
        command -v emapper.py || true
        emapper.py --version || true
        command -v diamond || true
        diamond version || true
        command -v prodigal || true
        prodigal -v || true
    } >> "\$STATUS" 2>&1

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

        "\$EGGNOG_ENV/bin/eggnog-mapper-fixurl" \
            >> "\$STATUS" 2>&1
    else
        echo "EggNOG-mapper URL fixer disabled by --eggnog_fixurl false" >> "\$STATUS"
    fi

    echo "bioconda_2.1.13_fixurl" > "\$EGGNOG_INSTALL_MARKER"

    if ! check_env_full "\$EGGNOG_ENV"; then
        echo "ERROR: EggNOG environment failed final checks after URL fixer setup." >> "\$STATUS"
        exit 1
    fi

    has_required_database() {
        local candidate="\$1"
        local sqlite_db
        local diamond_db

        [[ -n "\$candidate" && -d "\$candidate" ]] || return 1

        sqlite_db="\$(find "\$candidate" -type f -name 'eggnog.db' -size +0c -print -quit 2>/dev/null || true)"
        diamond_db="\$(find "\$candidate" -type f -name '*.dmnd' -size +0c -print -quit 2>/dev/null || true)"

        [[ -n "\$sqlite_db" ]] || return 1
        [[ "${params.eggnog_method}" != "diamond" || -n "\$diamond_db" ]]
    }

    if [[ -n "\$CONFIGURED_EGGNOG_DATA_DIR" ]]; then
        EGGNOG_DATA_DIR="\$CONFIGURED_EGGNOG_DATA_DIR"
        echo "Using user-supplied EggNOG data directory: \$EGGNOG_DATA_DIR" >> "\$STATUS"
    elif has_required_database "\$MODULE6_EGGNOG_DATA_DIR"; then
        EGGNOG_DATA_DIR="\$MODULE6_EGGNOG_DATA_DIR"
        echo "Reusing valid Module 6 EggNOG database: \$EGGNOG_DATA_DIR" >> "\$STATUS"
    else
        EGGNOG_DATA_DIR="\$DEFAULT_EGGNOG_DATA_DIR"
        echo "No valid Module 6 EggNOG database found; using AuxModule 1 database directory: \$EGGNOG_DATA_DIR" >> "\$STATUS"
    fi

    mkdir -p "\$EGGNOG_DATA_DIR"

    export EGGNOG_DATA_DIR
    export EGGNOG_DATA_PATH="\$EGGNOG_DATA_DIR"

    find_database_file() {
        local pattern="\$1"

        find "\$EGGNOG_DATA_DIR" \\
            -type f \\
            -name "\$pattern" \\
            -print \\
            2>/dev/null \\
            | sort \\
            | head -n 1 || true
    }

    EGGNOG_DB_PATH="\$(find_database_file 'eggnog.db')"
    EGGNOG_DIAMOND_DB="\$(find_database_file '*.dmnd')"

    echo "Existing EggNOG SQLite database: \${EGGNOG_DB_PATH:-not found}" >> "\$STATUS"
    echo "Existing EggNOG DIAMOND database: \${EGGNOG_DIAMOND_DB:-not found}" >> "\$STATUS"

    NEED_DOWNLOAD="false"

    if [[ -z "\$EGGNOG_DB_PATH" ]]; then
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

    EGGNOG_DB_PATH="\$(find_database_file 'eggnog.db')"
    EGGNOG_DIAMOND_DB="\$(find_database_file '*.dmnd')"

    if [[ -z "\$EGGNOG_DB_PATH" || ! -s "\$EGGNOG_DB_PATH" ]]; then
        echo "ERROR: eggnog.db was not found under: \$EGGNOG_DATA_DIR" >> "\$STATUS"
        exit 1
    fi

    EGGNOG_DB_PATH="\$(readlink -f "\$EGGNOG_DB_PATH")"

    if [[ "${params.eggnog_method}" == "diamond" ]]; then
        if [[ -z "\$EGGNOG_DIAMOND_DB" || ! -s "\$EGGNOG_DIAMOND_DB" ]]; then
            echo "ERROR: A DIAMOND database was not found under: \$EGGNOG_DATA_DIR" >> "\$STATUS"
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

    echo "Final EggNOG data directory: \$EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "Final EggNOG SQLite database: \$EGGNOG_DB_PATH" >> "\$STATUS"
    echo "Final EggNOG DIAMOND database: \${EGGNOG_DIAMOND_DB:-not used}" >> "\$STATUS"
    echo "Final EggNOG MMseqs database: \${MMSEQS_DB:-not used}" >> "\$STATUS"
    echo "EGGNOG_DATA_DIR=\$EGGNOG_DATA_DIR" >> "\$STATUS"
    echo "EGGNOG_DB_PATH=\$EGGNOG_DB_PATH" >> "\$STATUS"
    echo "EGGNOG_DIAMOND_DB=\$EGGNOG_DIAMOND_DB" >> "\$STATUS"
    echo "EGGNOG_MMSEQS_DB=\$MMSEQS_DB" >> "\$STATUS"
    
    echo "EggNOG-mapper setup finished: \$(date)" >> "\$STATUS"
    """
}

process RUN_EGGNOG {

    tag "eggnog_assembly_scaffolds"

    publishDir "${params.outdir}/eggnog", mode: params.publish_tool_outputs_mode, pattern: "eggnog_out/*", saveAs: { filename ->
        filename.tokenize('/').last()
    }

    publishDir "${params.outdir}/logs", mode: "copy", pattern: "eggnog.log"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "eggnog_status.tsv"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "eggnog_scaffold_manifest.tsv"

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
    path "eggnog_status.tsv", emit: status

    path "eggnog.log", emit: log_file

    path "eggnog_scaffold_manifest.tsv", emit: scaffold_manifest

    path "eggnog_out/*", emit: eggnog_out

    script:
    """
    set -euo pipefail

    LOG="eggnog.log"

    EGGNOG_ENV="\$(grep '^EGGNOG_ENV=' "${setup_status}" \
        | tail -n 1 \
        | cut -d= -f2- || true)"

    EGGNOG_DATA_DIR="\$(grep '^EGGNOG_DATA_DIR=' "${setup_status}" \
        | tail -n 1 \
        | cut -d= -f2- || true)"

    EGGNOG_DB_PATH="\$(grep '^EGGNOG_DB_PATH=' "${setup_status}" \
        | tail -n 1 \
        | cut -d= -f2- || true)"

    EGGNOG_MMSEQS_DB="\$(grep '^EGGNOG_MMSEQS_DB=' "${setup_status}" \
        | tail -n 1 \
        | cut -d= -f2- || true)"

    LOG="eggnog.log"
    : > "\$LOG"
    echo "EggNOG-mapper task initialized: \$(date)" >> "\$LOG"
    echo "Setup status input: ${setup_status}" >> "\$LOG"
    echo "Resolved EggNOG environment: \${EGGNOG_ENV:-not found}" >> "\$LOG"
    echo "Resolved EggNOG data directory: \${EGGNOG_DATA_DIR:-not found}" >> "\$LOG"
    echo "Resolved EggNOG SQLite database: \${EGGNOG_DB_PATH:-not found}" >> "\$LOG"

    if [[ -z "\$EGGNOG_ENV" || ! -d "\$EGGNOG_ENV" ]]; then
        echo "ERROR: Invalid EggNOG environment: \$EGGNOG_ENV" | tee -a "\$LOG" >&2
        exit 1
    fi

    if [[ -z "\$EGGNOG_DATA_DIR" || ! -d "\$EGGNOG_DATA_DIR" ]]; then
        echo "ERROR: Invalid EggNOG data directory: \$EGGNOG_DATA_DIR" | tee -a "\$LOG" >&2
        exit 1
    fi

    if [[ -z "\$EGGNOG_DB_PATH" || ! -s "\$EGGNOG_DB_PATH" ]]; then
        echo "ERROR: Invalid EggNOG SQLite database: \$EGGNOG_DB_PATH" | tee -a "\$LOG" >&2
        exit 1
    fi

    export PATH="\$EGGNOG_ENV/bin:\$PATH"
    export EGGNOG_DATA_DIR
    export EGGNOG_DATA_PATH="\$EGGNOG_DATA_DIR"

    echo "EggNOG-mapper started: \$(date)" >> "\$LOG"
    echo "Input FASTA: ${combined_fasta}" >> "\$LOG"
    echo "Input scaffold manifest: ${input_scaffold_manifest}" >> "\$LOG"
    echo "EggNOG environment: \$EGGNOG_ENV" >> "\$LOG"
    echo "EggNOG data directory: \$EGGNOG_DATA_DIR" >> "\$LOG"
    echo "EggNOG SQLite database: \$EGGNOG_DB_PATH" >> "\$LOG"
    echo "EggNOG MMseqs database: \${EGGNOG_MMSEQS_DB:-not used}" >> "\$LOG"
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

    printf 'tool\\tstatus\\texit_status\\tminimum_scaffold_bp\\tinput_scaffolds\\tinput_bp\\toutput_dir\\tmessage\\n' \
        > eggnog_status.tsv

    if [[ ! -s "${combined_fasta}" ]]; then
        echo "ERROR: Combined EggNOG input FASTA is missing or empty." \
            >> "\$LOG"

        printf 'eggnog\\tfailed\\t1\\t%s\\t0\\t0\\t%s\\tInput FASTA missing or empty\\n' \
            "${params.min_scaffold_bp}" \
            "${params.outdir}/eggnog" \
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

        printf 'eggnog\\tfailed\\t1\\t%s\\t0\\t0\\t%s\\tInput FASTA contains no records\\n' \
            "${params.min_scaffold_bp}" \
            "${params.outdir}/eggnog" \
            >> eggnog_status.tsv

        exit 1
    fi

    SEARCH_DB_ARGS=()

    if [[ "${params.eggnog_method}" == "mmseqs" ]]; then

        if [[ -z "\$EGGNOG_MMSEQS_DB" ||
              ! -s "\$EGGNOG_MMSEQS_DB" ]]; then

            echo "ERROR: MMseqs mode was selected, but no valid MMseqs database was found." \
                >> "\$LOG"

            echo "Expected MMseqs database:" >> "\$LOG"
            echo "  \$EGGNOG_DATA_DIR/mmseqs/mmseqs.db" >> "\$LOG"

            printf 'eggnog\\tfailed\\t1\\t%s\\t%s\\t%s\\t%s\\tMMseqs database missing\\n' \
                "${params.min_scaffold_bp}" \
                "\$INPUT_SCAFFOLDS" \
                "\$INPUT_BP" \
                "${params.outdir}/eggnog" \
                >> eggnog_status.tsv

            exit 1
        fi

        SEARCH_DB_ARGS=(
            --mmseqs_db
            "\$EGGNOG_MMSEQS_DB"
        )

    elif [[ "${params.eggnog_method}" == "diamond" ]]; then

        # DIAMOND uses the database discovered by EggNOG-mapper under
        # EGGNOG_DATA_DIR. No --db or --dmnd_db argument is supplied.

        SEARCH_DB_ARGS=()

    else

        echo "ERROR: Unsupported EggNOG method: ${params.eggnog_method}" \
            >> "\$LOG"

        echo "Supported methods are: diamond, mmseqs" >> "\$LOG"

        printf 'eggnog\\tfailed\\t1\\t%s\\t%s\\t%s\\t%s\\tUnsupported EggNOG method\\n' \
            "${params.min_scaffold_bp}" \
            "\$INPUT_SCAFFOLDS" \
            "\$INPUT_BP" \
            "${params.outdir}/eggnog" \
            >> eggnog_status.tsv

        exit 1
    fi

    echo "Running EggNOG-mapper." >> "\$LOG"

    echo "Search database arguments: \${SEARCH_DB_ARGS[*]:-none}" \
        >> "\$LOG"

    echo "Command:" >> "\$LOG"

    echo "emapper.py -m ${params.eggnog_method} --cpu ${task.cpus} -i ${combined_fasta} --itype ${params.eggnog_itype} --genepred ${params.eggnog_genepred} --trans_table ${params.eggnog_trans_table} --data_dir \$EGGNOG_DATA_DIR \${SEARCH_DB_ARGS[*]:-} --output ${params.eggnog_output_prefix} --output_dir eggnog_out --excel ${params.eggnog_extra_args}" \
        >> "\$LOG"

    set +e

    emapper.py \
        -m "${params.eggnog_method}" \
        --cpu "${task.cpus}" \
        -i "${combined_fasta}" \
        --itype "${params.eggnog_itype}" \
        --genepred "${params.eggnog_genepred}" \
        --trans_table "${params.eggnog_trans_table}" \
        --data_dir "\$EGGNOG_DATA_DIR" \
        "\${SEARCH_DB_ARGS[@]}" \
        --output "${params.eggnog_output_prefix}" \
        --output_dir eggnog_out \
        --excel \
        ${params.eggnog_extra_args} \
        >> "\$LOG" 2>&1

    EGGNOG_EXIT="\$?"

    set -e

    if [[ "\$EGGNOG_EXIT" -ne 0 ]]; then
        echo "ERROR: EggNOG-mapper failed with exit status \$EGGNOG_EXIT." \
            >> "\$LOG"

        if [[ "${params.eggnog_fail_nonfatal}" == "true" ]]; then
            printf 'eggnog\\tfailed_nonfatal\\t%s\\t%s\\t%s\\t%s\\t%s\\tEggNOG failed; workflow continued\\n' \
                "\$EGGNOG_EXIT" \
                "${params.min_scaffold_bp}" \
                "\$INPUT_SCAFFOLDS" \
                "\$INPUT_BP" \
                "${params.outdir}/eggnog" \
                >> eggnog_status.tsv

            exit 0
        fi

        printf 'eggnog\\tfailed\\t%s\\t%s\\t%s\\t%s\\t%s\\tEggNOG-mapper failed\\n' \
            "\$EGGNOG_EXIT" \
            "${params.min_scaffold_bp}" \
            "\$INPUT_SCAFFOLDS" \
            "\$INPUT_BP" \
            "${params.outdir}/eggnog" \
            >> eggnog_status.tsv

        exit "\$EGGNOG_EXIT"
    fi

    OUTPUT_FILES="\$(find eggnog_out -type f | wc -l | tr -d ' ')"

    if [[ "\$OUTPUT_FILES" -eq 0 ]]; then
        echo "ERROR: EggNOG completed but produced no output files." \
            >> "\$LOG"

        printf 'eggnog\\tfailed\\t1\\t%s\\t%s\\t%s\\t%s\\tEggNOG produced no output files\\n' \
            "${params.min_scaffold_bp}" \
            "\$INPUT_SCAFFOLDS" \
            "\$INPUT_BP" \
            "${params.outdir}/eggnog" \
            >> eggnog_status.tsv

        exit 1
    fi

    printf 'eggnog\\tcompleted\\t0\\t%s\\t%s\\t%s\\t%s\\tEggNOG completed; output_files=%s\\n' \
        "${params.min_scaffold_bp}" \
        "\$INPUT_SCAFFOLDS" \
        "\$INPUT_BP" \
        "${params.outdir}/eggnog" \
        "\$OUTPUT_FILES" \
        >> eggnog_status.tsv

    echo "EggNOG output files: \$OUTPUT_FILES" >> "\$LOG"
    echo "EggNOG output listing:" >> "\$LOG"

    find eggnog_out \
        -maxdepth 2 \
        -type f \
        -printf '%p\\n' \
        | sort \
        >> "\$LOG"

    echo "EggNOG-mapper finished: \$(date)" >> "\$LOG"
    """
}

/*
 * Write a combined run summary.
 */
process WRITE_ANNOTATION_SUMMARY {

    tag "write_assembly_annotation_summary"

    publishDir "${params.outdir}/summary", mode: "copy", pattern: "assembly_annotation_run_summary.tsv"

    input:
    path filter_stats
    path assembly_manifest
    path scaffold_manifest
    path eggnog_status

    output:
    path "assembly_annotation_run_summary.tsv", emit: summary

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
