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

params.working_dir = null
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
params.normalized_reads_dir = null
params.normalized_reads_compression = 4

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
params.mvp_version = "1.1.4"
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
params.ads_seqid = 30
params.functional_rdrp = false
params.rdrp_evalue = 0.01
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
    : (
        params.output_dir
            ? absPath(params.output_dir)
            : absPath(".")
    )

params.module1_outdir =
    "${params.results_dir}/module_1_readtrimming"

params.module2_outdir =
    "${params.results_dir}/module_2_readassembly"

params.module2b_outdir =
    "${params.results_dir}/module_2b_coassembly"

params.module5_outdir =
    "${params.results_dir}/module_5_subtractiveassembly"

params.outdir =
    "${params.results_dir}/AuxModule_2_mvp"

params.mvp_normalized_reads_dir =
    params.normalized_reads_dir
        ? absPath(params.normalized_reads_dir)
        : "${params.outdir}/normalized_reads"

/*
 * Parse a comma-, semicolon-, or whitespace-separated module list.
 */
def parseMvpModules(value) {
    def aliases = [
        "00"  : "0",
        "01"  : "1",
        "02"  : "2",
        "03"  : "3",
        "04"  : "4",
        "05"  : "5",
        "06"  : "6",
        "07"  : "7",
        "099" : "99"
    ]

    return value
        .toString()
        .split(/[,;\s]+/)
        .collect { item -> item.trim() }
        .findAll { item -> item }
        .collect { item -> aliases.containsKey(item) ? aliases[item] : item }
        .unique()
}


workflow {

    def supported_modules = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "99",
        "100"
    ]

    def canonical_order = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "99",
        "100"
    ]

    def selected_modules =
        parseMvpModules(params.mvp_modules)

    if (!selected_modules) {
        error(
            """
            No MVP modules were selected.

            Use, for example:

              --mvp_modules "0,1,2,3,6"
            """.stripIndent()
        )
    }

    def invalid_modules =
        selected_modules.findAll { module ->
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

    def ordered_modules =
        canonical_order.findAll { module ->
            module in selected_modules
        }

    if (!(params.read_type.toString() in ["short", "long"])) {
        error(
            """
            Invalid --read_type value: ${params.read_type}

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

    if (
        !(
            params.filtration.toString() in
            ["relaxed", "conservative"]
        )
    ) {
        error(
            """
            Invalid --filtration value: ${params.filtration}

            Supported values:
              relaxed
              conservative
            """.stripIndent()
        )
    }

    if (
        params.genomad_relaxed.toString().toBoolean() &&
        params.genomad_conservative.toString().toBoolean()
    ) {
        error(
            """
            --genomad_relaxed and --genomad_conservative cannot
            both be enabled.
            """.stripIndent()
        )
    }

    if ("99" in selected_modules) {

        if (!params.miuvig_identifier) {
            error(
                """
                MVP Module 99 was selected, but no MIUViG identifier
                was supplied.

                Provide:

                  --miuvig_identifier IDENTIFIER
                """.stripIndent()
            )
        }

        if (
            !params.miuvig_step ||
            !(
                params.miuvig_step.toString() in
                ["setup_metadata", "prep_submission"]
            )
        ) {
            error(
                """
                MVP Module 99 requires:

                  --miuvig_step setup_metadata

                or:

                  --miuvig_step prep_submission
                """.stripIndent()
            )
        }

        if (
            params.miuvig_step.toString() == "prep_submission" &&
            !params.miuvig_template
        ) {
            error(
                """
                MVP Module 99 prep_submission requires:

                  --miuvig_template /path/to/template.tsv
                """.stripIndent()
            )
        }
    }

    def include_individual =
        params.include_individual_assemblies
            .toString()
            .toBoolean()

    def include_coassemblies =
        params.include_coassemblies
            .toString()
            .toBoolean()

    def include_subtractive =
        params.include_subtractive_assemblies
            .toString()
            .toBoolean()

    if (
        !include_individual &&
        !include_coassemblies &&
        !include_subtractive
    ) {
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

    def assembly_manifest_file =
        params.assembly_manifest
            ? absPath(params.assembly_manifest)
            : "${params.module2_outdir}/summary/assembly_manifest.tsv"

    def trimmed_manifest_file =
        params.trimmed_manifest
            ? absPath(params.trimmed_manifest)
            : "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    def coassembly_manifest_file =
        params.coassembly_manifest
            ? absPath(params.coassembly_manifest)
            : "${params.module2b_outdir}/summary/assembly_manifest.tsv"

    def coassembly_reads_file =
        params.coassembly_trimmed_manifest
            ? absPath(params.coassembly_trimmed_manifest)
            : "${params.module2b_outdir}/summary/coassembly_trimmed_manifest.tsv"

    def subtractive_manifest_file =
        params.subtractive_assembly_manifest
            ? absPath(params.subtractive_assembly_manifest)
            : "${params.module5_outdir}/summary/subtractive_assembly_manifest.tsv"

    def subtractive_assembly_dir =
        params.subtractive_assembly_dir
            ? absPath(params.subtractive_assembly_dir)
            : "${params.module5_outdir}/assemblies"

    def subtractive_unmapped_reads_dir =
        params.subtractive_unmapped_reads_dir
            ? absPath(params.subtractive_unmapped_reads_dir)
            : "${params.module5_outdir}/unmapped_reads"

    log.info(
        "Auxiliary Module 2 MVP output directory: " +
        "${params.outdir}"
    )

    log.info(
        "Selected MVP modules: " +
        "${ordered_modules.join(',')}"
    )

    log.info(
        "MVP modules will execute in this order: " +
        "${ordered_modules.join(' -> ')}"
    )

    log.info(
        "Include Module 2 individual/rarefied assemblies: " +
        "${include_individual}"
    )

    log.info(
        "Include Module 2b coassemblies: " +
        "${include_coassemblies}"
    )

    log.info(
        "Include Module 5 subtractive assemblies: " +
        "${include_subtractive}"
    )

    log.info(
        "Module 2 assembly manifest: " +
        "${assembly_manifest_file}"
    )

    log.info(
        "Module 1 trimmed-read manifest: " +
        "${trimmed_manifest_file}"
    )

    log.info(
        "Module 2b assembly manifest: " +
        "${coassembly_manifest_file}"
    )

    log.info(
        "Module 2b trimmed-read manifest: " +
        "${coassembly_reads_file}"
    )

    log.info(
        "Module 5 subtractive assembly manifest: " +
        "${subtractive_manifest_file}"
    )

    log.info(
        "Module 5 subtractive assembly directory: " +
        "${subtractive_assembly_dir}"
    )

    log.info(
        "Module 5 unmapped-read directory: " +
        "${subtractive_unmapped_reads_dir}"
    )

    log.info("MVP version: ${params.mvp_version}")
    log.info("MVP threads: ${params.threads}")
    log.info("MVP read type: ${params.read_type}")

    if (!("0" in selected_modules)) {
        log.warn(
            "MVP Module 00 was not selected. Its required setup " +
            "and directory outputs must already exist."
        )
    }

    if (
        ("2" in selected_modules) &&
        !("1" in selected_modules)
    ) {
        log.warn(
            "MVP Module 02 was selected without Module 01. " +
            "Module 01 outputs must already exist."
        )
    }

    if (
        ("3" in selected_modules) &&
        !("2" in selected_modules)
    ) {
        log.warn(
            "MVP Module 03 was selected without Module 02. " +
            "Module 02 outputs must already exist."
        )
    }

    if (
        ("4" in selected_modules) &&
        !("3" in selected_modules)
    ) {
        log.warn(
            "MVP Module 04 was selected without Module 03. " +
            "Module 03 outputs and the mapping index must already exist."
        )
    }

    if (
        ("5" in selected_modules) &&
        !("4" in selected_modules)
    ) {
        log.warn(
            "MVP Module 05 was selected without Module 04. " +
            "Module 04 coverage outputs must already exist."
        )
    }

    if (
        ("6" in selected_modules) &&
        !("3" in selected_modules)
    ) {
        log.warn(
            "MVP Module 06 was selected without Module 03. " +
            "Module 03 protein and sequence outputs must already exist."
        )
    }

    if (
        ("7" in selected_modules) &&
        !("3" in selected_modules)
    ) {
        log.warn(
            "MVP Module 07 was selected without Module 03. " +
            "Module 03 outputs must already exist."
        )
    }

    log.info(
        "MVP normalized Module 1 read directory: " +
        "${params.mvp_normalized_reads_dir}"
    )

    log.info(
        "All MVP read inputs will be treated as interleaved."
    )

    NORMALIZE_MVP_READS(
        channel.value(include_individual),
        channel.value(trimmed_manifest_file)
    )

    SETUP_AUXMODULE2_MVP()

    PREPARE_MVP_METADATA(
        channel.value(include_individual),
        channel.value(include_coassemblies),
        channel.value(include_subtractive),
        channel.value(assembly_manifest_file),
        NORMALIZE_MVP_READS.out.normalized_manifest,
        channel.value(coassembly_manifest_file),
        channel.value(coassembly_reads_file),
        channel.value(subtractive_manifest_file),
        channel.value(subtractive_assembly_dir),
        channel.value(subtractive_unmapped_reads_dir),
        NORMALIZE_MVP_READS.out.status
    )

    RUN_MVP(
        PREPARE_MVP_METADATA.out.metadata,
        PREPARE_MVP_METADATA.out.input_summary,
        SETUP_AUXMODULE2_MVP.out.status
    )
}

/*
 * Normalize all Module 1 reads to interleaved FASTQ for MVP.
 *
 * Paired R1/R2 files are interleaved record-by-record. Existing
 * interleaved files are validated and recompressed into the same
 * persistent normalized-read directory.
 *
 * Module 2b coassembly reads and Module 5 subtractive reads are already
 * interleaved and are handled later by PREPARE_MVP_METADATA.
 */
 
process NORMALIZE_MVP_READS {

    tag "normalize_mvp_reads_to_interleaved"

    publishDir "${params.outdir}/metadata",
        mode: "copy",
        pattern: "mvp_normalized_trimmed_manifest.tsv"

    publishDir "${params.outdir}/summary",
        mode: "copy",
        pattern: "mvp_read_normalization_stats.tsv"

    publishDir "${params.outdir}/logs",
        mode: "copy",
        pattern: "normalize_mvp_reads.log"

    input:
    val include_individual
    val trimmed_manifest

    output:
    path "mvp_normalized_trimmed_manifest.tsv",
        emit: normalized_manifest

    path "mvp_read_normalization_stats.tsv",
        emit: stats

    path "mvp_read_normalization_status.env",
        emit: status

    path "normalize_mvp_reads.log",
        emit: log_file

    script:
    """
    set -euo pipefail

    LOG="normalize_mvp_reads.log"
    NORMALIZED_MANIFEST="mvp_normalized_trimmed_manifest.tsv"
    STATS="mvp_read_normalization_stats.tsv"
    STATUS="mvp_read_normalization_status.env"

    NORMALIZED_DIR="${params.mvp_normalized_reads_dir}"

    echo "MVP read normalization started: \$(date)" > "\$LOG"
    echo "Include Module 2 assemblies: ${include_individual}" >> "\$LOG"
    echo "Input Module 1 manifest: ${trimmed_manifest}" >> "\$LOG"
    echo "Persistent normalized-read directory: \$NORMALIZED_DIR" >> "\$LOG"
    echo "Compression level: ${params.normalized_reads_compression}" >> "\$LOG"
    echo "Launch directory: ${workflow.launchDir}" >> "\$LOG"
    echo "Task directory: \$(pwd -P)" >> "\$LOG"
    echo "----------------------------------------" >> "\$LOG"

    mkdir -p "\$NORMALIZED_DIR"

    python3 - \\
        "${include_individual}" \\
        "${trimmed_manifest}" \\
        "\$NORMALIZED_DIR" \\
        "${params.normalized_reads_compression}" \\
        "${workflow.launchDir}" \\
        "\$NORMALIZED_MANIFEST" \\
        "\$STATS" \\
        "\$LOG" <<'PY'
import csv
import gzip
import os
import re
import sys
import tempfile
from itertools import zip_longest
from pathlib import Path


(
    include_individual,
    trimmed_manifest,
    normalized_dir,
    compression_level,
    launch_dir,
    output_manifest,
    output_stats,
    log_file,
) = sys.argv[1:]


include_individual = include_individual.lower() == "true"
trimmed_manifest = Path(trimmed_manifest)
normalized_dir = Path(normalized_dir).resolve()
compression_level = int(compression_level)
launch_dir = Path(launch_dir).resolve()
output_manifest = Path(output_manifest)
output_stats = Path(output_stats)
log_file = Path(log_file)


if not 0 <= compression_level <= 9:
    raise RuntimeError(
        "normalized_reads_compression must be between 0 and 9"
    )


normalized_dir.mkdir(
    parents=True,
    exist_ok=True,
)


def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)


def safe_name(value):
    value = str(value or "").strip()
    value = re.sub(
        r"[^A-Za-z0-9._-]+",
        "_",
        value,
    )
    value = value.strip("._-")

    return value or "unnamed"


def resolve_path(value):
    value = str(value or "").strip()

    if not value:
        return None

    path = Path(value)

    if not path.is_absolute():
        path = launch_dir / path

    return path.resolve()


def open_fastq(path):
    path = Path(path)

    if path.name.endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            errors="strict",
        )

    return path.open(
        "rt",
        errors="strict",
    )


def read_fastq_records(path):
    with open_fastq(path) as handle:
        record_number = 0

        while True:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            separator = handle.readline()
            quality = handle.readline()

            record_number += 1

            if not sequence or not separator or not quality:
                raise RuntimeError(
                    f"Incomplete FASTQ record {record_number} "
                    f"in {path}"
                )

            if not header.startswith("@"):
                raise RuntimeError(
                    f"Invalid FASTQ header at record "
                    f"{record_number} in {path}: "
                    f"{header.rstrip()}"
                )

            if not separator.startswith("+"):
                raise RuntimeError(
                    f"Invalid FASTQ separator at record "
                    f"{record_number} in {path}: "
                    f"{separator.rstrip()}"
                )

            sequence_length = len(
                sequence.rstrip("\\r\\n")
            )

            quality_length = len(
                quality.rstrip("\\r\\n")
            )

            if sequence_length != quality_length:
                raise RuntimeError(
                    f"Sequence/quality length mismatch at "
                    f"record {record_number} in {path}"
                )

            yield (
                header,
                sequence,
                separator,
                quality,
            )


def temporary_output(destination):
    descriptor, name = tempfile.mkstemp(
        prefix=f".{destination.name}.",
        suffix=".tmp.fastq.gz",
        dir=str(destination.parent),
    )

    os.close(descriptor)

    return Path(name)


def validate_interleaved_fastq(path):
    record_count = 0

    for _ in read_fastq_records(path):
        record_count += 1

    if record_count % 2 != 0:
        raise RuntimeError(
            f"Interleaved FASTQ contains an odd number "
            f"of records ({record_count}): {path}"
        )

    return record_count, record_count // 2


def interleave_paired_fastq(
    read1,
    read2,
    destination,
):
    temporary = temporary_output(destination)
    pair_count = 0

    try:
        with gzip.open(
            temporary,
            "wt",
            compresslevel=compression_level,
        ) as output:
            iterator1 = read_fastq_records(read1)
            iterator2 = read_fastq_records(read2)

            for pair_number, pair in enumerate(
                zip_longest(
                    iterator1,
                    iterator2,
                ),
                start=1,
            ):
                record1, record2 = pair

                if record1 is None or record2 is None:
                    raise RuntimeError(
                        f"Paired FASTQ files have unequal "
                        f"record counts at pair {pair_number}: "
                        f"{read1}, {read2}"
                    )

                output.writelines(record1)
                output.writelines(record2)
                pair_count += 1

        record_count, validated_pairs = (
            validate_interleaved_fastq(temporary)
        )

        if validated_pairs != pair_count:
            raise RuntimeError(
                f"Internal pair-count mismatch while "
                f"interleaving {read1} and {read2}"
            )

        os.replace(
            temporary,
            destination,
        )

    except Exception:
        if temporary.exists():
            temporary.unlink()

        raise

    return record_count, pair_count


def normalize_existing_interleaved(
    source,
    destination,
):
    temporary = temporary_output(destination)
    record_count = 0

    try:
        with gzip.open(
            temporary,
            "wt",
            compresslevel=compression_level,
        ) as output:
            for record in read_fastq_records(source):
                output.writelines(record)
                record_count += 1

        if record_count % 2 != 0:
            raise RuntimeError(
                f"Interleaved FASTQ contains an odd "
                f"number of records ({record_count}): "
                f"{source}"
            )

        validated_records, pair_count = (
            validate_interleaved_fastq(temporary)
        )

        if validated_records != record_count:
            raise RuntimeError(
                f"Internal record-count mismatch while "
                f"normalizing {source}"
            )

        os.replace(
            temporary,
            destination,
        )

    except Exception:
        if temporary.exists():
            temporary.unlink()

        raise

    return record_count, pair_count


manifest_fields = [
    "sample_id",
    "safe_sample_id",
    "layout",
    "read1",
    "read2",
    "interleaved",
    "merged",
    "fastp_html",
    "fastp_json",
]

stats_fields = [
    "sample_id",
    "safe_sample_id",
    "source_layout",
    "source_read1",
    "source_read2",
    "source_interleaved",
    "normalized_interleaved",
    "fastq_records",
    "pairs_or_fragments",
    "normalization_action",
]

normalized_rows = []
stats_rows = []
used_output_names = set()


if include_individual:
    if not trimmed_manifest.is_file():
        raise RuntimeError(
            f"Required Module 1 trimmed manifest does not exist: "
            f"{trimmed_manifest}"
        )

    with trimmed_manifest.open() as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\\t",
        )

        available_fields = set(
            reader.fieldnames or []
        )

        missing = set(manifest_fields) - available_fields

        if missing:
            raise RuntimeError(
                f"Module 1 trimmed manifest is missing "
                f"required columns: {sorted(missing)}"
            )

        for row in reader:
            if not any(
                str(value or "").strip()
                for value in row.values()
            ):
                continue

            sample_id = str(
                row.get("sample_id", "")
            ).strip()

            safe_sample_id = safe_name(
                row.get("safe_sample_id", "")
                or sample_id
            )

            source_layout = str(
                row.get("layout", "")
            ).strip()

            if not sample_id:
                raise RuntimeError(
                    "Module 1 trimmed manifest contains "
                    "an empty sample_id"
                )

            output_name = (
                f"{safe_sample_id}."
                f"mvp_interleaved.fastq.gz"
            )

            if output_name in used_output_names:
                raise RuntimeError(
                    f"Multiple Module 1 rows would create "
                    f"the same normalized FASTQ: "
                    f"{output_name}"
                )

            used_output_names.add(output_name)

            destination = (
                normalized_dir / output_name
            )

            if source_layout == "paired":
                read1 = resolve_path(
                    row.get("read1", "")
                )

                read2 = resolve_path(
                    row.get("read2", "")
                )

                if read1 is None or read2 is None:
                    raise RuntimeError(
                        f"Paired sample {sample_id} is "
                        f"missing read1 or read2"
                    )

                if not read1.is_file():
                    raise RuntimeError(
                        f"Read 1 does not exist for "
                        f"{sample_id}: {read1}"
                    )

                if not read2.is_file():
                    raise RuntimeError(
                        f"Read 2 does not exist for "
                        f"{sample_id}: {read2}"
                    )

                records, pairs = (
                    interleave_paired_fastq(
                        read1,
                        read2,
                        destination,
                    )
                )

                action = "paired_to_interleaved"
                source_read1 = str(read1)
                source_read2 = str(read2)
                source_interleaved = ""

            elif source_layout == "interleaved":
                source_interleaved_path = resolve_path(
                    row.get("interleaved", "")
                )

                if source_interleaved_path is None:
                    raise RuntimeError(
                        f"Interleaved sample {sample_id} "
                        f"has no interleaved read path"
                    )

                if not source_interleaved_path.is_file():
                    raise RuntimeError(
                        f"Interleaved reads do not exist "
                        f"for {sample_id}: "
                        f"{source_interleaved_path}"
                    )

                records, pairs = (
                    normalize_existing_interleaved(
                        source_interleaved_path,
                        destination,
                    )
                )

                action = (
                    "normalized_existing_interleaved"
                )

                source_read1 = ""
                source_read2 = ""
                source_interleaved = str(
                    source_interleaved_path
                )

            else:
                raise RuntimeError(
                    f"Unsupported Module 1 read layout "
                    f"for {sample_id}: {source_layout}"
                )

            if records == 0:
                raise RuntimeError(
                    f"Normalized FASTQ contains zero "
                    f"records for {sample_id}: "
                    f"{destination}"
                )

            normalized_rows.append(
                {
                    "sample_id": sample_id,
                    "safe_sample_id": safe_sample_id,
                    "layout": "interleaved",
                    "read1": "",
                    "read2": "",
                    "interleaved": str(
                        destination.resolve()
                    ),
                    "merged": str(
                        row.get("merged", "")
                    ).strip(),
                    "fastp_html": str(
                        row.get("fastp_html", "")
                    ).strip(),
                    "fastp_json": str(
                        row.get("fastp_json", "")
                    ).strip(),
                }
            )

            stats_rows.append(
                {
                    "sample_id": sample_id,
                    "safe_sample_id": safe_sample_id,
                    "source_layout": source_layout,
                    "source_read1": source_read1,
                    "source_read2": source_read2,
                    "source_interleaved": (
                        source_interleaved
                    ),
                    "normalized_interleaved": str(
                        destination.resolve()
                    ),
                    "fastq_records": records,
                    "pairs_or_fragments": pairs,
                    "normalization_action": action,
                }
            )

            log(
                f"Normalized sample {sample_id}: "
                f"source_layout={source_layout}, "
                f"records={records}, "
                f"pairs_or_fragments={pairs}, "
                f"output={destination}"
            )


with output_manifest.open("w") as output:
    writer = csv.DictWriter(
        output,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=manifest_fields,
    )

    writer.writeheader()
    writer.writerows(normalized_rows)


with output_stats.open("w") as output:
    writer = csv.DictWriter(
        output,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=stats_fields,
    )

    writer.writeheader()
    writer.writerows(stats_rows)


log("----------------------------------------")
log(
    f"Normalized Module 1 samples: "
    f"{len(normalized_rows)}"
)
log(
    f"Normalized-read directory: "
    f"{normalized_dir}"
)
PY

    printf 'NORMALIZED_READS_DIR=%s\\n' \
        "\$NORMALIZED_DIR" \
        > "\$STATUS"

    printf 'NORMALIZED_MANIFEST=%s\\n' \
        "\$(pwd -P)/\$NORMALIZED_MANIFEST" \
        >> "\$STATUS"

    printf 'NORMALIZED_SAMPLE_COUNT=%s\\n' \
        "\$(tail -n +2 "\$NORMALIZED_MANIFEST" |
            awk 'NF > 0' |
            wc -l |
            tr -d ' ')" \
        >> "\$STATUS"

    printf 'STATUS=complete\\n' \
        >> "\$STATUS"

    echo "MVP read normalization finished: \$(date)" \
        >> "\$LOG"
    """
}

/*
 * Install and validate MVP.
 */
process SETUP_AUXMODULE2_MVP {

    tag "setup_auxmodule2_mvp"

    publishDir "${params.outdir}/setup",
        mode: "copy",
        pattern: "auxmodule2_mvp_tools_status.env"

    output:
    path "auxmodule2_mvp_tools_status.env",
        emit: status

    script:
    def env_dir = params.tool_env_dir
        ? absPath(params.tool_env_dir)
        : absPath(
            "${params.outdir}/conda_envs/mvp"
        )

    """
    set -euo pipefail

    STATUS_FILE="auxmodule2_mvp_tools_status.env"
    TOOL_ENV="${env_dir}"

    echo "Auxiliary Module 2 MVP setup started: \$(date)" \
        > "\$STATUS_FILE"

    echo "Requested MVP version: ${params.mvp_version}" \
        >> "\$STATUS_FILE"

    echo "Requested tool environment: \$TOOL_ENV" \
        >> "\$STATUS_FILE"

    echo "----------------------------------------" \
        >> "\$STATUS_FILE"

    check_mvp() {
        local prefix="\$1"

        if [[ "\$prefix" == "SYSTEM" ]]; then
            if ! command -v mvip >/dev/null 2>&1; then
                return 1
            fi

            mvip -h \
                >> "\$STATUS_FILE" 2>&1 ||
                return 1

            return 0
        fi

        if [[ ! -x "\$prefix/bin/mvip" ]]; then
            return 1
        fi

        "\$prefix/bin/mvip" -h \
            >> "\$STATUS_FILE" 2>&1 ||
            return 1

        return 0
    }

    if [[ -d "\$TOOL_ENV" ]]; then
        echo "Existing MVP environment detected: \$TOOL_ENV" \
            >> "\$STATUS_FILE"

        if check_mvp "\$TOOL_ENV"; then
            echo "Existing MVP environment passed validation." \
                >> "\$STATUS_FILE"

            echo "TOOL_ENV=\$TOOL_ENV" \
                >> "\$STATUS_FILE"

            echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
                >> "\$STATUS_FILE"

            exit 0
        fi

        echo "Existing MVP environment is incomplete or broken." \
            >> "\$STATUS_FILE"

        echo "Removing incomplete MVP environment." \
            >> "\$STATUS_FILE"

        rm -rf "\$TOOL_ENV"
    fi

    echo "Checking system/runtime MVP installation..." \
        >> "\$STATUS_FILE"

    if check_mvp "SYSTEM"; then
        echo "MVP is available from the system/runtime PATH." \
            >> "\$STATUS_FILE"

        echo "TOOL_ENV=SYSTEM" \
            >> "\$STATUS_FILE"

        echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
            >> "\$STATUS_FILE"

        exit 0
    fi

    echo "MVP is not available from the system/runtime PATH." \
        >> "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: MVP is unavailable and auto_install is false." \
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

    echo "Creating MVP environment: \$TOOL_ENV" \
        >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \
        -p "\$TOOL_ENV" \
        -c conda-forge \
        -c bioconda \
        "mvip=${params.mvp_version}" \
        >> "\$STATUS_FILE" 2>&1

    if ! check_mvp "\$TOOL_ENV"; then
        echo "ERROR: Newly created MVP environment failed validation." \
            >> "\$STATUS_FILE"

        exit 1
    fi

    echo "TOOL_ENV=\$TOOL_ENV" \
        >> "\$STATUS_FILE"

    echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
        >> "\$STATUS_FILE"
    """
}

process PREPARE_MVP_METADATA {

    tag "prepare_all_mvp_metadata"

    publishDir "${params.outdir}/metadata",
        mode: "copy",
        pattern: "mvp_metadata.tsv"

    publishDir "${params.outdir}/metadata",
        mode: "copy",
        pattern: "mvp_input_summary.tsv"

    publishDir "${params.outdir}/logs",
        mode: "copy",
        pattern: "prepare_mvp_metadata.log"

    input:
    val include_individual
    val include_coassemblies
    val include_subtractive

    val assembly_manifest
    path trimmed_manifest

    val coassembly_manifest
    val coassembly_trimmed_manifest

    val subtractive_manifest
    val subtractive_assembly_dir
    val subtractive_unmapped_reads_dir

    path normalization_status

    output:
    path "mvp_metadata.tsv",
        emit: metadata

    path "mvp_input_summary.tsv",
        emit: input_summary

    path "prepare_mvp_metadata.log",
        emit: log_file

    script:
    """
    set -euo pipefail

    if ! grep -q '^STATUS=complete\$' \
        "${normalization_status}"; then

        echo "ERROR: MVP read normalization did not complete successfully." \
            >&2

        cat "${normalization_status}" \
            >&2 || true

        exit 1
    fi

    if [[ "${include_individual}" == "true" ]]; then
        NORMALIZED_COUNT="\$(
            grep '^NORMALIZED_SAMPLE_COUNT=' \
                "${normalization_status}" |
            tail -n 1 |
            cut -d= -f2- ||
            true
        )"

        if [[ -z "\$NORMALIZED_COUNT" ||
              "\$NORMALIZED_COUNT" -eq 0 ]]; then

            echo "ERROR: Module 2 assemblies are enabled, but no " \
                 "normalized Module 1 read records were generated." \
                 >&2

            cat "${normalization_status}" \
                >&2 || true

            exit 1
        fi
    fi

    python3 - \\
        "${include_individual}" \\
        "${include_coassemblies}" \\
        "${include_subtractive}" \\
        "${assembly_manifest}" \\
        "${trimmed_manifest}" \\
        "${coassembly_manifest}" \\
        "${coassembly_trimmed_manifest}" \\
        "${subtractive_manifest}" \\
        "${subtractive_assembly_dir}" \\
        "${subtractive_unmapped_reads_dir}" \\
        "${workflow.launchDir}" \\
        "mvp_metadata.tsv" \\
        "mvp_input_summary.tsv" \\
        "prepare_mvp_metadata.log" <<'PY'
import csv
import re
import sys
from pathlib import Path

(
    include_individual,
    include_coassemblies,
    include_subtractive,
    assembly_manifest,
    trimmed_manifest,
    coassembly_manifest,
    coassembly_trimmed_manifest,
    subtractive_manifest,
    subtractive_assembly_dir,
    subtractive_unmapped_reads_dir,
    launch_dir,
    output_metadata,
    output_summary,
    log_file,
) = sys.argv[1:]


include_individual = include_individual.lower() == "true"
include_coassemblies = include_coassemblies.lower() == "true"
include_subtractive = include_subtractive.lower() == "true"

launch_dir = Path(launch_dir).resolve()

assembly_manifest = Path(assembly_manifest)
trimmed_manifest = Path(trimmed_manifest)

coassembly_manifest = Path(coassembly_manifest)
coassembly_trimmed_manifest = Path(
    coassembly_trimmed_manifest
)

subtractive_manifest = Path(subtractive_manifest)
subtractive_assembly_dir = Path(
    subtractive_assembly_dir
).resolve()

subtractive_unmapped_reads_dir = Path(
    subtractive_unmapped_reads_dir
).resolve()

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
        path = launch_dir / path

    return str(path.resolve())


def safe_name(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("._-")

    return value or "unnamed"


def read_tsv(path, required_columns, required=True):
    path = Path(path)

    if not path.is_file():
        if required:
            raise RuntimeError(
                f"Required manifest does not exist: {path}"
            )

        log(
            f"Optional manifest does not exist; skipping: "
            f"{path}"
        )

        return []

    with path.open() as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\\t",
        )

        fields = set(reader.fieldnames or [])
        missing = set(required_columns) - fields

        if missing:
            raise RuntimeError(
                f"Manifest {path} is missing required columns: "
                f"{sorted(missing)}"
            )

        return [
            row
            for row in reader
            if any(str(value or "").strip() for value in row.values())
        ]


def build_read_lookup(rows, label):
    lookup = {}

    for row in rows:
        sample_id = str(
            row.get("sample_id", "")
        ).strip()

        safe_sample_id = str(
            row.get("safe_sample_id", "")
        ).strip()

        layout = str(
            row.get("layout", "")
        ).strip()

        if not sample_id:
            raise RuntimeError(
                f"{label} contains an empty sample_id"
            )

        if layout == "paired":
            read_path = resolve_path(
                row.get("read1", "")
            )

            read2_path = resolve_path(
                row.get("read2", "")
            )

            if not read_path or not read2_path:
                raise RuntimeError(
                    f"Paired sample {sample_id} is missing "
                    f"read1 or read2"
                )

            if not Path(read_path).is_file():
                raise RuntimeError(
                    f"Read 1 does not exist for sample "
                    f"{sample_id}: {read_path}"
                )

            if not Path(read2_path).is_file():
                raise RuntimeError(
                    f"Read 2 does not exist for sample "
                    f"{sample_id}: {read2_path}"
                )

        elif layout == "interleaved":
            read_path = resolve_path(
                row.get("interleaved", "")
            )

            if not read_path:
                raise RuntimeError(
                    f"Interleaved sample {sample_id} has no "
                    f"read path"
                )

            if not Path(read_path).is_file():
                raise RuntimeError(
                    f"Interleaved reads do not exist for sample "
                    f"{sample_id}: {read_path}"
                )

        else:
            raise RuntimeError(
                f"Unsupported read layout for sample "
                f"{sample_id}: {layout}"
            )

        record = {
            "sample_id": sample_id,
            "safe_sample_id": safe_sample_id,
            "layout": layout,
            "read_path": read_path,
        }

        for key in {sample_id, safe_sample_id}:
            if not key:
                continue

            if key in lookup:
                previous = lookup[key]

                if previous["sample_id"] != sample_id:
                    raise RuntimeError(
                        f"Ambiguous read lookup key in "
                        f"{label}: {key}"
                    )

            lookup[key] = record

    return lookup


required_read_columns = [
    "sample_id",
    "safe_sample_id",
    "layout",
    "read1",
    "read2",
    "interleaved",
]

required_standard_assembly_columns = [
    "sample_id",
    "safe_sample_id",
    "assembly_sample_id",
    "assembler",
    "assembly_mode",
    "rarefaction_label",
    "assembly_strategy",
    "renamed_fasta",
]

required_subtractive_columns = [
    "sample_id",
    "safe_sample_id",
    "assembly_sample_id",
    "assembler",
    "assembly_mode",
    "renamed_fasta",
]


log("Preparing MVP metadata from all assembly modules")
log(f"Launch directory: {launch_dir}")
log(f"Include individual assemblies: {include_individual}")
log(f"Include coassemblies: {include_coassemblies}")
log(f"Include subtractive assemblies: {include_subtractive}")
log(f"Module 2 assembly manifest: {assembly_manifest}")
log(f"Module 1 trimmed manifest: {trimmed_manifest}")
log(f"Module 2b assembly manifest: {coassembly_manifest}")
log(
    "Module 2b trimmed manifest: "
    f"{coassembly_trimmed_manifest}"
)
log(
    "Module 5 subtractive manifest: "
    f"{subtractive_manifest}"
)
log(
    "Module 5 subtractive assembly directory: "
    f"{subtractive_assembly_dir}"
)
log(
    "Module 5 unmapped-read directory: "
    f"{subtractive_unmapped_reads_dir}"
)


ordinary_assemblies = []
ordinary_read_lookup = {}

if include_individual:
    ordinary_assemblies = read_tsv(
        assembly_manifest,
        required_standard_assembly_columns,
        required=True,
    )

    ordinary_read_rows = read_tsv(
        trimmed_manifest,
        required_read_columns,
        required=True,
    )

    ordinary_read_lookup = build_read_lookup(
        ordinary_read_rows,
        "Module 1 trimmed-read manifest",
    )


coassemblies = []
coassembly_read_lookup = {}

if include_coassemblies:
    coassembly_available = (
        coassembly_manifest.is_file() and
        coassembly_trimmed_manifest.is_file()
    )

    if coassembly_available:
        coassemblies = read_tsv(
            coassembly_manifest,
            required_standard_assembly_columns,
            required=True,
        )

        coassembly_read_rows = read_tsv(
            coassembly_trimmed_manifest,
            required_read_columns,
            required=True,
        )

        coassembly_read_lookup = build_read_lookup(
            coassembly_read_rows,
            "Module 2b coassembly read manifest",
        )
    else:
        log(
            "Module 2b outputs were not found. "
            "No coassemblies will be included."
        )


subtractive_assemblies = []

if include_subtractive:
    if subtractive_manifest.is_file():
        subtractive_assemblies = read_tsv(
            subtractive_manifest,
            required_subtractive_columns,
            required=True,
        )
    else:
        log(
            "Module 5 subtractive manifest was not found. "
            "No subtractive assemblies will be included."
        )


records = []
used_mvp_names = set()
skipped_assemblies = []


def add_record(
    row,
    source,
    assembly_path,
    read_path,
    read_layout,
):
    sample_id = str(
        row.get("sample_id", "")
    ).strip()

    safe_sample_id = str(
        row.get("safe_sample_id", "")
    ).strip()

    assembly_sample_id = str(
        row.get("assembly_sample_id", "")
    ).strip()

    assembler = str(
        row.get("assembler", "")
    ).strip()

    assembly_mode = str(
        row.get("assembly_mode", "")
    ).strip()

    rarefaction_label = str(
        row.get("rarefaction_label", "")
    ).strip()

    assembly_strategy = str(
        row.get("assembly_strategy", "")
    ).strip()

    assembly_path = resolve_path(assembly_path)
    read_path = resolve_path(read_path)

    if not sample_id:
        raise RuntimeError(
            f"An assembly record from {source} has no sample_id"
        )

    if not assembly_path:
        raise RuntimeError(
            f"Missing assembly path for {source} sample "
            f"{sample_id}"
        )

    assembly_file = Path(assembly_path)

    if not assembly_file.is_file():
        raise RuntimeError(
            f"Assembly file does not exist for {sample_id}: "
            f"{assembly_path}"
        )

    if assembly_file.stat().st_size == 0:
        skipped_assemblies.append(
            (
                source,
                sample_id,
                assembler,
                assembly_mode,
                assembly_path,
                "empty assembly file",
            )
        )

        log(
            f"WARNING: skipping empty assembly for "
            f"{sample_id}: {assembly_path}"
        )

        return

    if not read_path:
        raise RuntimeError(
            f"Missing read path for {source} sample "
            f"{sample_id}"
        )

    if not Path(read_path).is_file():
        raise RuntimeError(
            f"Read file does not exist for {source} sample "
            f"{sample_id}: {read_path}"
        )

    components = [
        source,
        safe_sample_id or sample_id,
        assembler or "assembler",
        assembly_mode or "assembly",
    ]

    if rarefaction_label:
        components.append(rarefaction_label)

    if assembly_strategy:
        components.append(assembly_strategy)

    base_name = safe_name("_".join(components))
    mvp_name = base_name
    suffix = 1

    while mvp_name in used_mvp_names:
        suffix += 1
        mvp_name = f"{base_name}_{suffix}"

    used_mvp_names.add(mvp_name)

    records.append(
        {
            "Sample_number": len(records) + 1,
            "Sample": mvp_name,
            "Assembly_Path": assembly_path,
            "Read_Path": read_path,
            "Variable": source,
            "source_sample_id": sample_id,
            "safe_sample_id": safe_sample_id,
            "assembly_sample_id": assembly_sample_id,
            "assembler": assembler,
            "assembly_mode": assembly_mode,
            "rarefaction_label": rarefaction_label,
            "assembly_strategy": assembly_strategy,
            "read_layout": read_layout,
        }
    )


for row in ordinary_assemblies:
    sample_id = str(
        row.get("sample_id", "")
    ).strip()

    safe_sample_id = str(
        row.get("safe_sample_id", "")
    ).strip()

    read_record = (
        ordinary_read_lookup.get(sample_id) or
        ordinary_read_lookup.get(safe_sample_id)
    )

    if read_record is None:
        raise RuntimeError(
            f"No Module 1 trimmed-read record matched "
            f"Module 2 assembly "
            f"{sample_id!r}/{safe_sample_id!r}"
        )

    add_record(
        row=row,
        source="individual",
        assembly_path=row.get("renamed_fasta", ""),
        read_path=read_record["read_path"],
        read_layout=read_record["layout"],
    )


for row in coassemblies:
    sample_id = str(
        row.get("sample_id", "")
    ).strip()

    safe_sample_id = str(
        row.get("safe_sample_id", "")
    ).strip()

    read_record = (
        coassembly_read_lookup.get(sample_id) or
        coassembly_read_lookup.get(safe_sample_id)
    )

    if read_record is None:
        raise RuntimeError(
            f"No Module 2b coassembly read record matched "
            f"coassembly {sample_id!r}/{safe_sample_id!r}"
        )

    add_record(
        row=row,
        source="coassembly",
        assembly_path=row.get("renamed_fasta", ""),
        read_path=read_record["read_path"],
        read_layout=read_record["layout"],
    )


for row in subtractive_assemblies:
    sample_id = str(
        row.get("sample_id", "")
    ).strip()

    safe_sample_id = str(
        row.get("safe_sample_id", "")
    ).strip()

    assembler = str(
        row.get("assembler", "")
    ).strip() or "megahit"

    expected_assembly = (
        subtractive_assembly_dir /
        f"{safe_sample_id}_{assembler}_subtractive.renamed.fa"
    )

    if not expected_assembly.is_file():
        candidates = sorted(
            subtractive_assembly_dir.glob(
                f"{safe_sample_id}_*_subtractive.renamed.fa"
            )
        )

        if len(candidates) == 1:
            expected_assembly = candidates[0]
        elif len(candidates) == 0:
            manifest_path = resolve_path(
                row.get("renamed_fasta", "")
            )

            if manifest_path and Path(manifest_path).is_file():
                expected_assembly = Path(manifest_path)
            else:
                raise RuntimeError(
                    f"Could not find subtractive assembly for "
                    f"{sample_id}: expected {expected_assembly}"
                )
        else:
            raise RuntimeError(
                f"Multiple subtractive assemblies matched "
                f"{safe_sample_id}: "
                f"{[str(path) for path in candidates]}"
            )

    unmapped_reads = (
        subtractive_unmapped_reads_dir /
        f"{safe_sample_id}.unmapped_interleaved.fastq.gz"
    )

    if not unmapped_reads.is_file():
        raise RuntimeError(
            f"Could not find subtractive unmapped reads for "
            f"{sample_id}: {unmapped_reads}"
        )

    add_record(
        row=row,
        source="subtractive",
        assembly_path=str(expected_assembly),
        read_path=str(unmapped_reads),
        read_layout="interleaved",
    )


if not records:
    raise RuntimeError(
        "No usable MVP input records were generated from "
        "Module 2, Module 2b, or Module 5."
    )


with output_metadata.open("w") as out:
    fields = [
        "Sample_number",
        "Sample",
        "Assembly_Path",
        "Read_Path",
        "Variable",
    ]

    writer = csv.DictWriter(
        out,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=fields,
    )

    writer.writeheader()

    for record in records:
        writer.writerow(
            {
                key: record[key]
                for key in fields
            }
        )


with output_summary.open("w") as out:
    fields = [
        "Sample_number",
        "Sample",
        "source_sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "read_layout",
        "Assembly_Path",
        "Read_Path",
        "Variable",
    ]

    writer = csv.DictWriter(
        out,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=fields,
    )

    writer.writeheader()
    writer.writerows(records)


source_counts = {
    source: sum(
        record["Variable"] == source
        for record in records
    )
    for source in [
        "individual",
        "coassembly",
        "subtractive",
    ]
}

read_layouts = sorted(
    {
        record["read_layout"]
        for record in records
    }
)

log(f"MVP input records written: {len(records)}")
log(
    "Individual/rarefied assembly records: "
    f"{source_counts['individual']}"
)
log(
    "Coassembly records: "
    f"{source_counts['coassembly']}"
)
log(
    "Subtractive assembly records: "
    f"{source_counts['subtractive']}"
)
log(
    "Read layouts represented: "
    f"{','.join(read_layouts)}"
)
log(
    "Empty assembly records skipped: "
    f"{len(skipped_assemblies)}"
)

for skipped in skipped_assemblies:
    log(
        "Skipped assembly: "
        + "\\t".join(str(value) for value in skipped)
    )
PY
    """
}

process RUN_MVP {

    tag "run_auxmodule2_mvp"

    stageInMode "symlink"

    /*
     * MVP writes directly into params.outdir.
     *
     * Do not add mvp_work as a publishDir output. The direct-write
     * behavior preserves MVP's required directory structure and retains
     * partial results if a later module fails.
     */
    publishDir "${params.outdir}/logs",
        mode: "copy",
        pattern: "mvp_commands.log"

    publishDir "${params.outdir}/logs",
        mode: "copy",
        pattern: "mvp_run.log"

    publishDir "${params.outdir}/summary",
        mode: "copy",
        pattern: "mvp_input_summary.tsv"

    publishDir "${params.outdir}/summary",
        mode: "copy",
        pattern: "mvp_complete.txt"

    cpus {
        params.threads as int
    }

    memory {
        def gb = params.memory_gb as int

        return gb > 0
            ? "${gb} GB"
            : null
    }

    input:

    path metadata
    path input_summary
    path tools_status

    output:

    path "mvp_commands.log",
        emit: commands

    path "mvp_run.log",
        emit: log_file

    path "mvp_input_summary.tsv",
        emit: input_summary

    path "mvp_complete.txt",
        emit: completion_marker

    script:

    def selected_modules =
        parseMvpModules(params.mvp_modules)

    def canonical_order = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "99",
        "100"
    ]

    def ordered_modules =
        canonical_order.findAll { module ->
            module in selected_modules
        }

    def run_module_0 =
        "0" in selected_modules

    def run_module_1 =
        "1" in selected_modules

    def run_module_2 =
        "2" in selected_modules

    def run_module_3 =
        "3" in selected_modules

    def run_module_4 =
        "4" in selected_modules

    def run_module_5 =
        "5" in selected_modules

    def run_module_6 =
        "6" in selected_modules

    def run_module_7 =
        "7" in selected_modules

    def run_module_99 =
        "99" in selected_modules

    def run_module_100 =
        "100" in selected_modules

    def install_databases =
        params.install_databases
            .toString()
            .toBoolean()

    def force =
        params.force
            .toString()
            .toBoolean()

    def mvp_work =
        absPath(params.outdir)

    def database_root = params.mvp_database_dir
        ? absPath(params.mvp_database_dir)
        : "${mvp_work}/databases"

    def genomad_database = params.genomad_db_path
        ? absPath(params.genomad_db_path)
        : "${database_root}/genomad_db"

    def checkv_database = params.checkv_db_path
        ? absPath(params.checkv_db_path)
        : "${database_root}/checkv_db"

    def module00_database_flag =
        install_databases
            ? ""
            : "--skip_install_databases"

    def skip_check_errors_argument =
        params.skip_check_errors
            .toString()
            .toBoolean()
            ? "--skip_check_errors"
            : ""

    def min_seq_argument =
        (params.min_seq_size as int) > 0
            ? "--min_seq_size ${params.min_seq_size}"
            : ""

    def genomad_filter_argument = ""

    if (
        params.genomad_relaxed
            .toString()
            .toBoolean()
    ) {
        genomad_filter_argument =
            "--genomad_relaxed"
    }
    else if (
        params.genomad_conservative
            .toString()
            .toBoolean()
    ) {
        genomad_filter_argument =
            "--genomad_conservative"
    }

    def skip_modify_headers_argument =
        params.skip_modify_headers
            .toString()
            .toBoolean()
            ? "--skip_modify_headers"
            : ""

    def unfiltered_protein_argument =
        params.unfiltered_protein_file
            .toString()
            .toBoolean()
            ? "--Unfiltered_protein_file"
            : ""

    def force_module01_argument =
        force
            ? "--force_genomad --force_checkv"
            : ""

    def interleaved_argument = "--interleaved"

    def mapping_delete_argument =
        params.delete_mapping_intermediates
            .toString()
            .toBoolean()
            ? "--delete_files"
            : ""

    def force_mapping_argument =
        force
            ? "--force_read_mapping"
            : ""

    def covered_fraction_argument =
        params.covered_fraction != null
            ? "--covered_fraction ${params.covered_fraction}"
            : ""

    def functional_ads_argument =
        params.functional_ads
            .toString()
            .toBoolean()
            ? "--ADS"
            : ""

    def functional_rdrp_argument =
        params.functional_rdrp
            .toString()
            .toBoolean()
            ? "--RdRP"
            : ""

    def functional_dram_argument =
        params.functional_dram
            .toString()
            .toBoolean()
            ? "--DRAM"
            : ""

    def functional_delete_argument =
        params.delete_functional_intermediates
            .toString()
            .toBoolean()
            ? "--delete_files"
            : ""

    def force_functional_argument =
        force
            ? "--force_prodigal --force_PHROGS " +
              "--force_PFAM --force_outputs"
            : ""

    if (
        force &&
        params.functional_ads
            .toString()
            .toBoolean()
    ) {
        force_functional_argument +=
            " --force_ADS"
    }

    if (
        force &&
        params.functional_rdrp
            .toString()
            .toBoolean()
    ) {
        force_functional_argument +=
            " --force_RdRP"
    }

    def keep_bam_argument =
        params.keep_bam
            .toString()
            .toBoolean()
            ? "--keep_bam"
            : ""

    def binning_delete_argument =
        params.delete_binning_intermediates
            .toString()
            .toBoolean()
            ? "--delete_files"
            : ""

    def binning_sample_group_argument =
        params.binning_sample_group != null
            ? "--binning_sample_group ${params.binning_sample_group}"
            : ""

    def read_mapping_sample_group_argument =
        params.read_mapping_sample_group != null
            ? "--read_mapping_sample_group ${params.read_mapping_sample_group}"
            : ""

    def force_binning_argument =
        force
            ? "--force_vrhyme --force_checkv " +
              "--force_read_mapping --force_outputs"
            : ""

    def summary_force_argument =
        force
            ? "--force"
            : ""

    def miuvig_template_argument =
        params.miuvig_template
            ? "-t ${absPath(params.miuvig_template)}"
            : ""

    """
    set -euo pipefail

    cp "${input_summary}" mvp_input_summary.tsv

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" |
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

    METADATA="\$(readlink -f "${metadata}")"

    MVP_WORK="${mvp_work}"

    GENOMAD_DB="${genomad_database}"
    CHECKV_DB="${checkv_database}"

    mkdir -p "\$MVP_WORK"
    mkdir -p "${database_root}"
    mkdir -p "\$(dirname "\$GENOMAD_DB")"
    mkdir -p "\$(dirname "\$CHECKV_DB")"

    : > mvp_commands.log
    : > mvp_run.log

    echo "Auxiliary Module 2 MVP started: \$(date)" \
        >> mvp_run.log

    echo "MVP executable: \$(command -v mvip)" \
        >> mvp_run.log

    echo "MVP working directory: \$MVP_WORK" \
        >> mvp_run.log

    echo "MVP metadata: \$METADATA" \
        >> mvp_run.log

    echo "Selected MVP modules: ${ordered_modules.join(',')}" \
        >> mvp_run.log

    echo "geNomad database: \$GENOMAD_DB" \
        >> mvp_run.log

    echo "CheckV database: \$CHECKV_DB" \
        >> mvp_run.log

    echo "Threads: ${task.cpus}" \
        >> mvp_run.log

    echo "----------------------------------------" \
        >> mvp_run.log

    run_command() {
        printf '%q ' "\$@" >> mvp_commands.log
        printf '\\n' >> mvp_commands.log

        echo "==================================================" \
            >> mvp_run.log

        echo "Started: \$(date)" \
            >> mvp_run.log

        printf 'Command: ' >> mvp_run.log
        printf '%q ' "\$@" >> mvp_run.log
        printf '\\n' >> mvp_run.log

        "\$@" >> mvp_run.log 2>&1

        echo "Finished: \$(date)" \
            >> mvp_run.log
    }

    if [[ "${run_module_0}" == "true" ]]; then

        run_command \
            mvip MVP_00_set_up_MVP \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${module00_database_flag} \
            ${skip_check_errors_argument} \
            --genomad_db_path "\$GENOMAD_DB" \
            --checkv_db_path "\$CHECKV_DB"
    fi

    if [[ "${run_module_1}" == "true" ]]; then

        run_command \
            mvip MVP_01_run_genomad_checkv \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --genomad_db_path "\$GENOMAD_DB" \
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
            --host_viral_genes_ratio \
                ${params.host_viral_genes_ratio}
    fi

    if [[ "${run_module_3}" == "true" ]]; then

        run_command \
            mvip MVP_03_do_clustering \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --min_ani ${params.min_ani} \
            --min_tcov ${params.min_tcov} \
            --min_qcov ${params.min_qcov} \
            --read_type ${params.read_type} \
            ${unfiltered_protein_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_4}" == "true" ]]; then

        run_command \
            mvip MVP_04_do_read_mapping \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --read_type ${params.read_type} \
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
            --host_viral_genes_ratio \
                ${params.host_viral_genes_ratio}
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
            --ADS_seqid ${params.ads_seqid} \
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
            --read_type ${params.read_type} \
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

    echo "==================================================" \
        >> mvp_run.log

    echo "Auxiliary Module 2 MVP finished: \$(date)" \
        >> mvp_run.log

    cat > mvp_complete.txt <<EOF
status=complete
finished=\$(date --iso-8601=seconds 2>/dev/null || date)
mvp_work=\$MVP_WORK
metadata=\$METADATA
selected_modules=${ordered_modules.join(',')}
EOF
    """
}