#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Module 1: Read trimming / quality control with fastp, followed by FastQC.
 *
 * Input:
 *   The read_manifest.tsv produced by module 0.
 *
 * This module auto-installs fastp and FastQC if missing, using mamba or conda.
 */

/*
 * Default parameters. Command-line values will override.
 */
params.input_manifest          = null
params.module0_outdir          = "./results/module_0_readprocess"

params.outdir                  = "./results/module_1_readtrimming"
params.output_dir              = null

params.fastp_version           = "0.23.4"
params.fastqc_version          = "0.12.1"

params.auto_install            = true
params.tool_env_dir            = null

params.fastp_threads           = 4
params.fastqc_threads          = 2
params.threads                 = null

/*
 * To avoid duplicating large trimmed FASTQ files, the default is symlink.
 *
 * Options:
 *   symlink  = lowest disk use, but results links depend on work/ remaining
 *   copy     = durable results, but duplicates large FASTQ files
 *   move     = generally not recommended here because trimmed reads are used downstream by FastQC
 */
params.publish_trimmed_mode    = "symlink"

/*
 * fastp compression level for gzip output.
 */
params.compression             = 4

/*
 * fastp trimming/filtering defaults.
 *
 * Adapter trimming is enabled by fastp by default.
 * For paired-end data, --detect_adapter_for_pe improves adapter detection.
 */
params.detect_adapter_for_pe   = true
params.enable_correction       = true

params.cut_front               = true
params.cut_tail                = true
params.cut_window_size         = 4
params.cut_mean_quality        = 20

params.qualified_quality_phred = 15
params.unqualified_percent     = 40
params.n_base_limit            = 5
params.length_required         = 15

params.trim_poly_g             = false
params.trim_poly_x             = true

workflow {
    def manifest_file = params.input_manifest ?: "${params.module0_outdir}/naming/read_manifest.tsv"

    def manifest_ch = channel.fromPath(
        manifest_file,
        type: 'file',
        checkIfExists: true
    )

    def samples_ch = manifest_ch
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row.sample_id.toString()
            def safe_id   = sample_id.replaceAll('[^A-Za-z0-9._-]+', '_')
            def layout    = row.layout.toString()

            if( layout == 'paired' ) {
                return tuple(
                    sample_id,
                    safe_id,
                    layout,
                    file(row.read1),
                    file(row.read2),
                    null
                )
            }
            else if( layout == 'interleaved' ) {
                return tuple(
                    sample_id,
                    safe_id,
                    layout,
                    null,
                    null,
                    file(row.interleaved)
                )
            }
            else {
                error "Unsupported layout in manifest for sample '${sample_id}': ${layout}"
            }
        }

    def paired_reads_ch = samples_ch
        .filter { _sample_id, _safe_id, layout, _read1, _read2, _interleaved ->
            layout == 'paired'
        }
        .map { sample_id, safe_id, _layout, read1, read2, _interleaved ->
            tuple(sample_id, safe_id, read1, read2)
        }

    def interleaved_reads_ch = samples_ch
        .filter { _sample_id, _safe_id, layout, _read1, _read2, _interleaved ->
            layout == 'interleaved'
        }
        .map { sample_id, safe_id, _layout, _read1, _read2, interleaved ->
            tuple(sample_id, safe_id, interleaved)
        }

    SETUP_MODULE1_TOOLS()

    FASTP_PAIRED(paired_reads_ch.combine(SETUP_MODULE1_TOOLS.out.status))
    FASTP_INTERLEAVED(interleaved_reads_ch.combine(SETUP_MODULE1_TOOLS.out.status))

    /*
     * Generate per-sample read count / trimming count summaries.
     */
    TRIMMING_STATS_PAIRED(FASTP_PAIRED.out.trimmed_reads)
    TRIMMING_STATS_INTERLEAVED(FASTP_INTERLEAVED.out.trimmed_reads)

    def all_trimming_stats_ch = TRIMMING_STATS_PAIRED.out.stats_file.mix(TRIMMING_STATS_INTERLEAVED.out.stats_file)

    WRITE_TRIMMING_STATS_SUMMARY(all_trimming_stats_ch.collect())

    /*
     * Build FastQC input channel from trimmed primary reads.
     *
     * Paired:
     *   Run FastQC on trimmed R1 and trimmed R2.
     *
     * Interleaved:
     *   Run FastQC on trimmed interleaved file.
     */
    def paired_fastqc_reads_ch = FASTP_PAIRED.out.trimmed_reads
        .flatMap { sample_id, _safe_id, read1_trimmed, read2_trimmed, _fastp_json ->
            return [
                tuple(sample_id, read1_trimmed),
                tuple(sample_id, read2_trimmed)
            ]
        }

    def interleaved_fastqc_reads_ch = FASTP_INTERLEAVED.out.trimmed_reads
        .map { sample_id, _safe_id, interleaved_trimmed, _fastp_json ->
            tuple(sample_id, interleaved_trimmed)
        }

    def trimmed_fastqc_input_ch = paired_fastqc_reads_ch.mix(interleaved_fastqc_reads_ch)

    RUN_FASTQC_TRIMMED(trimmed_fastqc_input_ch.combine(SETUP_MODULE1_TOOLS.out.status))

    /*
     * Write a public manifest pointing to published module 1 output paths.
     *
     * Each FASTP process writes one small manifest-record TSV.
     */
    def all_manifest_record_files_ch = FASTP_PAIRED.out.manifest_record.mix(FASTP_INTERLEAVED.out.manifest_record)

    WRITE_TRIMMED_MANIFEST(all_manifest_record_files_ch.collect())
}

process SETUP_MODULE1_TOOLS {
    tag "setup_fastp_fastqc"

    publishDir "${params.output_dir ?: params.outdir}/setup",
        mode: 'copy',
        pattern: "module1_tools_status.env"

    output:
    path "module1_tools_status.env", emit: status

    script:
    def base_outdir = params.output_dir ?: params.outdir
    def env_dir = params.tool_env_dir ?: "${base_outdir}/conda_envs/module1_tools"

    """
    set -euo pipefail

    STATUS_FILE="module1_tools_status.env"
    TOOL_ENV="${env_dir}"

    echo "Module 1 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested fastp version: ${params.fastp_version}" >> "\$STATUS_FILE"
    echo "Requested FastQC version: ${params.fastqc_version}" >> "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    if [[ -x "\$TOOL_ENV/bin/fastp" && -x "\$TOOL_ENV/bin/fastqc" ]]; then
        echo "Existing module-local environment detected." >> "\$STATUS_FILE"
        echo "ENV_DIR=\$TOOL_ENV" >> "\$STATUS_FILE"
        "\$TOOL_ENV/bin/fastp" --version >> "\$STATUS_FILE" 2>&1 || true
        "\$TOOL_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1 || true
        echo "Module 1 tool setup finished: \$(date)" >> "\$STATUS_FILE"
        exit 0
    fi

    if command -v fastp >/dev/null 2>&1 && command -v fastqc >/dev/null 2>&1; then
        echo "System/runtime fastp and FastQC detected." >> "\$STATUS_FILE"
        echo "fastp path: \$(command -v fastp)" >> "\$STATUS_FILE"
        echo "fastqc path: \$(command -v fastqc)" >> "\$STATUS_FILE"
        echo "ENV_DIR=SYSTEM" >> "\$STATUS_FILE"
        fastp --version >> "\$STATUS_FILE" 2>&1 || true
        fastqc --version >> "\$STATUS_FILE" 2>&1 || true
        echo "Module 1 tool setup finished: \$(date)" >> "\$STATUS_FILE"
        exit 0
    fi

    echo "fastp and/or FastQC not found in PATH." >> "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: Auto-install is disabled." >> "\$STATUS_FILE"
        echo "Install fastp/FastQC manually or rerun with --auto_install true." >> "\$STATUS_FILE"
        exit 1
    fi

    INSTALLER=""

    if command -v mamba >/dev/null 2>&1; then
        INSTALLER="mamba"
        echo "Using mamba: \$(command -v mamba)" >> "\$STATUS_FILE"
    elif command -v conda >/dev/null 2>&1; then
        INSTALLER="conda"
        echo "Using conda: \$(command -v conda)" >> "\$STATUS_FILE"
    else
        echo "ERROR: Neither mamba nor conda was found in PATH." >> "\$STATUS_FILE"
        echo "Please install mamba/conda or install fastp/FastQC manually." >> "\$STATUS_FILE"
        exit 1
    fi

    mkdir -p "\$(dirname "\$TOOL_ENV")"

    echo "Creating module-local environment:" >> "\$STATUS_FILE"
    echo "  \$TOOL_ENV" >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \\
        -p "\$TOOL_ENV" \\
        -c conda-forge \\
        -c bioconda \\
        "fastp=${params.fastp_version}" \\
        "fastqc=${params.fastqc_version}" \\
        >> "\$STATUS_FILE" 2>&1

    if [[ ! -x "\$TOOL_ENV/bin/fastp" ]]; then
        echo "ERROR: fastp was not found after installation." >> "\$STATUS_FILE"
        exit 1
    fi

    if [[ ! -x "\$TOOL_ENV/bin/fastqc" ]]; then
        echo "ERROR: FastQC was not found after installation." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "ENV_DIR=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "fastp installed at: \$TOOL_ENV/bin/fastp" >> "\$STATUS_FILE"
    echo "FastQC installed at: \$TOOL_ENV/bin/fastqc" >> "\$STATUS_FILE"

    "\$TOOL_ENV/bin/fastp" --version >> "\$STATUS_FILE" 2>&1 || true
    "\$TOOL_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1 || true

    echo "Module 1 tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}

process FASTP_PAIRED {
    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/trimmed_reads",
        mode: params.publish_trimmed_mode,
        pattern: "*.fastq.gz"

    publishDir "${params.output_dir ?: params.outdir}/fastp_reports",
        mode: 'copy',
        pattern: "*_fastp.*"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.fastp_threads as int
    }

    input:
    tuple val(sample_id), val(safe_id), path(read1), path(read2), path(tools_status)

    output:
    tuple val(sample_id),
          val(safe_id),
          path("${safe_id}_R1_trimmed.fastq.gz"),
          path("${safe_id}_R2_trimmed.fastq.gz"),
          path("${safe_id}_fastp.json"),
          emit: trimmed_reads

    path "${safe_id}_trimmed_manifest_record.tsv", emit: manifest_record

    path "${safe_id}_fastp.html", emit: html
    path "${safe_id}_fastp.log",  emit: log

    script:
    """
    set -euo pipefail

    TOOLS_ENV="\$(grep '^ENV_DIR=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOLS_ENV" && "\$TOOLS_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOLS_ENV/bin:\$PATH"
    fi

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH after setup." >&2
        echo "Tool setup status:" >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    FASTP_ARGS=()

    if [[ "${params.detect_adapter_for_pe.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--detect_adapter_for_pe)
    fi

    if [[ "${params.enable_correction.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--correction)
    fi

    if [[ "${params.cut_front.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--cut_front)
    fi

    if [[ "${params.cut_tail.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--cut_tail)
    fi

    if [[ "${params.trim_poly_g.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--trim_poly_g)
    fi

    if [[ "${params.trim_poly_x.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--trim_poly_x)
    fi

    echo "Running fastp paired-end trimming for sample: ${sample_id}"
    echo "Input R1: ${read1}"
    echo "Input R2: ${read2}"
    echo "Threads: ${task.cpus}"
    echo "NOTE: Read merging is disabled. Output will be trimmed R1 and trimmed R2 only."

    fastp \\
        -i "${read1}" \\
        -I "${read2}" \\
        -o "${safe_id}_R1_trimmed.fastq.gz" \\
        -O "${safe_id}_R2_trimmed.fastq.gz" \\
        --html "${safe_id}_fastp.html" \\
        --json "${safe_id}_fastp.json" \\
        --report_title "fastp report: ${sample_id}" \\
        --thread ${task.cpus} \\
        --compression ${params.compression} \\
        --cut_window_size ${params.cut_window_size} \\
        --cut_mean_quality ${params.cut_mean_quality} \\
        --qualified_quality_phred ${params.qualified_quality_phred} \\
        --unqualified_percent_limit ${params.unqualified_percent} \\
        --n_base_limit ${params.n_base_limit} \\
        --length_required ${params.length_required} \\
        "\${FASTP_ARGS[@]}" \\
        > "${safe_id}_fastp.log" 2>&1

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" \\
        "${safe_id}" \\
        "paired" \\
        "${params.output_dir ?: params.outdir}/trimmed_reads/${safe_id}_R1_trimmed.fastq.gz" \\
        "${params.output_dir ?: params.outdir}/trimmed_reads/${safe_id}_R2_trimmed.fastq.gz" \\
        "" \\
        "" \\
        "${params.output_dir ?: params.outdir}/fastp_reports/${safe_id}_fastp.html" \\
        "${params.output_dir ?: params.outdir}/fastp_reports/${safe_id}_fastp.json" \\
        > "${safe_id}_trimmed_manifest_record.tsv"
    """
}

process FASTP_INTERLEAVED {
    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/trimmed_reads",
        mode: params.publish_trimmed_mode,
        pattern: "*.fastq.gz"

    publishDir "${params.output_dir ?: params.outdir}/fastp_reports",
        mode: 'copy',
        pattern: "*_fastp.*"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.fastp_threads as int
    }

    input:
    tuple val(sample_id), val(safe_id), path(interleaved), path(tools_status)

    output:
    tuple val(sample_id),
          val(safe_id),
          path("${safe_id}_interleaved_trimmed.fastq.gz"),
          path("${safe_id}_fastp.json"),
          emit: trimmed_reads

    path "${safe_id}_trimmed_manifest_record.tsv", emit: manifest_record

    path "${safe_id}_fastp.html", emit: html
    path "${safe_id}_fastp.log",  emit: log

    script:
    """
    set -euo pipefail

    TOOLS_ENV="\$(grep '^ENV_DIR=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOLS_ENV" && "\$TOOLS_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOLS_ENV/bin:\$PATH"
    fi

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH after setup." >&2
        echo "Tool setup status:" >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    FASTP_ARGS=()

    # Do NOT add --detect_adapter_for_pe or --correction here.
    #
    # With interleaved input, fastp may try to open a missing/empty R2 filename
    # when --detect_adapter_for_pe is used, causing:
    #
    #   ERROR: Failed to open file:
    #
    # Those options are kept only in FASTP_PAIRED.

    if [[ "${params.cut_front.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--cut_front)
    fi

    if [[ "${params.cut_tail.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--cut_tail)
    fi

    if [[ "${params.trim_poly_g.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--trim_poly_g)
    fi

    if [[ "${params.trim_poly_x.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS+=(--trim_poly_x)
    fi

    echo "Running fastp interleaved trimming for sample: ${sample_id}" >&2
    echo "Input interleaved: ${interleaved}" >&2
    echo "Threads: ${task.cpus}" >&2
    echo "NOTE: --detect_adapter_for_pe and --correction are disabled for interleaved input." >&2

    fastp \\
        -i "${interleaved}" \\
        --interleaved_in \\
        --stdout \\
        --html "${safe_id}_fastp.html" \\
        --json "${safe_id}_fastp.json" \\
        --report_title "fastp report: ${sample_id}" \\
        --thread ${task.cpus} \\
        --cut_window_size ${params.cut_window_size} \\
        --cut_mean_quality ${params.cut_mean_quality} \\
        --qualified_quality_phred ${params.qualified_quality_phred} \\
        --unqualified_percent_limit ${params.unqualified_percent} \\
        --n_base_limit ${params.n_base_limit} \\
        --length_required ${params.length_required} \\
        "\${FASTP_ARGS[@]}" \\
        2> "${safe_id}_fastp.log" \\
        | gzip -${params.compression} > "${safe_id}_interleaved_trimmed.fastq.gz"

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" \\
        "${safe_id}" \\
        "interleaved" \\
        "" \\
        "" \\
        "${params.output_dir ?: params.outdir}/trimmed_reads/${safe_id}_interleaved_trimmed.fastq.gz" \\
        "" \\
        "${params.output_dir ?: params.outdir}/fastp_reports/${safe_id}_fastp.html" \\
        "${params.output_dir ?: params.outdir}/fastp_reports/${safe_id}_fastp.json" \\
        > "${safe_id}_trimmed_manifest_record.tsv"
    """
}

process TRIMMING_STATS_PAIRED {
    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/summary",
        mode: 'copy',
        pattern: "*_trimming_stats.tsv"

    input:
    tuple val(sample_id),
          val(safe_id),
          path(read1_trimmed),
          path(read2_trimmed),
          path(fastp_json)

    output:
    path "${safe_id}_trimming_stats.tsv", emit: stats_file

    script:
    """
    set -euo pipefail

    python3 - "${sample_id}" "${safe_id}" "paired" "${read1_trimmed}" "${read2_trimmed}" "" "" "${fastp_json}" "${safe_id}_trimming_stats.tsv" <<'PY'
import gzip
import json
import sys
from pathlib import Path

sample_id, safe_id, layout, read1_path, read2_path, merged_path, interleaved_path, json_path, out_tsv = sys.argv[1:]

def count_fastq_records(path):
    if path is None or str(path).strip() == "":
        return 0

    path = Path(path)

    if not path.exists():
        return 0

    if path.stat().st_size == 0:
        return 0

    opener = gzip.open if str(path).endswith(".gz") else open

    lines = 0
    with opener(path, "rt", errors="replace") as handle:
        for _ in handle:
            lines += 1

    if lines % 4 != 0:
        raise RuntimeError(f"FASTQ line count is not divisible by 4 for {path}: {lines} lines")

    return lines // 4

def get_nested(data, keys, default="NA"):
    value = data

    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]

    return value

with open(json_path) as handle:
    data = json.load(handle)

raw_reads = get_nested(data, ["summary", "before_filtering", "total_reads"])
after_filtering_reads = get_nested(data, ["summary", "after_filtering", "total_reads"])
adapter_trimmed_reads = get_nested(data, ["adapter_cutting", "adapter_trimmed_reads"])

try:
    raw_pairs = int(raw_reads) // 2
except Exception:
    raw_pairs = "NA"

try:
    reads_removed_by_filters = int(raw_reads) - int(after_filtering_reads)
except Exception:
    reads_removed_by_filters = "NA"

final_r1_reads = count_fastq_records(read1_path)
final_r2_reads = count_fastq_records(read2_path)

final_interleaved_reads = 0
final_merged_reads = 0

final_fastq_records_total = final_r1_reads + final_r2_reads

# Standard paired-end output. This is the number of complete paired records.
final_pair_or_fragment_records = min(final_r1_reads, final_r2_reads)

columns = [
    "sample_id",
    "safe_sample_id",
    "layout",
    "raw_reads",
    "raw_pairs",
    "reads_after_filtering_fastp",
    "reads_removed_by_filters",
    "reads_with_adapters_trimmed",
    "final_read1_reads",
    "final_read2_reads",
    "final_interleaved_reads",
    "final_merged_reads",
    "final_fastq_records_total",
    "final_pair_or_fragment_records",
    "fastp_json"
]

values = [
    sample_id,
    safe_id,
    layout,
    raw_reads,
    raw_pairs,
    after_filtering_reads,
    reads_removed_by_filters,
    adapter_trimmed_reads,
    final_r1_reads,
    final_r2_reads,
    final_interleaved_reads,
    final_merged_reads,
    final_fastq_records_total,
    final_pair_or_fragment_records,
    json_path
]

with open(out_tsv, "w") as out:
    print(*columns, sep="\\t", file=out)
    print(*values, sep="\\t", file=out)
PY
    """
}

process TRIMMING_STATS_INTERLEAVED {
    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/summary",
        mode: 'copy',
        pattern: "*_trimming_stats.tsv"

    input:
    tuple val(sample_id),
          val(safe_id),
          path(interleaved_trimmed),
          path(fastp_json)

    output:
    path "${safe_id}_trimming_stats.tsv", emit: stats_file

    script:
    """
    set -euo pipefail

    python3 - "${sample_id}" "${safe_id}" "interleaved" "" "" "" "${interleaved_trimmed}" "${fastp_json}" "${safe_id}_trimming_stats.tsv" <<'PY'
import gzip
import json
import sys
from pathlib import Path

sample_id, safe_id, layout, read1_path, read2_path, merged_path, interleaved_path, json_path, out_tsv = sys.argv[1:]

def count_fastq_records(path):
    if path is None or str(path).strip() == "":
        return 0

    path = Path(path)

    if not path.exists():
        return 0

    if path.stat().st_size == 0:
        return 0

    opener = gzip.open if str(path).endswith(".gz") else open

    lines = 0
    with opener(path, "rt", errors="replace") as handle:
        for _ in handle:
            lines += 1

    if lines % 4 != 0:
        raise RuntimeError(f"FASTQ line count is not divisible by 4 for {path}: {lines} lines")

    return lines // 4

def get_nested(data, keys, default="NA"):
    value = data

    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return default
        value = value[key]

    return value

with open(json_path) as handle:
    data = json.load(handle)

raw_reads = get_nested(data, ["summary", "before_filtering", "total_reads"])
after_filtering_reads = get_nested(data, ["summary", "after_filtering", "total_reads"])
adapter_trimmed_reads = get_nested(data, ["adapter_cutting", "adapter_trimmed_reads"])

try:
    raw_pairs = int(raw_reads) // 2
except Exception:
    raw_pairs = "NA"

try:
    reads_removed_by_filters = int(raw_reads) - int(after_filtering_reads)
except Exception:
    reads_removed_by_filters = "NA"

final_interleaved_reads = count_fastq_records(interleaved_path)

final_r1_reads = 0
final_r2_reads = 0
final_merged_reads = 0

final_fastq_records_total = final_interleaved_reads

try:
    final_pair_or_fragment_records = final_interleaved_reads // 2
except Exception:
    final_pair_or_fragment_records = "NA"

columns = [
    "sample_id",
    "safe_sample_id",
    "layout",
    "raw_reads",
    "raw_pairs",
    "reads_after_filtering_fastp",
    "reads_removed_by_filters",
    "reads_with_adapters_trimmed",
    "final_read1_reads",
    "final_read2_reads",
    "final_interleaved_reads",
    "final_merged_reads",
    "final_fastq_records_total",
    "final_pair_or_fragment_records",
    "fastp_json"
]

values = [
    sample_id,
    safe_id,
    layout,
    raw_reads,
    raw_pairs,
    after_filtering_reads,
    reads_removed_by_filters,
    adapter_trimmed_reads,
    final_r1_reads,
    final_r2_reads,
    final_interleaved_reads,
    final_merged_reads,
    final_fastq_records_total,
    final_pair_or_fragment_records,
    json_path
]

with open(out_tsv, "w") as out:
    print(*columns, sep="\\t", file=out)
    print(*values, sep="\\t", file=out)
PY
    """
}

process RUN_FASTQC_TRIMMED {
    tag { read_file.simpleName }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/fastqc_trimmed",
        mode: 'copy'

    cpus {
        params.threads != null
            ? params.threads as int
            : params.fastqc_threads as int
    }

    input:
    tuple val(sample_id), path(read_file), path(tools_status)

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    set -euo pipefail

    TOOLS_ENV="\$(grep '^ENV_DIR=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOLS_ENV" && "\$TOOLS_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOLS_ENV/bin:\$PATH"
    fi

    if ! command -v fastqc >/dev/null 2>&1; then
        echo "ERROR: FastQC is not available in PATH after setup." >&2
        echo "Tool setup status:" >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    echo "Running FastQC on trimmed read file: ${read_file}"
    echo "Sample: ${sample_id}"
    echo "FastQC threads: ${task.cpus}"

    fastqc \\
        -t ${task.cpus} \\
        "${read_file}" \\
        --outdir .
    """
}

process WRITE_TRIMMED_MANIFEST {
    tag "write_trimmed_manifest"

    publishDir "${params.output_dir ?: params.outdir}/summary",
        mode: 'copy',
        pattern: "trimmed_manifest.tsv"

    input:
    path manifest_records

    output:
    path "trimmed_manifest.tsv", emit: manifest

    script:
    def files = manifest_records.collect { record_file -> record_file.name }.join(' ')

    """
    set -euo pipefail

    printf 'sample_id\\tsafe_sample_id\\tlayout\\tread1\\tread2\\tinterleaved\\tmerged\\tfastp_html\\tfastp_json\\n' > trimmed_manifest.tsv

    if [[ -z "${files}" ]]; then
        echo "ERROR: No trimmed manifest record files were received." >&2
        exit 1
    fi

    for f in ${files}; do
        cat "\$f" >> trimmed_manifest.tsv
    done

    echo "Wrote trimmed manifest:"
    cat trimmed_manifest.tsv
    """
}

process WRITE_TRIMMING_STATS_SUMMARY {
    tag "write_trimming_stats_summary"

    publishDir "${params.output_dir ?: params.outdir}/summary",
        mode: 'copy',
        pattern: "trimming_stats_summary.tsv"

    input:
    path stats_files

    output:
    path "trimming_stats_summary.tsv", emit: summary

    script:
    def files = stats_files.collect { stats_file -> stats_file.name }.join(' ')

    """
    set -euo pipefail

    first=1
    : > trimming_stats_summary.tsv

    for f in ${files}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" >> trimming_stats_summary.tsv
            first=0
        else
            tail -n +2 "\$f" >> trimming_stats_summary.tsv
        fi
    done
    """
}