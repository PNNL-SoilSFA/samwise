#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Module 1: Read trimming / quality control with fastp, followed by FastQC.
 *
 * Input:
 *   The read_manifest.tsv produced by module 0.
 *
 * Example:
 *   nextflow run module_1_readtrimming.nf \
 *     --input_manifest ./results/module_0_readprocess/naming/read_manifest.tsv \
 *     -with-conda
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

/*
 * Paired reads are run with fastp merge mode.
 * The merged output is produced, but FastQC is not run on it by default because
 * merged output can sometimes be empty or very small depending on overlap.
 */
params.fastqc_merged           = false

workflow {
    def manifest_file = params.input_manifest ?: "${params.module0_outdir}/naming/read_manifest.tsv"
    def outbase       = params.output_dir ?: params.outdir

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

    FASTP_PAIRED(paired_reads_ch)
    FASTP_INTERLEAVED(interleaved_reads_ch)

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
     *
     * Optional:
     *   If --fastqc_merged true, also run FastQC on paired merged output.
     */
    def paired_fastqc_reads_ch = FASTP_PAIRED.out.trimmed_reads
        .flatMap { sample_id, _safe_id, read1_trimmed, read2_trimmed, merged_trimmed, _fastp_json ->
            def reads = [
                tuple(sample_id, read1_trimmed),
                tuple(sample_id, read2_trimmed)
            ]

            if( params.fastqc_merged.toString().toBoolean() ) {
                reads << tuple(sample_id, merged_trimmed)
            }

            return reads
        }

    def interleaved_fastqc_reads_ch = FASTP_INTERLEAVED.out.trimmed_reads
        .map { sample_id, _safe_id, interleaved_trimmed, _fastp_json ->
            tuple(sample_id, interleaved_trimmed)
        }

    def trimmed_fastqc_input_ch = paired_fastqc_reads_ch.mix(interleaved_fastqc_reads_ch)

    RUN_FASTQC_TRIMMED(trimmed_fastqc_input_ch)

    /*
     * Write a public manifest pointing to published module 1 output paths.
     */
    def paired_manifest_records_ch = FASTP_PAIRED.out.manifest_record
        .map { sample_id, safe_id, layout, read1_name, read2_name, _interleaved_name, merged_name, html_name, json_name ->
            tuple(
                sample_id,
                safe_id,
                layout,
                "${outbase}/trimmed_reads/${read1_name}",
                "${outbase}/trimmed_reads/${read2_name}",
                "",
                "${outbase}/trimmed_reads/${merged_name}",
                "${outbase}/fastp_reports/${html_name}",
                "${outbase}/fastp_reports/${json_name}"
            )
        }

    def interleaved_manifest_records_ch = FASTP_INTERLEAVED.out.manifest_record
        .map { sample_id, safe_id, layout, _read1_name, _read2_name, interleaved_name, _merged_name, html_name, json_name ->
            tuple(
                sample_id,
                safe_id,
                layout,
                "",
                "",
                "${outbase}/trimmed_reads/${interleaved_name}",
                "",
                "${outbase}/fastp_reports/${html_name}",
                "${outbase}/fastp_reports/${json_name}"
            )
        }

    def all_manifest_records_ch = paired_manifest_records_ch.mix(interleaved_manifest_records_ch)

    WRITE_TRIMMED_MANIFEST(all_manifest_records_ch.collect())
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

    conda "bioconda::fastp=${params.fastp_version}"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.fastp_threads as int
    }

    input:
    tuple val(sample_id), val(safe_id), path(read1), path(read2)

    output:
    tuple val(sample_id),
          val(safe_id),
          path("${safe_id}_R1_trimmed.fastq.gz"),
          path("${safe_id}_R2_trimmed.fastq.gz"),
          path("${safe_id}_merged_trimmed.fastq.gz"),
          path("${safe_id}_fastp.json"),
          emit: trimmed_reads

    tuple val(sample_id),
          val(safe_id),
          val("paired"),
          val("${safe_id}_R1_trimmed.fastq.gz"),
          val("${safe_id}_R2_trimmed.fastq.gz"),
          val(""),
          val("${safe_id}_merged_trimmed.fastq.gz"),
          val("${safe_id}_fastp.html"),
          val("${safe_id}_fastp.json"),
          emit: manifest_record

    path "${safe_id}_fastp.html", emit: html
    path "${safe_id}_fastp.log",  emit: log

    script:
    """
    set -euo pipefail

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH." >&2
        echo "Use '-with-conda', configure conda/mamba support, or install fastp in the runtime environment." >&2
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

    fastp \\
        -i "${read1}" \\
        -I "${read2}" \\
        -o "${safe_id}_R1_trimmed.fastq.gz" \\
        -O "${safe_id}_R2_trimmed.fastq.gz" \\
        -m \\
        --merged_out "${safe_id}_merged_trimmed.fastq.gz" \\
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

    # Be defensive: if fastp completes successfully but produces no merged reads,
    # ensure the declared merged output exists as a valid empty gzip FASTQ.
    if [[ ! -f "${safe_id}_merged_trimmed.fastq.gz" ]]; then
        gzip -c /dev/null > "${safe_id}_merged_trimmed.fastq.gz"
    fi
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

    conda "bioconda::fastp=${params.fastp_version}"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.fastp_threads as int
    }

    input:
    tuple val(sample_id), val(safe_id), path(interleaved)

    output:
    tuple val(sample_id),
          val(safe_id),
          path("${safe_id}_interleaved_trimmed.fastq.gz"),
          path("${safe_id}_fastp.json"),
          emit: trimmed_reads

    tuple val(sample_id),
          val(safe_id),
          val("interleaved"),
          val(""),
          val(""),
          val("${safe_id}_interleaved_trimmed.fastq.gz"),
          val(""),
          val("${safe_id}_fastp.html"),
          val("${safe_id}_fastp.json"),
          emit: manifest_record

    path "${safe_id}_fastp.html", emit: html
    path "${safe_id}_fastp.log",  emit: log

    script:
    """
    set -euo pipefail

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH." >&2
        echo "Use '-with-conda', configure conda/mamba support, or install fastp in the runtime environment." >&2
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

    echo "Running fastp interleaved trimming for sample: ${sample_id}" >&2
    echo "Input interleaved: ${interleaved}" >&2
    echo "Threads: ${task.cpus}" >&2

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
          path(merged_trimmed),
          path(fastp_json)

    output:
    path "${safe_id}_trimming_stats.tsv", emit: stats_file

    script:
    """
    set -euo pipefail

    python3 - "${sample_id}" "${safe_id}" "paired" "${read1_trimmed}" "${read2_trimmed}" "${merged_trimmed}" "" "${fastp_json}" "${safe_id}_trimming_stats.tsv" <<'PY'
import gzip
import json
import sys
from pathlib import Path

sample_id, safe_id, layout, read1_path, read2_path, merged_path, interleaved_path, json_path, out_tsv = sys.argv[1:]

def count_fastq_records(path):
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
passed_filter_reads = get_nested(data, ["filtering_result", "passed_filter_reads"], after_filtering_reads)
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
final_merged_reads = count_fastq_records(merged_path)
final_interleaved_reads = 0

final_fastq_records_total = final_r1_reads + final_r2_reads + final_merged_reads

# For paired output with merge mode:
#   each unmerged pair contributes one R1 record and one R2 record
#   each merged pair contributes one merged record
final_pair_or_fragment_records = min(final_r1_reads, final_r2_reads) + final_merged_reads

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
passed_filter_reads = get_nested(data, ["filtering_result", "passed_filter_reads"], after_filtering_reads)
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

    conda "bioconda::fastqc=${params.fastqc_version}"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.fastqc_threads as int
    }

    input:
    tuple val(sample_id), path(read_file)

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    set -euo pipefail

    if ! command -v fastqc >/dev/null 2>&1; then
        echo "ERROR: FastQC is not available in PATH." >&2
        echo "Use '-with-conda', configure conda/mamba support, or install FastQC in the runtime environment." >&2
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
    val records

    output:
    path "trimmed_manifest.tsv", emit: manifest

    script:
    def rows = records
        .collect { record -> record.join('\t') }
        .join('\n')

    """
    set -euo pipefail

    cat > trimmed_manifest.tsv <<'EOF'
sample_id\tsafe_sample_id\tlayout\tread1\tread2\tinterleaved\tmerged\tfastp_html\tfastp_json
${rows}
EOF
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