#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Module 1: Read trimming / quality control with fastp, followed by FastQC.
 */

/*
 * Parameters
 */

params.working_dir     = null
params.input_manifest  = null
params.output_dir      = null
params.fastp_version   = "0.23.4"
params.fastqc_version  = "0.12.1"
params.auto_install    = true
params.tool_env_dir    = null
params.fastp_threads   = 4
params.fastqc_threads  = 2
params.threads         = null
params.publish_trimmed_mode = "symlink"
params.compression          = 4
params.detect_adapter_for_pe   = true
params.enable_correction       = false
params.cut_front               = true
params.cut_tail                = true
params.cut_window_size         = 4
params.cut_mean_quality        = 30
params.qualified_quality_phred = 30
params.unqualified_percent     = 40
params.n_base_limit            = 5
params.length_required         = 75
params.trim_poly_g             = false
params.trim_poly_x             = false
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.module0_outdir = "${params.results_dir}/module_0_readprocess"
params.outdir = "${params.results_dir}/module_1_readtrimming"

workflow {

    def manifest_file = params.input_manifest ?: "${params.module0_outdir}/naming/read_manifest.tsv"

    log.info "Module 1 working_dir parameter: ${params.working_dir ?: 'not supplied'}"
    log.info "Module 1 results directory: ${params.results_dir}"
    log.info "Using Module 0 manifest: ${manifest_file}"
    log.info "Writing Module 1 outputs to: ${params.outdir}"

    def manifest_ch = channel.fromPath(
        manifest_file,
        type: 'file',
        checkIfExists: true
    )

    def rows_ch = manifest_ch.splitCsv(header: true, sep: '\t')

    def paired_reads_ch = rows_ch
        .filter { row ->
            row.layout.toString() == 'paired'
        }
        .map { row ->
            def sample_id = row.sample_id.toString()
            def safe_id   = sample_id.replaceAll('[^A-Za-z0-9._-]+', '_')

            tuple(
                sample_id,
                safe_id,
                file(row.read1),
                file(row.read2)
            )
        }

    def interleaved_reads_ch = rows_ch
        .filter { row ->
            row.layout.toString() == 'interleaved'
        }
        .map { row ->
            def sample_id = row.sample_id.toString()
            def safe_id   = sample_id.replaceAll('[^A-Za-z0-9._-]+', '_')

            tuple(
                sample_id,
                safe_id,
                file(row.interleaved)
            )
        }

    SETUP_MODULE1_TOOLS()

    FASTP_PAIRED(
        paired_reads_ch.combine(SETUP_MODULE1_TOOLS.out.status)
    )

    FASTP_INTERLEAVED(
        interleaved_reads_ch.combine(SETUP_MODULE1_TOOLS.out.status)
    )

    TRIMMING_STATS_PAIRED(
        FASTP_PAIRED.out.trimmed_reads
    )

    TRIMMING_STATS_INTERLEAVED(
        FASTP_INTERLEAVED.out.trimmed_reads
    )

    def all_stats_ch = TRIMMING_STATS_PAIRED.out.stats_file.mix(
        TRIMMING_STATS_INTERLEAVED.out.stats_file
    )

    WRITE_TRIMMING_STATS_SUMMARY(
        all_stats_ch.collect()
    )

    def paired_fastqc_reads_ch = FASTP_PAIRED.out.trimmed_reads
        .flatMap { sample_id, _safe_id, read1_trimmed, read2_trimmed, _fastp_json ->
            [
                tuple(sample_id, read1_trimmed),
                tuple(sample_id, read2_trimmed)
            ]
        }

    def interleaved_fastqc_reads_ch = FASTP_INTERLEAVED.out.trimmed_reads
        .map { sample_id, _safe_id, interleaved_trimmed, _fastp_json ->
            tuple(sample_id, interleaved_trimmed)
        }

    def fastqc_reads_ch = paired_fastqc_reads_ch.mix(interleaved_fastqc_reads_ch)

    RUN_FASTQC_TRIMMED(
        fastqc_reads_ch.combine(SETUP_MODULE1_TOOLS.out.status)
    )

    def manifest_records_ch = FASTP_PAIRED.out.manifest_record.mix(
        FASTP_INTERLEAVED.out.manifest_record
    )

    WRITE_TRIMMED_MANIFEST(
        manifest_records_ch.collect()
    )
}


process SETUP_MODULE1_TOOLS {

    tag "setup_fastp_fastqc"

    publishDir "${params.outdir}/setup",
        mode: 'copy',
        pattern: "module1_tools_status.env"

    output:
    path "module1_tools_status.env", emit: status

    script:
    def env_dir = params.tool_env_dir ?: "${params.outdir}/conda_envs/module1_tools"

    """
    set -euo pipefail

    STATUS_FILE="module1_tools_status.env"
    TOOL_ENV="${env_dir}"

    echo "Module 1 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested fastp version: ${params.fastp_version}" >> "\$STATUS_FILE"
    echo "Requested FastQC version: ${params.fastqc_version}" >> "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    if [[ -d "\$TOOL_ENV" ]]; then
        echo "Existing Module 1 environment detected: \$TOOL_ENV" >> "\$STATUS_FILE"

        if [[ -x "\$TOOL_ENV/bin/fastp" && -x "\$TOOL_ENV/bin/fastqc" ]]; then
            set +e
            "\$TOOL_ENV/bin/fastp" --version >> "\$STATUS_FILE" 2>&1
            FASTP_TEST=\$?
            "\$TOOL_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1
            FASTQC_TEST=\$?
            set -e

            if [[ "\$FASTP_TEST" -eq 0 && "\$FASTQC_TEST" -eq 0 ]]; then
                echo "Existing Module 1 environment passed checks." >> "\$STATUS_FILE"
                echo "FASTP_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
                echo "FASTQC_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
                echo "Module 1 tool setup finished: \$(date)" >> "\$STATUS_FILE"
                exit 0
            fi
        fi

        echo "Existing Module 1 environment is incomplete or broken. Removing." >> "\$STATUS_FILE"
        rm -rf "\$TOOL_ENV"
    fi

    SYSTEM_FASTP="false"
    SYSTEM_FASTQC="false"

    if command -v fastp >/dev/null 2>&1; then
        set +e
        fastp --version >> "\$STATUS_FILE" 2>&1
        FASTP_SYSTEM_TEST=\$?
        set -e

        if [[ "\$FASTP_SYSTEM_TEST" -eq 0 ]]; then
            SYSTEM_FASTP="true"
        fi
    fi

    if command -v fastqc >/dev/null 2>&1; then
        set +e
        fastqc --version >> "\$STATUS_FILE" 2>&1
        FASTQC_SYSTEM_TEST=\$?
        set -e

        if [[ "\$FASTQC_SYSTEM_TEST" -eq 0 ]]; then
            SYSTEM_FASTQC="true"
        fi
    fi

    if [[ "\$SYSTEM_FASTP" == "true" && "\$SYSTEM_FASTQC" == "true" ]]; then
        echo "System fastp and FastQC are both available." >> "\$STATUS_FILE"
        echo "FASTP_ENV=SYSTEM" >> "\$STATUS_FILE"
        echo "FASTQC_ENV=SYSTEM" >> "\$STATUS_FILE"
        echo "Module 1 tool setup finished: \$(date)" >> "\$STATUS_FILE"
        exit 0
    fi

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: fastp and/or FastQC missing, and auto_install is false." >> "\$STATUS_FILE"
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
        echo "ERROR: Neither mamba nor conda found in PATH." >> "\$STATUS_FILE"
        exit 1
    fi

    mkdir -p "\$(dirname "\$TOOL_ENV")"

    echo "Creating Module 1 environment:" >> "\$STATUS_FILE"
    echo "  \$TOOL_ENV" >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \\
        -p "\$TOOL_ENV" \\
        -c conda-forge \\
        -c bioconda \\
        "perl" \\
        "fastp=${params.fastp_version}" \\
        "fastqc=${params.fastqc_version}" \\
        >> "\$STATUS_FILE" 2>&1

    if [[ ! -x "\$TOOL_ENV/bin/fastp" ]]; then
        echo "ERROR: fastp not found after installation." >> "\$STATUS_FILE"
        exit 1
    fi

    if [[ ! -x "\$TOOL_ENV/bin/fastqc" ]]; then
        echo "ERROR: FastQC not found after installation." >> "\$STATUS_FILE"
        exit 1
    fi

    "\$TOOL_ENV/bin/fastp" --version >> "\$STATUS_FILE" 2>&1
    "\$TOOL_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1

    echo "FASTP_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "FASTQC_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Module 1 tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}

process FASTP_PAIRED {

    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.outdir}/trimmed_reads",
        mode: params.publish_trimmed_mode,
        pattern: "*.fastq.gz"

    publishDir "${params.outdir}/fastp_reports",
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

    FASTP_ENV="\$(grep '^FASTP_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$FASTP_ENV" && "\$FASTP_ENV" != "SYSTEM" ]]; then
        export PATH="\$FASTP_ENV/bin:\$PATH"
    fi

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH after setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    FASTP_ARGS=""

    if [[ "${params.detect_adapter_for_pe.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --detect_adapter_for_pe"
    fi

    if [[ "${params.enable_correction.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --correction"
    fi

    if [[ "${params.cut_front.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --cut_front"
    fi

    if [[ "${params.cut_tail.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --cut_tail"
    fi

    if [[ "${params.trim_poly_g.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --trim_poly_g"
    fi

    if [[ "${params.trim_poly_x.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --trim_poly_x"
    fi

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
        \${FASTP_ARGS} \\
        > "${safe_id}_fastp.log" 2>&1

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" \\
        "${safe_id}" \\
        "paired" \\
        "${params.outdir}/trimmed_reads/${safe_id}_R1_trimmed.fastq.gz" \\
        "${params.outdir}/trimmed_reads/${safe_id}_R2_trimmed.fastq.gz" \\
        "" \\
        "" \\
        "${params.outdir}/fastp_reports/${safe_id}_fastp.html" \\
        "${params.outdir}/fastp_reports/${safe_id}_fastp.json" \\
        > "${safe_id}_trimmed_manifest_record.tsv"
    """
}

process FASTP_INTERLEAVED {

    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.outdir}/trimmed_reads",
        mode: params.publish_trimmed_mode,
        pattern: "*.fastq.gz"

    publishDir "${params.outdir}/fastp_reports",
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

    FASTP_ENV="\$(grep '^FASTP_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$FASTP_ENV" && "\$FASTP_ENV" != "SYSTEM" ]]; then
        export PATH="\$FASTP_ENV/bin:\$PATH"
    fi

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH after setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    FASTP_ARGS=""

    if [[ "${params.cut_front.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --cut_front"
    fi

    if [[ "${params.cut_tail.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --cut_tail"
    fi

    if [[ "${params.trim_poly_g.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --trim_poly_g"
    fi

    if [[ "${params.trim_poly_x.toString().toLowerCase()}" == "true" ]]; then
        FASTP_ARGS="\${FASTP_ARGS} --trim_poly_x"
    fi

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
        \${FASTP_ARGS} \\
        2> "${safe_id}_fastp.log" \\
        | gzip -${params.compression} > "${safe_id}_interleaved_trimmed.fastq.gz"

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" \\
        "${safe_id}" \\
        "interleaved" \\
        "" \\
        "" \\
        "${params.outdir}/trimmed_reads/${safe_id}_interleaved_trimmed.fastq.gz" \\
        "" \\
        "${params.outdir}/fastp_reports/${safe_id}_fastp.html" \\
        "${params.outdir}/fastp_reports/${safe_id}_fastp.json" \\
        > "${safe_id}_trimmed_manifest_record.tsv"
    """
}

process TRIMMING_STATS_PAIRED {

    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.outdir}/summary",
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

    python3 - "${sample_id}" "${safe_id}" "paired" "${read1_trimmed}" "${read2_trimmed}" "${fastp_json}" "${safe_id}_trimming_stats.tsv" <<'PY'
import gzip
import json
import sys
from pathlib import Path

sample_id, safe_id, layout, r1, r2, json_path, out_tsv = sys.argv[1:]

def count_fastq_records(path):
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    lines = 0
    with opener(path, "rt", errors="replace") as handle:
        for _ in handle:
            lines += 1
    if lines % 4 != 0:
        raise RuntimeError(f"FASTQ line count not divisible by 4: {path}")
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
after_reads = get_nested(data, ["summary", "after_filtering", "total_reads"])
adapter_reads = get_nested(data, ["adapter_cutting", "adapter_trimmed_reads"])

try:
    raw_pairs = int(raw_reads) // 2
except Exception:
    raw_pairs = "NA"

try:
    removed = int(raw_reads) - int(after_reads)
except Exception:
    removed = "NA"

r1_count = count_fastq_records(r1)
r2_count = count_fastq_records(r2)

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
    after_reads,
    removed,
    adapter_reads,
    r1_count,
    r2_count,
    0,
    r1_count + r2_count,
    min(r1_count, r2_count),
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

    publishDir "${params.outdir}/summary",
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

    python3 - "${sample_id}" "${safe_id}" "interleaved" "${interleaved_trimmed}" "${fastp_json}" "${safe_id}_trimming_stats.tsv" <<'PY'
import gzip
import json
import sys
from pathlib import Path

sample_id, safe_id, layout, interleaved, json_path, out_tsv = sys.argv[1:]

def count_fastq_records(path):
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    lines = 0
    with opener(path, "rt", errors="replace") as handle:
        for _ in handle:
            lines += 1
    if lines % 4 != 0:
        raise RuntimeError(f"FASTQ line count not divisible by 4: {path}")
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
after_reads = get_nested(data, ["summary", "after_filtering", "total_reads"])
adapter_reads = get_nested(data, ["adapter_cutting", "adapter_trimmed_reads"])

try:
    raw_pairs = int(raw_reads) // 2
except Exception:
    raw_pairs = "NA"

try:
    removed = int(raw_reads) - int(after_reads)
except Exception:
    removed = "NA"

interleaved_count = count_fastq_records(interleaved)

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
    after_reads,
    removed,
    adapter_reads,
    0,
    0,
    interleaved_count,
    interleaved_count,
    interleaved_count // 2,
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

    publishDir "${params.outdir}/fastqc_reports",
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

    FASTQC_ENV="\$(grep '^FASTQC_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$FASTQC_ENV" && "\$FASTQC_ENV" != "SYSTEM" ]]; then
        export PATH="\$FASTQC_ENV/bin:\$PATH"
    fi

    if ! command -v fastqc >/dev/null 2>&1; then
        echo "ERROR: FastQC is not available in PATH after setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    fastqc \\
        -t ${task.cpus} \\
        "${read_file}" \\
        --outdir .
    """
}

process WRITE_TRIMMED_MANIFEST {

    tag "write_trimmed_manifest"

    publishDir "${params.outdir}/summary",
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
    """
}

process WRITE_TRIMMING_STATS_SUMMARY {

    tag "write_trimming_stats_summary"

    publishDir "${params.outdir}/summary",
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

    if [[ -z "${files}" ]]; then
        echo "ERROR: No trimming stats files were received." >&2
        exit 1
    fi

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