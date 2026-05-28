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

/*
 * If a pinned conda package installs but crashes, allow retrying with the latest
 * available unpinned package. This improves portability across macOS/Linux/HPC.
 */
params.allow_unpinned_tool_fallback = true

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
 * Conservative fastp defaults.
 *
 * fastp already performs basic adapter trimming and filtering by default.
 * These extra trimming options are disabled by default to avoid overly aggressive trimming.
 */
params.detect_adapter_for_pe   = false
params.enable_correction       = false

params.cut_front               = false
params.cut_tail                = false
params.cut_window_size         = 4
params.cut_mean_quality        = 20

params.qualified_quality_phred = 15
params.unqualified_percent     = 40
params.n_base_limit            = 5
params.length_required         = 15

params.trim_poly_g             = false
params.trim_poly_x             = false

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

    TRIMMING_STATS_PAIRED(FASTP_PAIRED.out.trimmed_reads)
    TRIMMING_STATS_INTERLEAVED(FASTP_INTERLEAVED.out.trimmed_reads)

    def all_trimming_stats_ch = TRIMMING_STATS_PAIRED.out.stats_file.mix(TRIMMING_STATS_INTERLEAVED.out.stats_file)

    WRITE_TRIMMING_STATS_SUMMARY(all_trimming_stats_ch.collect())

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
    def env_base = params.tool_env_dir ?: "${base_outdir}/conda_envs/module1_tools"

    """
    set -euo pipefail

    STATUS_FILE="module1_tools_status.env"

    TOOL_ENV_BASE="${env_base}"
    FASTP_ENV="\$TOOL_ENV_BASE/fastp_env"
    FASTQC_ENV="\$TOOL_ENV_BASE/fastqc_env"

    FASTP_ENV_VALUE=""
    FASTQC_ENV_VALUE=""

    echo "Module 1 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested fastp version: ${params.fastp_version}" >> "\$STATUS_FILE"
    echo "Requested FastQC version: ${params.fastqc_version}" >> "\$STATUS_FILE"
    echo "Allow unpinned fallback: ${params.allow_unpinned_tool_fallback}" >> "\$STATUS_FILE"
    echo "TOOL_ENV_BASE=\$TOOL_ENV_BASE" >> "\$STATUS_FILE"
    echo "FASTP_ENV_TARGET=\$FASTP_ENV" >> "\$STATUS_FILE"
    echo "FASTQC_ENV_TARGET=\$FASTQC_ENV" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    INSTALLER=""

    if command -v mamba >/dev/null 2>&1; then
        INSTALLER="mamba"
        echo "Using mamba: \$(command -v mamba)" >> "\$STATUS_FILE"
    elif command -v conda >/dev/null 2>&1; then
        INSTALLER="conda"
        echo "Using conda: \$(command -v conda)" >> "\$STATUS_FILE"
    else
        echo "No mamba/conda detected in PATH." >> "\$STATUS_FILE"
    fi

    ########################################
    # fastp setup
    ########################################

    echo "" >> "\$STATUS_FILE"
    echo "fastp setup" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    if [[ -d "\$FASTP_ENV" ]]; then
        echo "Existing fastp env detected: \$FASTP_ENV" >> "\$STATUS_FILE"

        if [[ -x "\$FASTP_ENV/bin/fastp" ]]; then
            echo "Testing existing fastp env..." >> "\$STATUS_FILE"

            set +e
            "\$FASTP_ENV/bin/fastp" --version >> "\$STATUS_FILE" 2>&1
            FASTP_TEST=\$?
            set -e

            if [[ "\$FASTP_TEST" -eq 0 ]]; then
                echo "Existing fastp env passed." >> "\$STATUS_FILE"
                FASTP_ENV_VALUE="\$FASTP_ENV"
            else
                echo "WARNING: Existing fastp env failed with code \$FASTP_TEST." >> "\$STATUS_FILE"
                echo "Removing broken fastp env." >> "\$STATUS_FILE"
                rm -rf "\$FASTP_ENV"
            fi
        else
            echo "WARNING: Existing fastp env is incomplete." >> "\$STATUS_FILE"
            echo "Removing incomplete fastp env." >> "\$STATUS_FILE"
            rm -rf "\$FASTP_ENV"
        fi
    fi

    if [[ -z "\$FASTP_ENV_VALUE" ]]; then
        if command -v fastp >/dev/null 2>&1; then
            echo "System/runtime fastp detected: \$(command -v fastp)" >> "\$STATUS_FILE"

            set +e
            fastp --version >> "\$STATUS_FILE" 2>&1
            FASTP_SYSTEM_TEST=\$?
            set -e

            if [[ "\$FASTP_SYSTEM_TEST" -eq 0 ]]; then
                echo "System/runtime fastp passed." >> "\$STATUS_FILE"
                FASTP_ENV_VALUE="SYSTEM"
            else
                echo "WARNING: System/runtime fastp exists but failed with code \$FASTP_SYSTEM_TEST." >> "\$STATUS_FILE"
            fi
        fi
    fi

    if [[ -z "\$FASTP_ENV_VALUE" ]]; then
        echo "fastp is not available as a working tool." >> "\$STATUS_FILE"

        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: Auto-install is disabled and fastp is missing." >> "\$STATUS_FILE"
            exit 1
        fi

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda was found in PATH." >> "\$STATUS_FILE"
            exit 1
        fi

        mkdir -p "\$(dirname "\$FASTP_ENV")"

        echo "Creating pinned fastp env:" >> "\$STATUS_FILE"
        echo "  \$FASTP_ENV" >> "\$STATUS_FILE"

        rm -rf "\$FASTP_ENV"

        "\$INSTALLER" create -y \\
            -p "\$FASTP_ENV" \\
            -c conda-forge \\
            -c bioconda \\
            "fastp=${params.fastp_version}" \\
            >> "\$STATUS_FILE" 2>&1

        if [[ ! -x "\$FASTP_ENV/bin/fastp" ]]; then
            echo "ERROR: fastp was not found after pinned installation." >> "\$STATUS_FILE"
            exit 1
        fi

        echo "Testing pinned fastp env..." >> "\$STATUS_FILE"

        set +e
        "\$FASTP_ENV/bin/fastp" --version >> "\$STATUS_FILE" 2>&1
        FASTP_PINNED_TEST=\$?
        set -e

        if [[ "\$FASTP_PINNED_TEST" -eq 0 ]]; then
            echo "Pinned fastp env passed." >> "\$STATUS_FILE"
            FASTP_ENV_VALUE="\$FASTP_ENV"
        else
            echo "WARNING: Pinned fastp env failed with code \$FASTP_PINNED_TEST." >> "\$STATUS_FILE"

            if [[ "${params.allow_unpinned_tool_fallback}" == "true" ]]; then
                echo "Attempting unpinned fastp fallback." >> "\$STATUS_FILE"
                rm -rf "\$FASTP_ENV"

                "\$INSTALLER" create -y \\
                    -p "\$FASTP_ENV" \\
                    -c conda-forge \\
                    -c bioconda \\
                    "fastp" \\
                    >> "\$STATUS_FILE" 2>&1

                if [[ ! -x "\$FASTP_ENV/bin/fastp" ]]; then
                    echo "ERROR: fastp was not found after unpinned fallback installation." >> "\$STATUS_FILE"
                    exit 1
                fi

                echo "Testing unpinned fastp fallback env..." >> "\$STATUS_FILE"
                "\$FASTP_ENV/bin/fastp" --version >> "\$STATUS_FILE" 2>&1

                echo "Unpinned fastp fallback passed." >> "\$STATUS_FILE"
                FASTP_ENV_VALUE="\$FASTP_ENV"
            else
                echo "ERROR: Pinned fastp failed and unpinned fallback is disabled." >> "\$STATUS_FILE"
                exit 1
            fi
        fi
    fi

    ########################################
    # FastQC setup
    ########################################

    echo "" >> "\$STATUS_FILE"
    echo "FastQC setup" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    if [[ -d "\$FASTQC_ENV" ]]; then
        echo "Existing FastQC env detected: \$FASTQC_ENV" >> "\$STATUS_FILE"

        if [[ -x "\$FASTQC_ENV/bin/fastqc" ]]; then
            echo "Testing existing FastQC env..." >> "\$STATUS_FILE"

            set +e
            "\$FASTQC_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1
            FASTQC_TEST=\$?
            set -e

            if [[ "\$FASTQC_TEST" -eq 0 ]]; then
                echo "Existing FastQC env passed." >> "\$STATUS_FILE"
                FASTQC_ENV_VALUE="\$FASTQC_ENV"
            else
                echo "WARNING: Existing FastQC env failed with code \$FASTQC_TEST." >> "\$STATUS_FILE"
                echo "Removing broken FastQC env." >> "\$STATUS_FILE"
                rm -rf "\$FASTQC_ENV"
            fi
        else
            echo "WARNING: Existing FastQC env is incomplete." >> "\$STATUS_FILE"
            echo "Removing incomplete FastQC env." >> "\$STATUS_FILE"
            rm -rf "\$FASTQC_ENV"
        fi
    fi

    if [[ -z "\$FASTQC_ENV_VALUE" ]]; then
        if command -v fastqc >/dev/null 2>&1; then
            echo "System/runtime FastQC detected: \$(command -v fastqc)" >> "\$STATUS_FILE"

            set +e
            fastqc --version >> "\$STATUS_FILE" 2>&1
            FASTQC_SYSTEM_TEST=\$?
            set -e

            if [[ "\$FASTQC_SYSTEM_TEST" -eq 0 ]]; then
                echo "System/runtime FastQC passed." >> "\$STATUS_FILE"
                FASTQC_ENV_VALUE="SYSTEM"
            else
                echo "WARNING: System/runtime FastQC exists but failed with code \$FASTQC_SYSTEM_TEST." >> "\$STATUS_FILE"
            fi
        fi
    fi

    if [[ -z "\$FASTQC_ENV_VALUE" ]]; then
        echo "FastQC is not available as a working tool." >> "\$STATUS_FILE"

        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: Auto-install is disabled and FastQC is missing." >> "\$STATUS_FILE"
            exit 1
        fi

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda was found in PATH." >> "\$STATUS_FILE"
            exit 1
        fi

        mkdir -p "\$(dirname "\$FASTQC_ENV")"

        echo "Creating pinned FastQC env:" >> "\$STATUS_FILE"
        echo "  \$FASTQC_ENV" >> "\$STATUS_FILE"

        rm -rf "\$FASTQC_ENV"

        "\$INSTALLER" create -y \\
            -p "\$FASTQC_ENV" \\
            -c conda-forge \\
            -c bioconda \\
            "perl" \\
            "fastqc=${params.fastqc_version}" \\
            >> "\$STATUS_FILE" 2>&1

        if [[ ! -x "\$FASTQC_ENV/bin/fastqc" ]]; then
            echo "ERROR: FastQC was not found after pinned installation." >> "\$STATUS_FILE"
            exit 1
        fi

        echo "Testing pinned FastQC env..." >> "\$STATUS_FILE"

        set +e
        "\$FASTQC_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1
        FASTQC_PINNED_TEST=\$?
        set -e

        if [[ "\$FASTQC_PINNED_TEST" -eq 0 ]]; then
            echo "Pinned FastQC env passed." >> "\$STATUS_FILE"
            FASTQC_ENV_VALUE="\$FASTQC_ENV"
        else
            echo "WARNING: Pinned FastQC env failed with code \$FASTQC_PINNED_TEST." >> "\$STATUS_FILE"

            if [[ "${params.allow_unpinned_tool_fallback}" == "true" ]]; then
                echo "Attempting unpinned FastQC fallback." >> "\$STATUS_FILE"
                rm -rf "\$FASTQC_ENV"

                "\$INSTALLER" create -y \\
                    -p "\$FASTQC_ENV" \\
                    -c conda-forge \\
                    -c bioconda \\
                    "perl" \\
                    "fastqc" \\
                    >> "\$STATUS_FILE" 2>&1

                if [[ ! -x "\$FASTQC_ENV/bin/fastqc" ]]; then
                    echo "ERROR: FastQC was not found after unpinned fallback installation." >> "\$STATUS_FILE"
                    exit 1
                fi

                echo "Testing unpinned FastQC fallback env..." >> "\$STATUS_FILE"
                "\$FASTQC_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1

                echo "Unpinned FastQC fallback passed." >> "\$STATUS_FILE"
                FASTQC_ENV_VALUE="\$FASTQC_ENV"
            else
                echo "ERROR: Pinned FastQC failed and unpinned fallback is disabled." >> "\$STATUS_FILE"
                exit 1
            fi
        fi
    fi

    echo "" >> "\$STATUS_FILE"
    echo "Final tool environment summary" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"
    echo "FASTP_ENV=\$FASTP_ENV_VALUE" >> "\$STATUS_FILE"
    echo "FASTQC_ENV=\$FASTQC_ENV_VALUE" >> "\$STATUS_FILE"
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

    FASTP_ENV="\$(grep '^FASTP_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$FASTP_ENV" && "\$FASTP_ENV" != "SYSTEM" ]]; then
        export PATH="\$FASTP_ENV/bin:\$PATH"
    fi

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH after setup." >&2
        echo "Tool setup status:" >&2
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

    echo "Running fastp paired-end trimming for sample: ${sample_id}"
    echo "Input R1: ${read1}"
    echo "Input R2: ${read2}"
    echo "Threads: ${task.cpus}"
    echo "Optional fastp args:\${FASTP_ARGS}"
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
        \${FASTP_ARGS} \\
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

    FASTP_ENV="\$(grep '^FASTP_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$FASTP_ENV" && "\$FASTP_ENV" != "SYSTEM" ]]; then
        export PATH="\$FASTP_ENV/bin:\$PATH"
    fi

    if ! command -v fastp >/dev/null 2>&1; then
        echo "ERROR: fastp is not available in PATH after setup." >&2
        echo "Tool setup status:" >&2
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

    echo "Running fastp interleaved trimming for sample: ${sample_id}" >&2
    echo "Input interleaved: ${interleaved}" >&2
    echo "Threads: ${task.cpus}" >&2
    echo "Optional fastp args:\${FASTP_ARGS}" >&2
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
        \${FASTP_ARGS} \\
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

    FASTQC_ENV="\$(grep '^FASTQC_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$FASTQC_ENV" && "\$FASTQC_ENV" != "SYSTEM" ]]; then
        export PATH="\$FASTQC_ENV/bin:\$PATH"
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