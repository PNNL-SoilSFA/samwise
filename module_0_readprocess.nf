#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Module 0: Read naming validation, optional FASTQ structure validation,
 * and raw-read FastQC.
 *
 * This module auto-installs FastQC if missing, using mamba or conda.
 */

/*
 * Default parameters.
 * Command-line values override these.
 */
params.input_dir       = null
params.outdir          = "./results/module_0_readprocess"
params.output_dir      = null
params.file_pattern    = "*.{fastq.gz,fq.gz,fastq,fq}"
params.fastqc_threads  = 2
params.threads         = null
params.skip_validate   = false

params.fastqc_version  = "0.12.1"
params.auto_install    = true
params.tool_env_dir    = null

workflow {
    if( !params.input_dir ) {
        error """
        Missing required parameter: --input_dir

        Example:
          nextflow run module_0_readprocess.nf --input_dir ./reads
        """.stripIndent()
    }

    def all_files_ch = channel.fromPath(
        "${params.input_dir}/${params.file_pattern}",
        type: 'file',
        checkIfExists: true
    )
    .map { path -> path.toAbsolutePath() }

    SETUP_MODULE0_TOOLS()

    CHECK_READ_NAMING(all_files_ch.collect())

    def valid_named_reads_ch = CHECK_READ_NAMING.out.manifest
        .splitCsv(header: true, sep: '\t')
        .flatMap { row ->
            if( row.layout == 'paired' ) {
                return [ file(row.read1), file(row.read2) ]
            }
            else if( row.layout == 'interleaved' ) {
                return [ file(row.interleaved) ]
            }
            else {
                return []
            }
        }

    def reads_ready_for_fastqc_ch

    if( params.skip_validate.toString().toBoolean() ) {
        SKIP_VALIDATE_READS(valid_named_reads_ch)
        reads_ready_for_fastqc_ch = SKIP_VALIDATE_READS.out
    }
    else {
        VALIDATE_READS(valid_named_reads_ch)
        reads_ready_for_fastqc_ch = VALIDATE_READS.out
    }

    def reads_ready_with_tools_ch = reads_ready_for_fastqc_ch.combine(SETUP_MODULE0_TOOLS.out.status)

    RUN_FASTQC(reads_ready_with_tools_ch)
}

process SETUP_MODULE0_TOOLS {
    tag "setup_fastqc"

    publishDir "${params.output_dir ?: params.outdir}/setup",
        mode: 'copy',
        pattern: "module0_tools_status.env"

    output:
    path "module0_tools_status.env", emit: status

    script:
    def base_outdir = params.output_dir ?: params.outdir
    def env_dir = params.tool_env_dir ?: "${base_outdir}/conda_envs/module0_tools"

    """
    set -euo pipefail

    STATUS_FILE="module0_tools_status.env"
    TOOL_ENV="${env_dir}"

    echo "Module 0 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested FastQC version: ${params.fastqc_version}" >> "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    if [[ -d "\$TOOL_ENV" ]]; then
        echo "Existing module-local environment detected: \$TOOL_ENV" >> "\$STATUS_FILE"

        if [[ -x "\$TOOL_ENV/bin/fastqc" ]]; then
            echo "Testing existing FastQC environment..." >> "\$STATUS_FILE"

            set +e
            "\$TOOL_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1
            FASTQC_TEST=\$?
            set -e

            if [[ "\$FASTQC_TEST" -eq 0 ]]; then
                echo "Existing module-local environment passed FastQC check." >> "\$STATUS_FILE"
                echo "ENV_DIR=\$TOOL_ENV" >> "\$STATUS_FILE"
                echo "Module 0 tool setup finished: \$(date)" >> "\$STATUS_FILE"
                exit 0
            else
                echo "WARNING: Existing module-local environment failed FastQC check." >> "\$STATUS_FILE"
                echo "FASTQC_TEST=\$FASTQC_TEST" >> "\$STATUS_FILE"
                echo "Removing broken environment and recreating it." >> "\$STATUS_FILE"
                rm -rf "\$TOOL_ENV"
            fi
        else
            echo "WARNING: Existing module-local environment is incomplete." >> "\$STATUS_FILE"
            echo "Removing incomplete environment and recreating it." >> "\$STATUS_FILE"
            rm -rf "\$TOOL_ENV"
        fi
    fi

    if command -v fastqc >/dev/null 2>&1; then
        echo "System/runtime FastQC detected." >> "\$STATUS_FILE"
        echo "fastqc path: \$(command -v fastqc)" >> "\$STATUS_FILE"

        set +e
        fastqc --version >> "\$STATUS_FILE" 2>&1
        FASTQC_TEST=\$?
        set -e

        if [[ "\$FASTQC_TEST" -eq 0 ]]; then
            echo "System/runtime FastQC passed check." >> "\$STATUS_FILE"
            echo "ENV_DIR=SYSTEM" >> "\$STATUS_FILE"
            echo "Module 0 tool setup finished: \$(date)" >> "\$STATUS_FILE"
            exit 0
        else
            echo "WARNING: System/runtime FastQC was detected but failed check." >> "\$STATUS_FILE"
            echo "Will attempt module-local installation." >> "\$STATUS_FILE"
        fi
    fi

    echo "FastQC is not available as a working tool." >> "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: Auto-install is disabled." >> "\$STATUS_FILE"
        echo "Install FastQC manually or rerun with --auto_install true." >> "\$STATUS_FILE"
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
        echo "Please install mamba/conda or install FastQC manually." >> "\$STATUS_FILE"
        exit 1
    fi

    mkdir -p "\$(dirname "\$TOOL_ENV")"

    echo "Creating module-local environment:" >> "\$STATUS_FILE"
    echo "  \$TOOL_ENV" >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \\
        -p "\$TOOL_ENV" \\
        -c conda-forge \\
        -c bioconda \\
        "perl" \\
        "fastqc=${params.fastqc_version}" \\
        >> "\$STATUS_FILE" 2>&1

    if [[ ! -x "\$TOOL_ENV/bin/fastqc" ]]; then
        echo "ERROR: FastQC was not found after installation." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "Testing newly created FastQC environment..." >> "\$STATUS_FILE"
    "\$TOOL_ENV/bin/fastqc" --version >> "\$STATUS_FILE" 2>&1

    echo "ENV_DIR=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "FastQC installed at: \$TOOL_ENV/bin/fastqc" >> "\$STATUS_FILE"
    echo "Module 0 tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}

process CHECK_READ_NAMING {
    tag "check_read_naming"

    publishDir "${params.output_dir ?: params.outdir}/naming",
        mode: 'copy',
        pattern: "*.{txt,tsv}"

    input:
    val read_files

    output:
    path "read_naming_report.txt", emit: report
    path "read_manifest.tsv", emit: manifest

    script:
    """
    set -euo pipefail

    cat > input_files.list <<'EOF'
${read_files.join('\n')}
EOF

    python3 - <<'PY'
import re
import sys
from pathlib import Path

report_path = Path("read_naming_report.txt")
manifest_path = Path("read_manifest.tsv")

r1_files = {}
r2_files = {}
interleaved_files = {}

r1_style = {}
r2_style = {}

errors = 0
warnings = 0

read1_patterns = [
    (re.compile(r'^(.+)_R1(?:_001)?\\.(fastq|fq)(\\.gz)?\\Z'), "R"),
    (re.compile(r'^(.+)_1\\.(fastq|fq)(\\.gz)?\\Z'), "numeric"),
]

read2_patterns = [
    (re.compile(r'^(.+)_R2(?:_001)?\\.(fastq|fq)(\\.gz)?\\Z'), "R"),
    (re.compile(r'^(.+)_2\\.(fastq|fq)(\\.gz)?\\Z'), "numeric"),
]

interleaved_re = re.compile(r'^(.+)_interleaved\\.(fastq|fq)(\\.gz)?\\Z')
fastq_like_re = re.compile(r'.*\\.(fastq|fq)(\\.gz)?\\Z')

def classify_read(name):
    for regex, style in read1_patterns:
        m = regex.match(name)
        if m:
            return "R1", m.group(1), style

    for regex, style in read2_patterns:
        m = regex.match(name)
        if m:
            return "R2", m.group(1), style

    m = interleaved_re.match(name)
    if m:
        return "interleaved", m.group(1), "interleaved"

    return None, None, None

with report_path.open("w") as report, manifest_path.open("w") as manifest:
    def log(message=""):
        print(message, file=report)

    print("sample_id", "layout", "read1", "read2", "interleaved", sep="\\t", file=manifest)

    log("Read naming and pairing report")
    log("Started")
    log("----------------------------------------")
    log("")

    with open("input_files.list") as handle:
        input_files = [line.strip() for line in handle if line.strip()]

    if not input_files:
        log("ERROR: No files were found in the input directory.")
        errors += 1

    log("Files detected:")
    for f in input_files:
        log("  {}".format(f))
    log("")

    for file_path in input_files:
        path = Path(file_path)
        name = path.name

        read_type, sample, style = classify_read(name)

        if read_type == "interleaved":
            if sample in interleaved_files:
                log("ERROR: Duplicate interleaved file for sample '{}': {}".format(sample, name))
                errors += 1
            interleaved_files[sample] = file_path
            continue

        if read_type == "R1":
            if sample in r1_files:
                log("ERROR: Duplicate R1 file for sample '{}': {}".format(sample, name))
                errors += 1
            r1_files[sample] = file_path
            r1_style[sample] = style
            continue

        if read_type == "R2":
            if sample in r2_files:
                log("ERROR: Duplicate R2 file for sample '{}': {}".format(sample, name))
                errors += 1
            r2_files[sample] = file_path
            r2_style[sample] = style
            continue

        if fastq_like_re.match(name):
            log("ERROR: FASTQ file has unsupported naming convention: {}".format(name))
            log("       Expected examples:")
            log("         sample_R1.fastq or sample_R1.fastq.gz")
            log("         sample_R2.fastq or sample_R2.fastq.gz")
            log("         sample_1.fastq  or sample_1.fastq.gz")
            log("         sample_2.fastq  or sample_2.fastq.gz")
            log("         sample_S1_L001_R1_001.fastq.gz")
            log("         sample_S1_L001_R2_001.fastq.gz")
            log("         sample_interleaved.fastq or sample_interleaved.fastq.gz")
            errors += 1
        else:
            log("ERROR: Non-FASTQ file detected: {}".format(name))
            errors += 1

    log("")
    log("Pair/interleaved checks:")
    log("----------------------------------------")

    emitted_files = []

    all_paired_samples = sorted(set(r1_files) | set(r2_files))

    for sample in all_paired_samples:
        has_r1 = sample in r1_files
        has_r2 = sample in r2_files

        if not has_r1:
            log("ERROR: Sample '{}' has R2 but no R1.".format(sample))
            errors += 1
            continue

        if not has_r2:
            log("ERROR: Sample '{}' has R1 but no R2.".format(sample))
            errors += 1
            continue

        if sample in interleaved_files:
            log("ERROR: Sample '{}' has both paired-end and interleaved files.".format(sample))
            errors += 1
            continue

        if r1_style[sample] != r2_style[sample]:
            log("ERROR: Sample '{}' mixes R1/R2 and 1/2 naming styles.".format(sample))
            errors += 1
            continue

        r1_src = str(Path(r1_files[sample]).resolve())
        r2_src = str(Path(r2_files[sample]).resolve())

        print(
            sample,
            "paired",
            r1_src,
            r2_src,
            "",
            sep="\\t",
            file=manifest
        )

        emitted_files.extend([r1_src, r2_src])

        log("PASS: Paired sample '{}'".format(sample))
        log("      R1: {}".format(r1_src))
        log("      R2: {}".format(r2_src))

    for sample in sorted(interleaved_files):
        if sample in r1_files or sample in r2_files:
            continue

        src = str(Path(interleaved_files[sample]).resolve())

        print(
            sample,
            "interleaved",
            "",
            "",
            src,
            sep="\\t",
            file=manifest
        )

        emitted_files.append(src)

        log("PASS: Interleaved sample '{}'".format(sample))
        log("      Interleaved: {}".format(src))

    log("")
    log("Summary:")
    log("----------------------------------------")
    log("Errors: {}".format(errors))
    log("Warnings: {}".format(warnings))
    log("Valid read files emitted: {}".format(len(emitted_files)))

    for p in emitted_files:
        log("  {}".format(p))

    if errors > 0:
        log("")
        log("FAIL: Read naming validation failed.")
        sys.exit(1)

    if not emitted_files:
        log("")
        log("FAIL: No valid read files were produced.")
        sys.exit(1)

    log("")
    log("PASS: Read naming validation completed successfully.")
PY
    """
}

process VALIDATE_READS {
    tag { read_file.simpleName }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/validation",
        mode: 'copy',
        pattern: "*_validation*.txt"

    input:
    path read_file

    output:
    tuple path(read_file), path("${read_file.simpleName}_validation.txt")

    script:
    """
    set -euo pipefail

    infile="${read_file}"
    report="${read_file.simpleName}_validation.txt"

    echo "Validation report for: ${read_file}" > "\$report"
    echo "Started: \$(date)" >> "\$report"
    echo "----------------------------------------" >> "\$report"

    case "\$infile" in
        *.fastq|*.fastq.gz|*.fq|*.fq.gz)
            echo "PASS: File extension appears valid." >> "\$report"
            ;;
        *)
            echo "ERROR: Unsupported file extension. Expected .fastq, .fastq.gz, .fq, or .fq.gz" >> "\$report"
            cat "\$report" >&2
            exit 1
            ;;
    esac

    if [[ "\$infile" == *.gz ]]; then
        READER="gzip -cd"
        echo "INFO: File is gzipped." >> "\$report"
    else
        READER="cat"
        echo "INFO: File is uncompressed." >> "\$report"
    fi

    if eval "\$READER" '"\$infile"' | awk '
    BEGIN {
        errors=0;
        rec=0;
    }
    {
        mod = NR % 4;

        if (mod == 1) {
            if (substr(\$0,1,1) != "@") {
                print "ERROR: Header line does not start with @ at line " NR;
                errors++;
            }
        }
        else if (mod == 2) {
            seq = \$0;
        }
        else if (mod == 3) {
            if (substr(\$0,1,1) != "+") {
                print "ERROR: Plus line does not start with + at line " NR;
                errors++;
            }
        }
        else if (mod == 0) {
            qual = \$0;
            rec++;

            if (length(seq) != length(qual)) {
                print "ERROR: Sequence and quality lengths differ for record ending at line " NR " (seq=" length(seq) ", qual=" length(qual) ")";
                errors++;
            }
        }
    }
    END {
        if (NR == 0) {
            print "ERROR: File is empty.";
            errors++;
        }

        if (NR % 4 != 0) {
            print "ERROR: Total number of lines is not a multiple of 4 (" NR " lines).";
            errors++;
        }

        if (errors == 0) {
            print "PASS: FASTQ structure appears valid.";
            print "INFO: Total records checked: " rec;
            exit 0;
        } else {
            print "FAIL: FASTQ validation failed with " errors " error(s).";
            exit 1;
        }
    }' >> "\$report" 2>&1
    then
        echo "Validation completed successfully." >> "\$report"
    else
        echo "Validation failed." >> "\$report"
        cat "\$report" >&2
        exit 1
    fi

    echo "Finished: \$(date)" >> "\$report"
    """
}

process SKIP_VALIDATE_READS {
    tag { read_file.simpleName }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/validation",
        mode: 'copy',
        pattern: "*_validation*.txt"

    input:
    path read_file

    output:
    tuple path(read_file), path("${read_file.simpleName}_validation_skipped.txt")

    script:
    """
    set -euo pipefail

    report="${read_file.simpleName}_validation_skipped.txt"

    echo "Validation skipped for: ${read_file}" > "\$report"
    echo "Started: \$(date)" >> "\$report"
    echo "----------------------------------------" >> "\$report"
    echo "INFO: FASTQ structure validation was skipped because --skip_validate true was used." >> "\$report"
    echo "INFO: File naming and extension checks were still performed by CHECK_READ_NAMING." >> "\$report"
    echo "INFO: This file was passed directly to FastQC." >> "\$report"
    echo "Finished: \$(date)" >> "\$report"
    """
}

process RUN_FASTQC {
    tag { read_file.simpleName }

    stageInMode 'symlink'

    publishDir "${params.output_dir ?: params.outdir}/fastqc",
        mode: 'copy'

    cpus {
        params.threads != null
            ? params.threads as int
            : params.fastqc_threads as int
    }

    input:
    tuple path(read_file), path(validation_report), path(tools_status)

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

    echo "Running FastQC on: ${read_file}"
    echo "Validation/skipping report: ${validation_report}"
    echo "FastQC threads: ${task.cpus}"

    fastqc \\
        -t ${task.cpus} \\
        "${read_file}" \\
        --outdir .
    """
}