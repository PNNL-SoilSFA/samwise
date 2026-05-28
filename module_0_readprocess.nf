#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Default parameters.
 * Command-line values override these.
 */

params.input_dir       = null
params.outdir          = "./results/module_0_readprocess"
params.output_dir      = null
params.auto_install    = true
params.file_pattern    = "*"
params.fastqc_threads  = 2
params.threads         = null
params.skip_validate   = false


workflow {
    if( !params.input_dir ) {
        error """
        Missing required parameter: --input_dir
        Example:
          nextflow run module_0_readprocess.nf --input_dir ./reads
        """.stripIndent()
    }

    all_files_ch = channel.fromPath(
        "${params.input_dir}/${params.file_pattern}",
        type: 'file',
        checkIfExists: true
    )

    INSTALL_FASTQC()

    CHECK_READ_NAMING(all_files_ch.collect())

    valid_named_reads_ch = CHECK_READ_NAMING.out.reads.flatten()

    if( params.skip_validate.toString() == 'true' ) {
        SKIP_VALIDATE_READS(valid_named_reads_ch)
        reads_ready_for_fastqc_ch = SKIP_VALIDATE_READS.out.combine(INSTALL_FASTQC.out)
    }
    else {
        VALIDATE_READS(valid_named_reads_ch)
        reads_ready_for_fastqc_ch = VALIDATE_READS.out.combine(INSTALL_FASTQC.out)
    }

    RUN_FASTQC(reads_ready_for_fastqc_ch)
}


process INSTALL_FASTQC {
    tag "check_fastqc"
    publishDir "${params.output_dir ?: params.outdir}/setup", mode: 'copy'

    output:
    path "fastqc_install_status.txt"

    script:
    """
    set -euo pipefail

    STATUS_FILE="fastqc_install_status.txt"

    echo "FastQC setup/check started: \$(date)" > "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    if command -v fastqc >/dev/null 2>&1; then
        echo "FastQC already installed: \$(command -v fastqc)" >> "\$STATUS_FILE"
        fastqc --version >> "\$STATUS_FILE" 2>&1 || true
        exit 0
    fi

    echo "FastQC not found in PATH." >> "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "Auto-install disabled." >> "\$STATUS_FILE"
        echo "Please install FastQC manually or run with --auto_install true" >> "\$STATUS_FILE"
        exit 1
    fi

    INSTALLER=""

    if command -v mamba >/dev/null 2>&1; then
        INSTALLER="mamba"
        echo "mamba detected: \$(command -v mamba)" >> "\$STATUS_FILE"
    elif command -v conda >/dev/null 2>&1; then
        INSTALLER="conda"
        echo "mamba not found." >> "\$STATUS_FILE"
        echo "conda detected: \$(command -v conda)" >> "\$STATUS_FILE"
    else
        echo "Neither mamba nor conda was found in PATH." >> "\$STATUS_FILE"
        echo "Please install FastQC manually, or install mamba/conda first." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "Attempting to install FastQC using \$INSTALLER..." >> "\$STATUS_FILE"

    if "\$INSTALLER" install -y -c conda-forge -c bioconda fastqc >> "\$STATUS_FILE" 2>&1; then
        echo "FastQC installation command completed successfully." >> "\$STATUS_FILE"
    else
        echo "FastQC installation failed using \$INSTALLER." >> "\$STATUS_FILE"
        exit 1
    fi

    hash -r || true

    if command -v fastqc >/dev/null 2>&1; then
        echo "FastQC path after installation: \$(command -v fastqc)" >> "\$STATUS_FILE"
        fastqc --version >> "\$STATUS_FILE" 2>&1 || true
    else
        echo "FastQC still not detected after installation." >> "\$STATUS_FILE"
        echo "The install may have succeeded, but the environment PATH may not have updated." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "FastQC setup/check finished: \$(date)" >> "\$STATUS_FILE"
    """
}

process CHECK_READ_NAMING {

    tag "check_read_naming"

    publishDir "${params.output_dir ?: params.outdir}/naming", mode: 'copy', pattern: "*.{txt,tsv}"

    input:
    path read_files

    output:
    path "read_naming_report.txt", emit: report
    path "read_manifest.tsv", emit: manifest
    path "valid_reads/*", emit: reads

    script:
    """
    set -euo pipefail

    mkdir -p valid_reads

    cat > staged_files.list <<'EOF'
${read_files.join('\n')}
EOF

    python3 - <<'PY'
import os
import re
import sys
from pathlib import Path

report_path = Path("read_naming_report.txt")
manifest_path = Path("read_manifest.tsv")
valid_reads_dir = Path("valid_reads")
valid_reads_dir.mkdir(exist_ok=True)

r1_files = {}
r2_files = {}
interleaved_files = {}
r1_style = {}
r2_style = {}

errors = 0
warnings = 0

read1_re = re.compile(r'^(.+)_(R1|1)\\.(fastq|fq)(\\.gz)?\\Z')
read2_re = re.compile(r'^(.+)_(R2|2)\\.(fastq|fq)(\\.gz)?\\Z')
interleaved_re = re.compile(r'^(.+)_interleaved\\.(fastq|fq)(\\.gz)?\\Z')
fastq_like_re = re.compile(r'.*\\.(fastq|fq)(\\.gz)?\\Z')


def normalize_read_name(filename):
    if filename.endswith(".fastq.gz"):
        return filename
    if filename.endswith(".fq.gz"):
        return filename[:-6] + ".fastq.gz"
    if filename.endswith(".fastq"):
        return filename
    if filename.endswith(".fq"):
        return filename[:-3] + ".fastq"
    return filename


def safe_symlink(src, dest):
    src_real = os.path.realpath(src)
    dest = Path(dest)

    if dest.exists() or dest.is_symlink():
        dest.unlink()

    os.symlink(src_real, dest)


with report_path.open("w") as report, manifest_path.open("w") as manifest:

    def log(message=""):
        print(message, file=report)

    print("sample_id", "layout", "read1", "read2", "interleaved", sep="\\t", file=manifest)

    log("Read naming and pairing report")
    log("Started")
    log("----------------------------------------")
    log("")

    with open("staged_files.list") as handle:
        staged_files = [line.strip() for line in handle if line.strip()]

    if not staged_files:
        log("ERROR: No files were found in the input directory.")
        errors += 1

    log("Files detected:")
    for f in staged_files:
        log("  {}".format(f))
    log("")

    for file_path in staged_files:
        path = Path(file_path)
        name = path.name

        m1 = read1_re.match(name)
        m2 = read2_re.match(name)
        mi = interleaved_re.match(name)

        if mi:
            sample = mi.group(1)

            if sample in interleaved_files:
                log("ERROR: Duplicate interleaved file for sample '{}': {}".format(sample, name))
                errors += 1

            interleaved_files[sample] = file_path
            continue

        if m1:
            sample = m1.group(1)
            style = m1.group(2)

            if sample in r1_files:
                log("ERROR: Duplicate R1 file for sample '{}': {}".format(sample, name))
                errors += 1

            r1_files[sample] = file_path
            r1_style[sample] = style
            continue

        if m2:
            sample = m2.group(1)
            style = m2.group(2)

            if sample in r2_files:
                log("ERROR: Duplicate R2 file for sample '{}': {}".format(sample, name))
                errors += 1

            r2_files[sample] = file_path
            r2_style[sample] = style
            continue

        if fastq_like_re.match(name):
            log("ERROR: FASTQ file has unsupported naming convention: {}".format(name))
            log("       Expected:")
            log("         sample_R1.fastq or sample_R1.fastq.gz")
            log("         sample_R2.fastq or sample_R2.fastq.gz")
            log("         sample_1.fastq  or sample_1.fastq.gz")
            log("         sample_2.fastq  or sample_2.fastq.gz")
            log("         sample_interleaved.fastq or sample_interleaved.fastq.gz")
            errors += 1
        else:
            log("ERROR: Non-FASTQ file detected: {}".format(name))
            errors += 1

    log("")
    log("Pair/interleaved checks:")
    log("----------------------------------------")

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

        if r1_style[sample] == "R1" and r2_style[sample] != "R2":
            log("ERROR: Sample '{}' mixes R1/R2 and 1/2 naming styles.".format(sample))
            errors += 1
            continue

        if r1_style[sample] == "1" and r2_style[sample] != "2":
            log("ERROR: Sample '{}' mixes 1/2 and R1/R2 naming styles.".format(sample))
            errors += 1
            continue

        r1_src = r1_files[sample]
        r2_src = r2_files[sample]

        r1_name = normalize_read_name(Path(r1_src).name)
        r2_name = normalize_read_name(Path(r2_src).name)

        r1_dest = valid_reads_dir / r1_name
        r2_dest = valid_reads_dir / r2_name

        safe_symlink(r1_src, r1_dest)
        safe_symlink(r2_src, r2_dest)

        print(
            sample,
            "paired",
            str(r1_dest),
            str(r2_dest),
            "",
            sep="\\t",
            file=manifest
        )

        log("PASS: Paired sample '{}'".format(sample))
        log("      R1: {} -> {}".format(Path(r1_src).name, r1_dest))
        log("      R2: {} -> {}".format(Path(r2_src).name, r2_dest))

    for sample in sorted(interleaved_files):
        if sample in r1_files or sample in r2_files:
            continue

        src = interleaved_files[sample]
        norm_name = normalize_read_name(Path(src).name)
        dest = valid_reads_dir / norm_name

        safe_symlink(src, dest)

        print(
            sample,
            "interleaved",
            "",
            "",
            str(dest),
            sep="\\t",
            file=manifest
        )

        log("PASS: Interleaved sample '{}'".format(sample))
        log("      Interleaved: {} -> {}".format(Path(src).name, dest))

    valid_outputs = list(valid_reads_dir.glob("*"))

    log("")
    log("Summary:")
    log("----------------------------------------")
    log("Errors: {}".format(errors))
    log("Warnings: {}".format(warnings))
    log("Valid read files emitted: {}".format(len(valid_outputs)))

    for p in valid_outputs:
        log("  {}".format(p))

    if errors > 0:
        log("")
        log("FAIL: Read naming validation failed.")
        sys.exit(1)

    if not valid_outputs:
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

    publishDir "${params.output_dir ?: params.outdir}/validation", mode: 'copy'

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
        *.fastq|*.fastq.gz)
            echo "PASS: File extension appears valid after normalization." >> "\$report"
            ;;
        *)
            echo "ERROR: Unsupported file extension after normalization. Expected .fastq or .fastq.gz" >> "\$report"
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

    publishDir "${params.output_dir ?: params.outdir}/validation", mode: 'copy'

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

    publishDir "${params.output_dir ?: params.outdir}/fastqc", mode: 'copy'

    cpus { (params.threads ?: params.fastqc_threads) as int }

    input:
    tuple path(read_file), path(validation_report), path(fastqc_status)

    output:
    path "*_fastqc.html"
    path "*_fastqc.zip"

    script:
    """
    set -euo pipefail

    if ! command -v fastqc >/dev/null 2>&1; then
        echo "ERROR: FastQC is not available in PATH."
        echo "FastQC install/check status file:"
        cat "${fastqc_status}" || true
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