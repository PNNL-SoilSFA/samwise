#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir    = params.input_dir ?: null
params.outdir       = params.outdir ?: "./results/module_0_readprocess"
params.auto_install = params.auto_install ?: true

/*
 * Use "*" so the workflow sees every file in the directory.
 * This allows it to catch non-FASTQ files or badly named FASTQ files.
 */
params.file_pattern = params.file_pattern ?: "*"

workflow {

    if( !params.input_dir ) {
        error """
        Missing required parameter: --input_dir

        Example:
          nextflow run module_0_readprocess.nf --input_dir ./reads
        """.stripIndent()
    }

    /*
     * Collect all files from the input directory.
     */
    all_files_ch = channel.fromPath(
        "${params.input_dir}/${params.file_pattern}",
        type: 'file',
        checkIfExists: true
    )

    /*
     * FastQC installation/check can happen independently,
     * but RUN_FASTQC will be forced to wait for it later.
     */
    INSTALL_FASTQC()

    /*
     * Check naming conventions, identify pairs/interleaved files,
     * and normalize .fq/.fq.gz to .fastq/.fastq.gz using symlinks.
     */
    CHECK_READ_NAMING(all_files_ch.collect())

    /*
     * The naming process emits valid normalized read files under valid_reads/.
     * Flatten is used because the process emits multiple files as a list.
     */
    valid_named_reads_ch = CHECK_READ_NAMING.out.reads.flatten()

    /*
     * Validate FASTQ structure only after naming normalization.
     */
    VALIDATE_READS(valid_named_reads_ch)

    /*
     * Force FastQC to wait until INSTALL_FASTQC completes.
     * combine() pairs every validated read with the FastQC install status file.
     */
    reads_ready_for_fastqc_ch = VALIDATE_READS.out.combine(INSTALL_FASTQC.out)

    RUN_FASTQC(reads_ready_for_fastqc_ch)
}


process INSTALL_FASTQC {

    tag "check_fastqc"

    publishDir "${params.outdir}/setup", mode: 'copy'

    output:
    path "fastqc_install_status.txt"

    script:
    """
    set -euo pipefail

    STATUS_FILE="fastqc_install_status.txt"

    if command -v fastqc >/dev/null 2>&1; then
        echo "FastQC already installed: \$(command -v fastqc)" > "\$STATUS_FILE"
        fastqc --version >> "\$STATUS_FILE" 2>&1 || true
        exit 0
    fi

    echo "FastQC not found in PATH." > "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "Auto-install disabled. Please install FastQC manually or run with --auto_install true" >> "\$STATUS_FILE"
        exit 1
    fi

    if ! command -v mamba >/dev/null 2>&1; then
        echo "mamba not found. Please install mamba first." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "Attempting to install FastQC using mamba..." >> "\$STATUS_FILE"

    if mamba install -y -c bioconda -c conda-forge fastqc >> "\$STATUS_FILE" 2>&1; then
        echo "FastQC installation successful." >> "\$STATUS_FILE"
    else
        echo "FastQC installation failed." >> "\$STATUS_FILE"
        exit 1
    fi

    if command -v fastqc >/dev/null 2>&1; then
        echo "FastQC path: \$(command -v fastqc)" >> "\$STATUS_FILE"
        fastqc --version >> "\$STATUS_FILE" 2>&1 || true
    else
        echo "FastQC still not detected after installation." >> "\$STATUS_FILE"
        exit 1
    fi
    """
}


process CHECK_READ_NAMING {
            log(f"ERROR: Sample '{sample}' mixes R1/R2 and 1/2 naming styles.")
            errors += 1
            continue

        if r1_style[sample] == "1" and r2_style[sample] != "2":
            log(f"ERROR: Sample '{sample}' mixes 1/2 and R1/R2 naming styles.")
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

        log(f"PASS: Paired sample '{sample}'")
        log(f"      R1: {Path(r1_src).name} -> {r1_dest}")
        log(f"      R2: {Path(r2_src).name} -> {r2_dest}")

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

        log(f"PASS: Interleaved sample '{sample}'")
        log(f"      Interleaved: {Path(src).name} -> {dest}")

    valid_outputs = list(valid_reads_dir.glob("*"))

    log("")
    log("Summary:")
    log("----------------------------------------")
    log(f"Errors: {errors}")
    log(f"Warnings: {warnings}")
    log(f"Valid read files emitted: {len(valid_outputs)}")

    for p in valid_outputs:
        log(f"  {p}")

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

    tag { file.simpleName }

    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path file

    output:
    tuple path(file), path("${file.simpleName}_validation.txt")

    script:
    """
    set -euo pipefail

    infile="${file}"
    report="${file.simpleName}_validation.txt"

    echo "Validation report for: ${file}" > "\$report"
    echo "Started: \$(date)" >> "\$report"
    echo "----------------------------------------" >> "\$report"

    /*
     * CHECK_READ_NAMING normalizes .fq/.fq.gz into .fastq/.fastq.gz.
     * Therefore, VALIDATE_READS should only receive .fastq or .fastq.gz.
     */
    case "\$infile" in
        *.fastq|*.fastq.gz)
            echo "PASS: File extension appears valid after normalization." >> "\$report"
            ;;
        *)
            echo "ERROR: Unsupported file extension after normalization. Expected .fastq or .fastq.gz" >> "\$report"
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


process RUN_FASTQC {

    tag { file.simpleName }

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    cpus params.fastqc_threads

    input:
    tuple path(file), path(report), path(fastqc_status)

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

    echo "Running FastQC on: ${file}"
    echo "FastQC threads: ${task.cpus}"

    fastqc \\
        -t ${task.cpus} \\
        "${file}" \\
        --outdir .
    """
}
