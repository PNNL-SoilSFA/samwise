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

    tag "check_read_naming"

    publishDir "${params.outdir}/naming", mode: 'copy', pattern: "*.{txt,tsv}"

    input:
    path read_files

    output:
    path "read_naming_report.txt", emit: report
    path "read_manifest.tsv", emit: manifest
    path "valid_reads/*", emit: reads

    script:
    """
    set -euo pipefail

    report="read_naming_report.txt"
    manifest="read_manifest.tsv"

    mkdir -p valid_reads

    echo "Read naming and pairing report" > "\$report"
    echo "Started: \$(date)" >> "\$report"
    echo "----------------------------------------" >> "\$report"
    echo "" >> "\$report"

    echo -e "sample_id\\tlayout\\tread1\\tread2\\tinterleaved" > "\$manifest"

    cat > staged_files.list <<'EOF'
${read_files.join('\n')}
EOF

    declare -A r1_files
    declare -A r2_files
    declare -A interleaved_files
    declare -A r1_style
    declare -A r2_style
    declare -A seen_keys

    errors=0
    warnings=0

    normalize_read_name() {
        local filename="\$1"

        if [[ "\$filename" == *.fastq.gz ]]; then
            echo "\$filename"
        elif [[ "\$filename" == *.fq.gz ]]; then
            echo "\${filename%.fq.gz}.fastq.gz"
        elif [[ "\$filename" == *.fastq ]]; then
            echo "\$filename"
        elif [[ "\$filename" == *.fq ]]; then
            echo "\${filename%.fq}.fastq"
        else
            echo "\$filename"
        fi
    }

    while IFS= read -r f; do

        if [[ -z "\$f" ]]; then
            continue
        fi

        if [[ ! -f "\$f" ]]; then
            continue
        fi

        base=\$(basename "\$f")

        echo "Checking file: \$base" >> "\$report"

        /*
         * Valid read 1 names:
         *   SampleID_R1.fastq
         *   SampleID_R1.fastq.gz
         *   SampleID_R1.fq
         *   SampleID_R1.fq.gz
         *   SampleID_1.fastq
         *   SampleID_1.fastq.gz
         *   SampleID_1.fq
         *   SampleID_1.fq.gz
         */
        if [[ "\$base" =~ ^(.+)_(R1|1)\\.(fastq|fq)(\\.gz)?\$ ]]; then

            sample="\${BASH_REMATCH[1]}"
            style="\${BASH_REMATCH[2]}"
            ext="\${BASH_REMATCH[3]}"
            key="\${sample}:R1"

            if [[ -n "\${seen_keys[\$key]:-}" ]]; then
                echo "ERROR: Duplicate read 1 file detected for sample '\$sample'." >> "\$report"
                echo "       Only one read 1 file should exist for each sample." >> "\$report"
                echo "       Accepted read 1 suffixes:" >> "\$report"
                echo "         _R1.fastq, _R1.fastq.gz, _R1.fq, _R1.fq.gz" >> "\$report"
                echo "         _1.fastq,  _1.fastq.gz,  _1.fq,  _1.fq.gz" >> "\$report"
                errors=\$((errors + 1))
            else
                seen_keys[\$key]=1
                r1_files[\$sample]="\$f"
                r1_style[\$sample]="\$style"

                if [[ "\$ext" == "fq" ]]; then
                    normalized=\$(normalize_read_name "\$base")
                    echo "PASS: Detected read 1 for sample '\$sample' using suffix '_\$style'." >> "\$report"
                    echo "INFO: File uses .fq extension and will be normalized downstream:" >> "\$report"
                    echo "      \$base -> \$normalized" >> "\$report"
                else
                    echo "PASS: Detected read 1 for sample '\$sample' using suffix '_\$style'." >> "\$report"
                fi
            fi

        /*
         * Valid read 2 names:
         *   SampleID_R2.fastq
         *   SampleID_R2.fastq.gz
         *   SampleID_R2.fq
         *   SampleID_R2.fq.gz
         *   SampleID_2.fastq
         *   SampleID_2.fastq.gz
         *   SampleID_2.fq
         *   SampleID_2.fq.gz
         */
        elif [[ "\$base" =~ ^(.+)_(R2|2)\\.(fastq|fq)(\\.gz)?\$ ]]; then

            sample="\${BASH_REMATCH[1]}"
            style="\${BASH_REMATCH[2]}"
            ext="\${BASH_REMATCH[3]}"
            key="\${sample}:R2"

            if [[ -n "\${seen_keys[\$key]:-}" ]]; then
                echo "ERROR: Duplicate read 2 file detected for sample '\$sample'." >> "\$report"
                echo "       Only one read 2 file should exist for each sample." >> "\$report"
                echo "       Accepted read 2 suffixes:" >> "\$report"
                echo "         _R2.fastq, _R2.fastq.gz, _R2.fq, _R2.fq.gz" >> "\$report"
                echo "         _2.fastq,  _2.fastq.gz,  _2.fq,  _2.fq.gz" >> "\$report"
                errors=\$((errors + 1))
            else
                seen_keys[\$key]=1
                r2_files[\$sample]="\$f"
                r2_style[\$sample]="\$style"

                if [[ "\$ext" == "fq" ]]; then
                    normalized=\$(normalize_read_name "\$base")
                    echo "PASS: Detected read 2 for sample '\$sample' using suffix '_\$style'." >> "\$report"
                    echo "INFO: File uses .fq extension and will be normalized downstream:" >> "\$report"
                    echo "      \$base -> \$normalized" >> "\$report"
                else
                    echo "PASS: Detected read 2 for sample '\$sample' using suffix '_\$style'." >> "\$report"
                fi
            fi

        /*
         * Valid interleaved names:
         *   SampleID_interleaved.fastq
         *   SampleID_interleaved.fastq.gz
         *   SampleID_interleaved.fq
         *   SampleID_interleaved.fq.gz
         */
        elif [[ "\$base" =~ ^(.+)_interleaved\\.(fastq|fq)(\\.gz)?\$ ]]; then

            sample="\${BASH_REMATCH[1]}"
            ext="\${BASH_REMATCH[2]}"
            key="\${sample}:interleaved"

            if [[ -n "\${seen_keys[\$key]:-}" ]]; then
                echo "ERROR: Duplicate interleaved file detected for sample '\$sample'." >> "\$report"
                errors=\$((errors + 1))
            else
                seen_keys[\$key]=1
                interleaved_files[\$sample]="\$f"

                if [[ "\$ext" == "fq" ]]; then
                    normalized=\$(normalize_read_name "\$base")
                    echo "PASS: Detected interleaved reads for sample '\$sample'." >> "\$report"
                    echo "INFO: File uses .fq extension and will be normalized downstream:" >> "\$report"
                    echo "      \$base -> \$normalized" >> "\$report"
                else
                    echo "PASS: Detected interleaved reads for sample '\$sample'." >> "\$report"
                fi
            fi

        elif [[ "\$base" =~ \\.(fastq|fq)(\\.gz)?\$ ]]; then

            echo "ERROR: FASTQ-like file has unsupported name: \$base" >> "\$report"
            echo "       Expected one of:" >> "\$report"
            echo "         SampleID_R1.fastq" >> "\$report"
            echo "         SampleID_R1.fastq.gz" >> "\$report"
            echo "         SampleID_R1.fq" >> "\$report"
            echo "         SampleID_R1.fq.gz" >> "\$report"
            echo "         SampleID_R2.fastq" >> "\$report"
            echo "         SampleID_R2.fastq.gz" >> "\$report"
            echo "         SampleID_R2.fq" >> "\$report"
            echo "         SampleID_R2.fq.gz" >> "\$report"
            echo "         SampleID_1.fastq" >> "\$report"
            echo "         SampleID_1.fastq.gz" >> "\$report"
            echo "         SampleID_1.fq" >> "\$report"
            echo "         SampleID_1.fq.gz" >> "\$report"
            echo "         SampleID_2.fastq" >> "\$report"
            echo "         SampleID_2.fastq.gz" >> "\$report"
            echo "         SampleID_2.fq" >> "\$report"
            echo "         SampleID_2.fq.gz" >> "\$report"
            echo "         SampleID_interleaved.fastq" >> "\$report"
            echo "         SampleID_interleaved.fastq.gz" >> "\$report"
            echo "         SampleID_interleaved.fq" >> "\$report"
            echo "         SampleID_interleaved.fq.gz" >> "\$report"
            errors=\$((errors + 1))

        else

            echo "ERROR: File does not appear to be a supported FASTQ file: \$base" >> "\$report"
            echo "       Supported extensions are .fastq, .fastq.gz, .fq, and .fq.gz." >> "\$report"
            errors=\$((errors + 1))

        fi

        echo "" >> "\$report"

    done < staged_files.list

    {
        printf "%s\\n" "\${!r1_files[@]}"
        printf "%s\\n" "\${!r2_files[@]}"
        printf "%s\\n" "\${!interleaved_files[@]}"
    } | sort -u > samples.tmp

    if [[ ! -s samples.tmp ]]; then
        echo "ERROR: No valid read files were found." >> "\$report"
        errors=\$((errors + 1))
    fi

    echo "" >> "\$report"
    echo "Pairing/interleaving checks" >> "\$report"
    echo "----------------------------------------" >> "\$report"

    while IFS= read -r sample; do

        has_r1=0
        has_r2=0
        has_interleaved=0

        if [[ -n "\${r1_files[\$sample]:-}" ]]; then
            has_r1=1
        fi

        if [[ -n "\${r2_files[\$sample]:-}" ]]; then
            has_r2=1
        fi

        if [[ -n "\${interleaved_files[\$sample]:-}" ]]; then
            has_interleaved=1
        fi

        if [[ "\$has_interleaved" -eq 1 ]] && { [[ "\$has_r1" -eq 1 ]] || [[ "\$has_r2" -eq 1 ]]; }; then

            echo "ERROR: Sample '\$sample' has both interleaved and paired read files." >> "\$report"
            echo "       Please provide either paired reads:" >> "\$report"
            echo "         \${sample}_R1.fastq.gz and \${sample}_R2.fastq.gz" >> "\$report"
            echo "       or:" >> "\$report"
            echo "         \${sample}_1.fastq.gz and \${sample}_2.fastq.gz" >> "\$report"
            echo "       or an interleaved file:" >> "\$report"
            echo "         \${sample}_interleaved.fastq.gz" >> "\$report"
            echo "       but not both." >> "\$report"
            errors=\$((errors + 1))

        elif [[ "\$has_interleaved" -eq 1 ]]; then

            interleaved="\${interleaved_files[\$sample]}"
            interleaved_base=\$(basename "\$interleaved")
            interleaved_norm=\$(normalize_read_name "\$interleaved_base")
            interleaved_real=\$(realpath "\$interleaved")

            echo "PASS: Sample '\$sample' is interleaved." >> "\$report"

            if [[ "\$interleaved_base" != "\$interleaved_norm" ]]; then
                echo "INFO: Normalized interleaved filename:" >> "\$report"
                echo "      \$interleaved_base -> \$interleaved_norm" >> "\$report"
            fi

            echo -e "\$sample\\tinterleaved\\tNA\\tNA\\t\$interleaved_norm" >> "\$manifest"
            ln -sf "\$interleaved_real" "valid_reads/\$interleaved_norm"

        elif [[ "\$has_r1" -eq 1 && "\$has_r2" -eq 1 ]]; then

            r1="\${r1_files[\$sample]}"
            r2="\${r2_files[\$sample]}"

            r1_base=\$(basename "\$r1")
            r2_base=\$(basename "\$r2")

            r1_norm=\$(normalize_read_name "\$r1_base")
            r2_norm=\$(normalize_read_name "\$r2_base")

            r1_real=\$(realpath "\$r1")
            r2_real=\$(realpath "\$r2")

            if [[ "\${r1_style[\$sample]}" == "R1" && "\${r2_style[\$sample]}" == "2" ]]; then
                echo "WARNING: Sample '\$sample' mixes naming styles: _R1 with _2." >> "\$report"
                echo "         This is allowed, but consider renaming consistently." >> "\$report"
                warnings=\$((warnings + 1))
            elif [[ "\${r1_style[\$sample]}" == "1" && "\${r2_style[\$sample]}" == "R2" ]]; then
                echo "WARNING: Sample '\$sample' mixes naming styles: _1 with _R2." >> "\$report"
                echo "         This is allowed, but consider renaming consistently." >> "\$report"
                warnings=\$((warnings + 1))
            fi

            echo "PASS: Sample '\$sample' has a valid read 1/read 2 pair." >> "\$report"

            if [[ "\$r1_base" != "\$r1_norm" ]]; then
                echo "INFO: Normalized read 1 filename:" >> "\$report"
                echo "      \$r1_base -> \$r1_norm" >> "\$report"
            fi

            if [[ "\$r2_base" != "\$r2_norm" ]]; then
                echo "INFO: Normalized read 2 filename:" >> "\$report"
                echo "      \$r2_base -> \$r2_norm" >> "\$report"
            fi

            echo -e "\$sample\\tpaired\\t\$r1_norm\\t\$r2_norm\\tNA" >> "\$manifest"

            ln -sf "\$r1_real" "valid_reads/\$r1_norm"
            ln -sf "\$r2_real" "valid_reads/\$r2_norm"

        elif [[ "\$has_r1" -eq 1 && "\$has_r2" -eq 0 ]]; then

            echo "ERROR: Sample '\$sample' has read 1 but is missing read 2." >> "\$report"
            echo "       Expected one of:" >> "\$report"
            echo "         \${sample}_R2.fastq" >> "\$report"
            echo "         \${sample}_R2.fastq.gz" >> "\$report"
            echo "         \${sample}_R2.fq" >> "\$report"
            echo "         \${sample}_R2.fq.gz" >> "\$report"
            echo "         \${sample}_2.fastq" >> "\$report"
            echo "         \${sample}_2.fastq.gz" >> "\$report"
            echo "         \${sample}_2.fq" >> "\$report"
            echo "         \${sample}_2.fq.gz" >> "\$report"
            errors=\$((errors + 1))

        elif [[ "\$has_r1" -eq 0 && "\$has_r2" -eq 1 ]]; then

            echo "ERROR: Sample '\$sample' has read 2 but is missing read 1." >> "\$report"
            echo "       Expected one of:" >> "\$report"
            echo "         \${sample}_R1.fastq" >> "\$report"
            echo "         \${sample}_R1.fastq.gz" >> "\$report"
            echo "         \${sample}_R1.fq" >> "\$report"
            echo "         \${sample}_R1.fq.gz" >> "\$report"
            echo "         \${sample}_1.fastq" >> "\$report"
            echo "         \${sample}_1.fastq.gz" >> "\$report"
            echo "         \${sample}_1.fq" >> "\$report"
            echo "         \${sample}_1.fq.gz" >> "\$report"
            errors=\$((errors + 1))

        else

            echo "ERROR: Unexpected read structure for sample '\$sample'." >> "\$report"
            errors=\$((errors + 1))

        fi

    done < samples.tmp

    echo "" >> "\$report"
    echo "Finished: \$(date)" >> "\$report"

    if [[ "\$errors" -gt 0 ]]; then
        echo "" >> "\$report"
        echo "FAIL: Read naming check failed with \$errors error(s) and \$warnings warning(s)." >> "\$report"
        echo "" >&2
        echo "Read naming check failed. Report:" >&2
        echo "----------------------------------------" >&2
        cat "\$report" >&2
        exit 1
    else
        echo "" >> "\$report"
        echo "PASS: All read names and pairings are valid." >> "\$report"
        echo "WARNINGS: \$warnings" >> "\$report"
    fi
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

    fastqc "${file}" --outdir .
    """
}