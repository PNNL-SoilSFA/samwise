#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Module 5: Subtractive assembly + second-pass binning + final joint MAG refinement.
 */

params.working_dir = null
params.output_dir = null
params.input_trimmed_manifest = null
params.input_original_binning_manifest = null
params.input_refined_manifest = null
params.auto_install = true
params.tool_env_dir = null
params.threads = null
params.mapping_threads = 4
params.assembly_threads = 4
params.bbmap_version = "39.81"
params.megahit = true
params.metaspades = false
params.megahit_version = "1.2.9"
params.spades_version = "4.2.0"
params.metaspades_memory_gb = 0
params.bbmap_extra_args = ""
params.bbmap_minid = 0.99
params.bbmap_ambig = "random"
params.bbmap_xmx = null
params.megahit_preset = "meta-large"
params.megahit_threads = null
params.run_second_pass_binning_refinement = true
params.secondpass_metabat2 = true
params.secondpass_quickbin = true
params.secondpass_maxbin2 = true
params.module3_script = "${projectDir}/module_3_binning.nf"
params.module4_script = "${projectDir}/module_4_binrefinement.nf"
params.nextflow_exe = "nextflow"
params.secondpass_working_dir = null
params.final_joint_working_dir = null
params.dependencies_dir = "${projectDir}/dependencies"
params.tigrfam_hmm = null
params.pfam_hmm = null
params.magscot_script = null
params.magscot_extra_args = ""
params.publish_reference_mode = "copy"
params.publish_final_mags_mode = "copy"
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.module1_outdir = "${params.results_dir}/module_1_readtrimming"
params.module3_outdir = "${params.results_dir}/module_3_binning"
params.module4_outdir = "${params.results_dir}/module_4_binrefinement"
params.outdir = "${params.results_dir}/module_5_subassembly"
params.secondpass_dir = params.secondpass_working_dir ?: "${params.outdir}/second_pass_binning"
params.final_joint_dir = params.final_joint_working_dir ?: "${params.outdir}/final_joint_refinement"

def absOrEmpty(value) {
    def s = value == null ? "" : value.toString().trim()

    if (!s || s == "null" || s == "NA") {
        return ""
    }

    return java.nio.file.Paths
        .get(s)
        .toAbsolutePath()
        .normalize()
        .toString()
}


workflow {

    def do_secondpass = params.run_second_pass_binning_refinement.toString().toBoolean()
    def use_megahit = params.megahit.toString().toBoolean()
    def use_metaspades = params.metaspades.toString().toBoolean()

    if (!use_megahit && !use_metaspades) {
        error(
            """
        No subtractive assembler selected.

        Please specify at least one of:
          --megahit true
          --metaspades true
        """.stripIndent()
        )
    }

    def use_secondpass_metabat2 = params.secondpass_metabat2.toString().toBoolean()
    def use_secondpass_quickbin = params.secondpass_quickbin.toString().toBoolean()
    def use_secondpass_maxbin2 = params.secondpass_maxbin2.toString().toBoolean()

    if (do_secondpass && !use_secondpass_metabat2 && !use_secondpass_quickbin && !use_secondpass_maxbin2) {
        error(
            """
        Second-pass binning/final refinement is enabled, but no second-pass binner is selected.

        Enable at least one:
          --secondpass_metabat2 true
          --secondpass_quickbin true
          --secondpass_maxbin2 true

        Or disable second pass:
          --run_second_pass_binning_refinement false
        """.stripIndent()
        )
    }

    def trimmed_manifest_file = params.input_trimmed_manifest ?: "${params.module1_outdir}/summary/trimmed_manifest.tsv"
    def original_binning_manifest_file = params.input_original_binning_manifest ?: "${params.module3_outdir}/summary/binning_manifest.tsv"
    def refined_manifest_file = params.input_refined_manifest ?: "${params.module4_outdir}/summary/magscot_refined_bins_manifest.tsv"

    log.info("Module 5 results directory: ${params.results_dir}")
    log.info("Using Module 1 trimmed manifest: ${trimmed_manifest_file}")
    log.info("Using original Module 3 binning manifest: ${original_binning_manifest_file}")
    log.info("Using Module 4 refined MAG manifest as subtractive reference: ${refined_manifest_file}")
    log.info("Writing Module 5 outputs to: ${params.outdir}")
    log.info("Subtractive assemblers: MEGAHIT=${use_megahit}, metaSPAdes=${use_metaspades}")
    log.info("Second-pass enabled: ${do_secondpass}")
    log.info("Second-pass binners: MetaBAT2=${use_secondpass_metabat2}, QuickBin=${use_secondpass_quickbin}, MaxBin2=${use_secondpass_maxbin2}")
    log.info("Second-pass working directory: ${params.secondpass_dir}")
    log.info("Final joint refinement working directory: ${params.final_joint_dir}")

    def assembler_list = []

    if (use_megahit) {
        assembler_list << "megahit"
    }

    if (use_metaspades) {
        assembler_list << "metaspades"
    }

    def trimmed_manifest_ch = channel.fromPath(
        trimmed_manifest_file,
        type: 'file',
        checkIfExists: true,
    )

    def original_binning_manifest_ch = channel.fromPath(
        original_binning_manifest_file,
        type: 'file',
        checkIfExists: true,
    )

    def refined_manifest_ch = channel.fromPath(
        refined_manifest_file,
        type: 'file',
        checkIfExists: true,
    )

    def reads_ch = trimmed_manifest_ch
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row.sample_id.toString()
            def safe_id = row.safe_sample_id.toString()
            def layout = row.layout.toString()

            def assembly_sample_id = sample_id.replaceAll('[^A-Za-z0-9]+', '')

            if (!assembly_sample_id) {
                assembly_sample_id = safe_id.replaceAll('[^A-Za-z0-9]+', '')
            }

            if (!assembly_sample_id) {
                error("Could not derive non-empty assembly SampleID from sample '${sample_id}'")
            }

            if (layout != "paired" && layout != "interleaved") {
                error("Unsupported layout in trimmed manifest for sample '${sample_id}': ${layout}")
            }

            def read1 = absOrEmpty(row.read1)
            def read2 = absOrEmpty(row.read2)
            def interleaved = absOrEmpty(row.interleaved)
            def read_files

            if (layout == "paired") {
                if (!read1 || !read2) {
                    error("Paired layout requires read1 and read2 for sample '${sample_id}'")
                }

                read_files = [
                    file(read1, checkIfExists: true),
                    file(read2, checkIfExists: true),
                ]
            } else {
                if (!interleaved) {
                    error("Interleaved layout requires an interleaved read file for sample '${sample_id}'")
                }

                read_files = [file(interleaved, checkIfExists: true)]
            }

            tuple(
                sample_id,
                safe_id,
                assembly_sample_id,
                layout,
                read_files,
            )
        }

    CLEAN_MODULE5_DURABLE_SUBTRACTIVE_OUTPUTS()

    SETUP_MODULE5_TOOLS()

    PREPARE_REFINED_MAG_REFERENCE(
        refined_manifest_ch
    )

    MAP_READS_TO_REFINED_MAGS(
        reads_ch
            .combine(PREPARE_REFINED_MAG_REFERENCE.out.reference_info)
            .combine(SETUP_MODULE5_TOOLS.out.status)
            .combine(CLEAN_MODULE5_DURABLE_SUBTRACTIVE_OUTPUTS.out.status)
    )

    WRITE_SUBTRACTIVE_MAPPING_SUMMARY(
        MAP_READS_TO_REFINED_MAGS.out.mapping_stats.collect()
    )

    def assembly_jobs_ch = MAP_READS_TO_REFINED_MAGS.out.unmapped_reads.flatMap { sample_id, safe_id, assembly_sample_id, _layout, unmapped_interleaved, unmapped_records, unmapped_pairs_or_fragments ->

        assembler_list.collect { assembler_name ->
            tuple(
                sample_id,
                safe_id,
                assembly_sample_id,
                unmapped_interleaved,
                unmapped_records,
                unmapped_pairs_or_fragments,
                assembler_name,
            )
        }
    }

    ASSEMBLE_SUBTRACTIVE(
        assembly_jobs_ch.combine(SETUP_MODULE5_TOOLS.out.status)
    )

    WRITE_SUBTRACTIVE_ASSEMBLY_SUMMARIES(
        ASSEMBLE_SUBTRACTIVE.out.manifest_record.collect(),
        ASSEMBLE_SUBTRACTIVE.out.stats_file.collect(),
    )

    WRITE_SUBTRACTIVE_MODULE3_INPUTS(
        MAP_READS_TO_REFINED_MAGS.out.trimmed_manifest_record.collect(),
        ASSEMBLE_SUBTRACTIVE.out.module3_manifest_record.collect(),
    )

    if (do_secondpass) {
        RUN_SECOND_PASS_BINNING(
            WRITE_SUBTRACTIVE_MODULE3_INPUTS.out.trimmed_manifest,
            WRITE_SUBTRACTIVE_MODULE3_INPUTS.out.assembly_manifest,
        )

        COMBINE_ORIGINAL_AND_SUBTRACTIVE_BINNING_MANIFESTS(
            original_binning_manifest_ch,
            RUN_SECOND_PASS_BINNING.out.status,
        )

        RUN_FINAL_JOINT_REFINEMENT(
            COMBINE_ORIGINAL_AND_SUBTRACTIVE_BINNING_MANIFESTS.out.combined_manifest,
            COMBINE_ORIGINAL_AND_SUBTRACTIVE_BINNING_MANIFESTS.out.stats,
        )

        BUILD_FINAL_MAG_DATABASE_FROM_JOINT_REFINEMENT(
            RUN_FINAL_JOINT_REFINEMENT.out.status
        )
    } else {
        WRITE_SECOND_PASS_DISABLED_STATUS()

        BUILD_FINAL_MAG_DATABASE_FROM_JOINT_REFINEMENT(
            WRITE_SECOND_PASS_DISABLED_STATUS.out.status
        )
    }
}


process CLEAN_MODULE5_DURABLE_SUBTRACTIVE_OUTPUTS {
    tag "clean_module5_durable_subtractive_outputs"
    cache false

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "module5_durable_output_cleanup_status.tsv"

    output:
    path "module5_durable_output_cleanup_status.tsv", emit: status

    script:
    """
    set -euo pipefail

    printf 'step\tstatus\tmessage\n' > module5_durable_output_cleanup_status.tsv

    for output_dir in \
        "${params.outdir}/assemblies" \
        "${params.outdir}/unmapped_reads"; do
        if [[ -e "\$output_dir" ]]; then
            rm -rf "\$output_dir"
        fi
        mkdir -p "\$output_dir"
    done

    printf 'module5_durable_subtractive_output_cleanup\tcompleted\tRemoved prior durable subtractive assemblies and unmapped reads\n' >> module5_durable_output_cleanup_status.tsv
    """
}

process SETUP_MODULE5_TOOLS {
    tag "setup_subtractive_assembly_tools"
    cache false

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "module5_tools_status.env"

    output:
    path "module5_tools_status.env", emit: status

    script:
    def env_dir = params.tool_env_dir ?: "${params.outdir}/conda_envs/module5_tools"

    def want_megahit = params.megahit.toString().toBoolean()
    def want_metaspades = params.metaspades.toString().toBoolean()

    def packages = []
    packages << "python"
    packages << "openjdk=17.*"
    packages << "bbmap=${params.bbmap_version}"

    if (want_megahit) {
        packages << "megahit=${params.megahit_version}"
    }

    if (want_metaspades) {
        def spades_version_for_conda = params.spades_version.toString().replaceFirst(/-\d+$/, '')
        packages << "spades=${spades_version_for_conda}"
    }

    def package_string = packages
        .collect { pkg ->
            "\"${pkg}\""
        }
        .join(" \\\n        ")

    """
    set -euo pipefail

    STATUS_FILE="module5_tools_status.env"
    TOOL_ENV="${env_dir}"

    WANT_MEGAHIT="${want_megahit}"
    WANT_METASPADES="${want_metaspades}"

    echo "Module 5 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Requested packages: ${packages.join(' ')}" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    find_installer() {
        if command -v mamba >/dev/null 2>&1; then
            echo "mamba"
        elif command -v conda >/dev/null 2>&1; then
            echo "conda"
        else
            echo ""
        fi
    }

    check_env() {
        local prefix="\$1"
        local ok="true"

        echo "Checking Module 5 environment: \$prefix" >> "\$STATUS_FILE"

        if [[ ! -x "\$prefix/bin/bbmap.sh" ]]; then
            echo "Missing bbmap.sh" >> "\$STATUS_FILE"
            ok="false"
        fi

        if [[ ! -x "\$prefix/bin/java" && -x "\$prefix/lib/jvm/bin/java" ]]; then
            mkdir -p "\$prefix/bin"
            ln -sfn "\$prefix/lib/jvm/bin/java" "\$prefix/bin/java"
        fi

        if [[ ! -x "\$prefix/bin/java" ]]; then
            echo "Missing java" >> "\$STATUS_FILE"
            ok="false"
        fi

        if [[ "\$WANT_MEGAHIT" == "true" && ! -x "\$prefix/bin/megahit" ]]; then
            echo "Missing megahit" >> "\$STATUS_FILE"
            ok="false"
        fi

        if [[ "\$WANT_METASPADES" == "true" && ! -x "\$prefix/bin/metaspades.py" ]]; then
            echo "Missing metaspades.py" >> "\$STATUS_FILE"
            ok="false"
        fi

        [[ "\$ok" == "true" ]]
    }

    if [[ -d "\$TOOL_ENV" ]]; then
        if check_env "\$TOOL_ENV"; then
            echo "Existing Module 5 environment passed checks." >> "\$STATUS_FILE"
            echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
            echo "Module 5 tool setup finished: \$(date)" >> "\$STATUS_FILE"
            exit 0
        else
            echo "Existing Module 5 environment failed checks. Removing." >> "\$STATUS_FILE"
            rm -rf "\$TOOL_ENV"
        fi
    fi

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: Module 5 environment missing/broken and --auto_install false." >> "\$STATUS_FILE"
        exit 1
    fi

    INSTALLER="\$(find_installer)"

    if [[ -z "\$INSTALLER" ]]; then
        echo "ERROR: Neither mamba nor conda found in PATH." >> "\$STATUS_FILE"
        exit 1
    fi

    mkdir -p "\$(dirname "\$TOOL_ENV")"

    "\$INSTALLER" create -y \\
        -p "\$TOOL_ENV" \\
        -c conda-forge \\
        -c bioconda \\
        ${package_string} \\
        >> "\$STATUS_FILE" 2>&1

    if ! check_env "\$TOOL_ENV"; then
        echo "ERROR: Newly created Module 5 environment failed checks." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Module 5 tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}


process PREPARE_REFINED_MAG_REFERENCE {
    tag "prepare_refined_mag_reference"

    publishDir "${params.outdir}/reference", mode: params.publish_reference_mode, pattern: "refined_mags_reference.fa"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "refined_mags_reference_stats.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "prepare_refined_mag_reference.log"

    input:
    path refined_manifest

    output:
    tuple path("refined_mags_reference.fa"), path("refined_mags_reference_stats.tsv"), emit: reference_info
    path "prepare_refined_mag_reference.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG_FILE="prepare_refined_mag_reference.log"
    REF_OUT="refined_mags_reference.fa"
    STATS_OUT="refined_mags_reference_stats.tsv"
    FASTA_LIST="refined_fasta_files.txt"

    echo "Preparing refined MAG reference: \$(date)" > "\$LOG_FILE"
    echo "Refined MAG manifest: ${refined_manifest}" >> "\$LOG_FILE"
    echo "Task working directory: \$(pwd -P)" >> "\$LOG_FILE"
    echo "Reference output file: \$REF_OUT" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    if [[ ! -s "${refined_manifest}" ]]; then
        echo "ERROR: refined MAG manifest is missing or empty: ${refined_manifest}" >> "\$LOG_FILE"
        exit 1
    fi

    python3 - "${refined_manifest}" "\$FASTA_LIST" <<'PY'
import csv
import sys
from pathlib import Path

manifest_path = Path(sys.argv[1]).resolve()
output_list = Path(sys.argv[2])

if not manifest_path.exists() or manifest_path.stat().st_size == 0:
    raise SystemExit(f"ERROR: manifest missing or empty: {manifest_path}")

selected = []
seen = set()

with manifest_path.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    if not reader.fieldnames:
        raise SystemExit(f"ERROR: manifest has no header: {manifest_path}")

    if "refined_bin_fasta" not in reader.fieldnames:
        raise SystemExit(
            "ERROR: refined MAG manifest does not contain required column "
            "'refined_bin_fasta'. Found columns: {}".format(
                ", ".join(reader.fieldnames)
            )
        )

    for row_number, row in enumerate(reader, start=2):
        value = (row.get("refined_bin_fasta") or "").strip()

        if not value:
            continue

        fasta = Path(value)

        if not fasta.is_absolute():
            fasta = manifest_path.parent / fasta

        fasta = fasta.resolve()

        if not fasta.exists():
            raise SystemExit(
                f"ERROR: refined MAG FASTA from manifest does not exist "
                f"(row {row_number}): {fasta}"
            )

        if fasta.stat().st_size == 0:
            raise SystemExit(
                f"ERROR: refined MAG FASTA from manifest is empty "
                f"(row {row_number}): {fasta}"
            )

        if fasta not in seen:
            selected.append(fasta)
            seen.add(fasta)

if not selected:
    raise SystemExit(
        f"ERROR: no non-empty refined_bin_fasta values found in {manifest_path}"
    )

with output_list.open("w") as out:
    for fasta in selected:
        print(fasta, file=out)
PY

    BINS_FOUND="\$(wc -l < "\$FASTA_LIST" | tr -d ' ')"

    echo "Refined MAG FASTAs selected from manifest: \$BINS_FOUND" >> "\$LOG_FILE"
    echo "Selected refined MAG FASTAs:" >> "\$LOG_FILE"
    cat "\$FASTA_LIST" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    : > "\$REF_OUT"

    while IFS= read -r fasta; do
        echo "Adding refined MAG FASTA to reference: \$fasta" >> "\$LOG_FILE"

        case "\$fasta" in
            *.gz)
                gzip -cd "\$fasta" >> "\$REF_OUT"
                ;;
            *)
                cat "\$fasta" >> "\$REF_OUT"
                ;;
        esac

        echo >> "\$REF_OUT"
    done < "\$FASTA_LIST"

    if [[ ! -s "\$REF_OUT" ]]; then
        echo "ERROR: reference FASTA was created but is empty: \$REF_OUT" >> "\$LOG_FILE"
        exit 1
    fi

    CONTIGS="\$(grep -c '^>' "\$REF_OUT" || true)"
    BP="\$(grep -v '^>' "\$REF_OUT" | tr -d '\\n[:space:]' | wc -c | tr -d ' ')"

    if [[ "\$CONTIGS" -eq 0 ]]; then
        echo "ERROR: reference FASTA contains zero contigs." >> "\$LOG_FILE"
        exit 1
    fi

    printf 'refined_manifest\\trefined_bins_used\\treference_contigs\\treference_bp\\treference_fasta\\n' > "\$STATS_OUT"

    printf '%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "\$(readlink -f "${refined_manifest}")" \\
        "\$BINS_FOUND" \\
        "\$CONTIGS" \\
        "\$BP" \\
        "\$REF_OUT" \\
        >> "\$STATS_OUT"

    echo "Reference contigs: \$CONTIGS" >> "\$LOG_FILE"
    echo "Reference bp: \$BP" >> "\$LOG_FILE"
    echo "Reference preparation finished: \$(date)" >> "\$LOG_FILE"
    """
}

process MAP_READS_TO_REFINED_MAGS {
    tag { sample_id }

    stageInMode 'symlink'

    publishDir "${params.outdir}/summary/per_sample_mapping", mode: 'copy', pattern: "*.subtractive_mapping_stats.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.bbmap_subtractive.log"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.mapping_threads as int
    }

    input:
    tuple val(sample_id), val(safe_id), val(assembly_sample_id), val(layout), path(read_files), path(refined_reference_fasta), path(reference_stats), path(tools_status), path(cleanup_status)

    output:
    tuple val(sample_id), val(safe_id), val(assembly_sample_id), val(layout), path("${safe_id}.unmapped_interleaved.fastq.gz"), val(0), val(0), emit: unmapped_reads
    path "${safe_id}.subtractive_mapping_stats.tsv", emit: mapping_stats
    path "${safe_id}.subtractive_trimmed_manifest_record.tsv", emit: trimmed_manifest_record
    path "${safe_id}.bbmap_subtractive.log", emit: log_file

    script:
    def durable_module5_outdir = absOrEmpty(params.outdir)
    def staged_read_files = read_files instanceof Collection ? read_files.toList() : [read_files]
    def staged_read_count = staged_read_files.size()

    if ((layout == 'paired' && staged_read_count != 2) ||
        (layout == 'interleaved' && staged_read_count != 1)) {
        error("Layout '${layout}' for sample '${sample_id}' requires ${layout == 'paired' ? 2 : 1} read file(s), but Nextflow staged ${staged_read_count}.")
    }

    def bbmap_xmx = params.bbmap_xmx == null ? '' : params.bbmap_xmx.toString().trim()
    if (bbmap_xmx.equalsIgnoreCase('null') || bbmap_xmx.equalsIgnoreCase('na')) {
        bbmap_xmx = ''
    }
    if (bbmap_xmx.startsWith('-Xmx')) {
        bbmap_xmx = bbmap_xmx.substring(4)
    }
    if (bbmap_xmx && !(bbmap_xmx ==~ /(?i)^\d+[kmgt]?$/)) {
        error("Invalid --bbmap_xmx value '${params.bbmap_xmx}'. Use a Java heap size such as '32g' or '8000m' (without -Xmx).")
    }

    def bbmap_xmx_arg = bbmap_xmx ? " -Xmx${bbmap_xmx}" : ''
    def bbmap_extra_args = params.bbmap_extra_args == null ? '' : params.bbmap_extra_args.toString().trim()
    if (bbmap_extra_args.equalsIgnoreCase('null') || bbmap_extra_args.equalsIgnoreCase('na')) {
        bbmap_extra_args = ''
    }

    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" && "\$TOOL_ENV" != "NOT_USED" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    TASK_DIR="\$(pwd -P)"

    LOG_FILE="${safe_id}.bbmap_subtractive.log"
    OUT_FASTQ="${safe_id}.unmapped_interleaved.fastq.gz"
    STATS_FILE="${safe_id}.subtractive_mapping_stats.tsv"
    RECORD_FILE="${safe_id}.subtractive_trimmed_manifest_record.tsv"

    REF_FASTA="${refined_reference_fasta}"

    echo "Subtractive BBMap mapping started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: ${sample_id}" >> "\$LOG_FILE"
    echo "Safe sample ID: ${safe_id}" >> "\$LOG_FILE"
    echo "Assembly sample ID: ${assembly_sample_id}" >> "\$LOG_FILE"
    echo "Layout: ${layout}" >> "\$LOG_FILE"
    echo "Reference FASTA staged in task: \$REF_FASTA" >> "\$LOG_FILE"
    echo "Reference stats file: ${reference_stats}" >> "\$LOG_FILE"
    echo "Task working directory: \$TASK_DIR" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    if [[ ! -s "\$REF_FASTA" ]]; then
        echo "ERROR: refined MAG reference FASTA is missing or empty: \$REF_FASTA" >> "\$LOG_FILE"
        exit 1
    fi

    REF_CONTIGS="\$(grep -c '^>' "\$REF_FASTA" || true)"

    echo "Reference contigs: \$REF_CONTIGS" >> "\$LOG_FILE"

    BBMAP_EXIT_STATUS=0
    MAPPING_STATUS="completed"

    if [[ "\$REF_CONTIGS" -eq 0 ]]; then
        echo "ERROR: refined MAG reference contains zero contigs: \$REF_FASTA" >> "\$LOG_FILE"
        echo "This should not happen if PREPARE_REFINED_MAG_REFERENCE succeeded." >> "\$LOG_FILE"
        exit 1
    fi

    if ! command -v bbmap.sh >/dev/null 2>&1; then
        echo "ERROR: bbmap.sh not available." >> "\$LOG_FILE"
        exit 1
    fi

    echo "BBMap executable: \$(command -v bbmap.sh)" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    BBMAP_ARGS="ref=\$REF_FASTA outu=\$OUT_FASTQ threads=${task.cpus} overwrite=t minid=${params.bbmap_minid} ambig=${params.bbmap_ambig}${bbmap_xmx_arg}"

    if [[ -n "${bbmap_extra_args}" ]]; then
        BBMAP_ARGS="\${BBMAP_ARGS} ${bbmap_extra_args}"
    fi

    echo "BBMap args: \$BBMAP_ARGS" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    if [[ "${layout}" == "paired" ]]; then
        ln -sfn "${staged_read_files[0]}" input_R1.fastq.gz
        ln -sfn "${staged_read_files[1]}" input_R2.fastq.gz

        set +e
        bbmap.sh \\
            \$BBMAP_ARGS \\
            in1=input_R1.fastq.gz \\
            in2=input_R2.fastq.gz \\
            interleaved=t \\
            >> "\$LOG_FILE" 2>&1
        BBMAP_EXIT_STATUS="\$?"
        set -e

    elif [[ "${layout}" == "interleaved" ]]; then
        ln -sfn "${staged_read_files[0]}" input_interleaved.fastq.gz

        set +e
        bbmap.sh \\
            \$BBMAP_ARGS \\
            in=input_interleaved.fastq.gz \\
            interleaved=t \\
            >> "\$LOG_FILE" 2>&1
        BBMAP_EXIT_STATUS="\$?"
        set -e

    else
        echo "ERROR: unsupported layout: ${layout}" >> "\$LOG_FILE"
        exit 1
    fi

    if [[ "\$BBMAP_EXIT_STATUS" -ne 0 ]]; then
        echo "ERROR: bbmap.sh failed with status \$BBMAP_EXIT_STATUS" >> "\$LOG_FILE"
        exit "\$BBMAP_EXIT_STATUS"
    fi

    if [[ ! -s "\$OUT_FASTQ" ]]; then
        echo "BBMap produced no unmapped reads. Creating empty gzipped FASTQ." >> "\$LOG_FILE"
        gzip -c /dev/null > "\$OUT_FASTQ"
    fi

    DURABLE_UNMAPPED_DIR="${durable_module5_outdir}/unmapped_reads"
    DURABLE_OUT_FASTQ="\$DURABLE_UNMAPPED_DIR/\$OUT_FASTQ"

    mkdir -p "\$DURABLE_UNMAPPED_DIR"
    cp -f "\$OUT_FASTQ" "\$DURABLE_OUT_FASTQ"

    if [[ ! -s "\$DURABLE_OUT_FASTQ" ]]; then
        echo "ERROR: failed to create durable unmapped FASTQ: \$DURABLE_OUT_FASTQ" >> "\$LOG_FILE"
        exit 1
    fi

    python3 - \\
        "\$OUT_FASTQ" \\
        "\$STATS_FILE" \\
        "${sample_id}" \\
        "${safe_id}" \\
        "${assembly_sample_id}" \\
        "${layout}" \\
        "\$REF_CONTIGS" \\
        "\$MAPPING_STATUS" \\
        "\$BBMAP_EXIT_STATUS" \\
        "\$DURABLE_OUT_FASTQ" <<'PY'
import gzip
import sys

(
    fastq,
    stats,
    sample_id,
    safe_id,
    assembly_sample_id,
    layout,
    ref_contigs,
    mapping_status,
    bbmap_exit_status,
    published,
) = sys.argv[1:]

lines = 0

with gzip.open(fastq, "rt", errors="replace") as handle:
    for _ in handle:
        lines += 1

if lines % 4 != 0:
    raise RuntimeError("FASTQ line count not divisible by 4: {}".format(lines))

records = lines // 4
pairs_or_fragments = records // 2

with open(stats, "w") as out:
    print(
        "sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "layout",
        "refined_reference_contigs",
        "unmapped_fastq_records",
        "unmapped_pairs_or_fragments",
        "mapping_status",
        "bbmap_exit_status",
        "unmapped_fastq",
        sep="\t",
        file=out,
    )

    print(
        sample_id,
        safe_id,
        assembly_sample_id,
        layout,
        ref_contigs,
        records,
        pairs_or_fragments,
        mapping_status,
        bbmap_exit_status,
        published,
        sep="\t",
        file=out,
    )
PY

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
        "${sample_id}" \\
        "${safe_id}" \\
        "interleaved" \\
        "" \\
        "" \\
        "\$DURABLE_OUT_FASTQ" \\
        "" \\
        "" \\
        "" \\
        > "\$RECORD_FILE"
    """
}

process ASSEMBLE_SUBTRACTIVE {
    tag { "${sample_id}:${assembler}" }

    stageInMode 'symlink'

    publishDir "${params.outdir}/header_maps", mode: 'copy', pattern: "*_subtractive.header_map.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*_subtractive_assembly.log"
    publishDir "${params.outdir}/summary/per_assembly_stats", mode: 'copy', pattern: "*_subtractive_assembly_stats.tsv"

    cpus {
        if (assembler == 'megahit' && params.megahit_threads != null) {
            return params.megahit_threads as int
        }

        return params.threads != null
            ? params.threads as int
            : params.assembly_threads as int
    }

    input:
    tuple val(sample_id), val(safe_id), val(assembly_sample_id), path(unmapped_interleaved), val(_unmapped_records_placeholder), val(_unmapped_pairs_placeholder), val(assembler), path(tools_status)

    output:
    path "*_subtractive.renamed.fa", emit: renamed_contigs
    path "*_subtractive.header_map.tsv", emit: header_map
    path "*_subtractive_assembly_stats.tsv", emit: stats_file
    path "*_subtractive_assembly_manifest_record.tsv", emit: manifest_record
    path "*_subtractive_module3_assembly_manifest_record.tsv", emit: module3_manifest_record
    path "*_subtractive_assembly.log", emit: log_file

    script:
    def durable_module5_outdir = absOrEmpty(params.outdir)
    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"
    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" && "\$TOOL_ENV" != "NOT_USED" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    SAMPLE_ID="${sample_id}"
    SAFE_ID="${safe_id}"
    ASSEMBLY_SAMPLE_ID="${assembly_sample_id}"
    ASSEMBLER="${assembler}"
    UNMAPPED_FASTQ="${unmapped_interleaved}"

    case "\$ASSEMBLER" in
        megahit)
            STRATEGY="E"
            RAW_OUT="megahit_subtractive_out"
            ;;
        metaspades)
            STRATEGY="F"
            RAW_OUT="metaspades_subtractive_out"
            ;;
        *)
            echo "ERROR: Unsupported subtractive assembler: \$ASSEMBLER" >&2
            exit 1
            ;;
    esac

    OUT_FASTA="\${SAFE_ID}_\${ASSEMBLER}_subtractive.renamed.fa"
    HEADER_MAP="\${SAFE_ID}_\${ASSEMBLER}_subtractive.header_map.tsv"
    STATS_FILE="\${SAFE_ID}_\${ASSEMBLER}_subtractive_assembly_stats.tsv"
    MANIFEST_RECORD="\${SAFE_ID}_\${ASSEMBLER}_subtractive_assembly_manifest_record.tsv"
    MODULE3_MANIFEST_RECORD="\${SAFE_ID}_\${ASSEMBLER}_subtractive_module3_assembly_manifest_record.tsv"
    LOG_FILE="\${SAFE_ID}_\${ASSEMBLER}_subtractive_assembly.log"

    echo "Subtractive assembly started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: \$SAMPLE_ID" >> "\$LOG_FILE"
    echo "Assembler: \$ASSEMBLER" >> "\$LOG_FILE"
    echo "Assembly strategy: \$STRATEGY" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"

    : > "\$OUT_FASTA"
    printf 'old_header\\tnew_header\\n' > "\$HEADER_MAP"
    : > "\$MANIFEST_RECORD"
    : > "\$MODULE3_MANIFEST_RECORD"

    UNMAPPED_RECORDS=0
    if [[ -s "\$UNMAPPED_FASTQ" ]]; then
        gzip -t "\$UNMAPPED_FASTQ" >> "\$LOG_FILE" 2>&1
        UNMAPPED_RECORDS="\$(gzip -cd "\$UNMAPPED_FASTQ" | awk 'END { print NR / 4 }')"
    fi
    UNMAPPED_PAIRS="\$((UNMAPPED_RECORDS / 2))"

    SRC_FASTA=""
    ASSEMBLER_EXIT_STATUS=0
    ASSEMBLY_STATUS="completed"
    ASSEMBLY_MESSAGE="Assembly completed"

    if [[ "\$UNMAPPED_RECORDS" -eq 0 || "\$UNMAPPED_PAIRS" -eq 0 ]]; then
        ASSEMBLY_STATUS="skipped_no_unmapped_reads"
        ASSEMBLY_MESSAGE="No unmapped reads available"
    elif [[ "\$ASSEMBLER" == "megahit" ]]; then
        if ! command -v megahit >/dev/null 2>&1; then
            ASSEMBLER_EXIT_STATUS=127
            ASSEMBLY_STATUS="tool_missing"
            ASSEMBLY_MESSAGE="megahit not available"
        else
            set +e
            megahit --12 "\$UNMAPPED_FASTQ" -t ${task.cpus} -o "\$RAW_OUT" --presets ${params.megahit_preset} >> "\$LOG_FILE" 2>&1
            ASSEMBLER_EXIT_STATUS="\$?"
            set -e
            if [[ "\$ASSEMBLER_EXIT_STATUS" -ne 0 ]]; then
                ASSEMBLY_STATUS="failed_nonfatal"
                ASSEMBLY_MESSAGE="MEGAHIT exited non-zero"
            elif [[ -s "\$RAW_OUT/final.contigs.fa" ]]; then
                SRC_FASTA="\$RAW_OUT/final.contigs.fa"
                ASSEMBLY_MESSAGE="MEGAHIT completed and produced contigs"
            else
                ASSEMBLY_STATUS="no_contigs"
                ASSEMBLY_MESSAGE="MEGAHIT produced no contigs"
            fi
        fi
    else
        if ! command -v metaspades.py >/dev/null 2>&1; then
            ASSEMBLER_EXIT_STATUS=127
            ASSEMBLY_STATUS="tool_missing"
            ASSEMBLY_MESSAGE="metaspades.py not available"
        else
            METASPADES_MEMORY_ARG=""
            if [[ "${params.metaspades_memory_gb}" -gt 0 ]]; then
                METASPADES_MEMORY_ARG="-m ${params.metaspades_memory_gb}"
            fi
            set +e
            metaspades.py --12 "\$UNMAPPED_FASTQ" -t ${task.cpus} \$METASPADES_MEMORY_ARG -o "\$RAW_OUT" >> "\$LOG_FILE" 2>&1
            ASSEMBLER_EXIT_STATUS="\$?"
            set -e
            if [[ "\$ASSEMBLER_EXIT_STATUS" -ne 0 ]]; then
                ASSEMBLY_STATUS="failed_nonfatal"
                ASSEMBLY_MESSAGE="metaSPAdes exited non-zero"
            elif [[ -s "\$RAW_OUT/scaffolds.fasta" ]]; then
                SRC_FASTA="\$RAW_OUT/scaffolds.fasta"
                ASSEMBLY_MESSAGE="metaSPAdes completed and produced scaffolds"
            elif [[ -s "\$RAW_OUT/contigs.fasta" ]]; then
                SRC_FASTA="\$RAW_OUT/contigs.fasta"
                ASSEMBLY_MESSAGE="metaSPAdes completed and produced contigs"
            else
                ASSEMBLY_STATUS="no_contigs"
                ASSEMBLY_MESSAGE="metaSPAdes produced no contigs"
            fi
        fi
    fi

    if [[ -n "\$SRC_FASTA" && -s "\$SRC_FASTA" ]]; then
        python3 - "\$SRC_FASTA" "\$OUT_FASTA" "\$HEADER_MAP" "\$ASSEMBLY_SAMPLE_ID" "\$STRATEGY" <<'PY2'
import re
import sys

src, out_fa, hmap, sample, strategy = sys.argv[1:]
with open(src) as inp, open(out_fa, 'w') as out, open(hmap, 'w') as hm:
    print('old_header', 'new_header', sep='\t', file=hm)
    for idx, line in enumerate(inp):
        if line.startswith('>'):
            old = line[1:].strip()
            token = old.split()[0] if old else f'contig_{idx}'
            suffix = re.sub(r'[^A-Za-z0-9._-]+', '_', token)
            new = f'{sample}_{strategy}_{suffix}'
            print(f'>{new}', file=out)
            print(old, new, sep='\t', file=hm)
        else:
            print(line.strip(), file=out)
PY2
    fi

    DURABLE_ASSEMBLIES_DIR="${durable_module5_outdir}/assemblies"
    DURABLE_OUT_FASTA="\$DURABLE_ASSEMBLIES_DIR/\$OUT_FASTA"
    mkdir -p "\$DURABLE_ASSEMBLIES_DIR"
    cp -f "\$OUT_FASTA" "\$DURABLE_OUT_FASTA"

    python3 - "\$OUT_FASTA" "\$STATS_FILE" "\$MANIFEST_RECORD" "\$MODULE3_MANIFEST_RECORD" "\$SAMPLE_ID" "\$SAFE_ID" "\$ASSEMBLY_SAMPLE_ID" "\$ASSEMBLER" "\$STRATEGY" "\$UNMAPPED_RECORDS" "\$UNMAPPED_PAIRS" "\$ASSEMBLER_EXIT_STATUS" "\$ASSEMBLY_STATUS" "\$ASSEMBLY_MESSAGE" "\$DURABLE_OUT_FASTA" <<'PY2'
import sys
from pathlib import Path

(fasta, stats, manifest, module3_manifest, sample, safe, assembly_sample, assembler, strategy,
 unmapped_records, unmapped_pairs, exit_status, status, message, published) = sys.argv[1:]
lengths, current, seen = [], 0, False
for line in Path(fasta).open():
    if line.startswith('>'):
        if seen: lengths.append(current)
        current, seen = 0, True
    else:
        current += len(line.strip())
if seen: lengths.append(current)
ordered = sorted(lengths, reverse=True)
half, running, n50 = sum(ordered) / 2, 0, 0
for length in ordered:
    running += length
    if running >= half:
        n50 = length
        break
with open(stats, 'w') as out:
    print('sample_id', 'safe_sample_id', 'assembly_sample_id', 'assembler', 'assembly_mode', 'unmapped_fastq_records', 'unmapped_pairs_or_fragments', 'contigs', 'total_bp', 'max_contig_bp', 'n50_bp', 'assembler_exit_status', 'assembly_status', 'assembly_message', 'renamed_fasta', sep='\t', file=out)
    print(sample, safe, assembly_sample, assembler, 'subtractive', unmapped_records, unmapped_pairs, len(lengths), sum(lengths), max(lengths) if lengths else 0, n50, exit_status, status, message, published, sep='\t', file=out)
with open(manifest, 'w') as out:
    print(sample, safe, assembly_sample, assembler, 'subtractive', published, sep='\t', file=out)
with open(module3_manifest, 'w') as out:
    if lengths:
        print(sample, safe, assembly_sample, assembler, 'subtractive', '', strategy, published, sep='\t', file=out)
PY2

    echo "Final assembly status: \$ASSEMBLY_STATUS" >> "\$LOG_FILE"
    echo "Output contigs: \$(grep -c '^>' "\$OUT_FASTA" || true)" >> "\$LOG_FILE"
    rm -rf "\$RAW_OUT"
    echo "Subtractive assembly finished: \$(date)" >> "\$LOG_FILE"
    """
}

process WRITE_SUBTRACTIVE_MAPPING_SUMMARY {
    tag "write_subtractive_mapping_summary"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "subtractive_mapping_summary.tsv"

    input:
    path stats_files

    output:
    path "subtractive_mapping_summary.tsv", emit: summary

    script:
    def stats_file_list = stats_files
        .collect { stats_file ->
            stats_file.name
        }
        .join(' ')

    """
    set -euo pipefail

    first=1
    : > subtractive_mapping_summary.tsv

    for f in ${stats_file_list}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" >> subtractive_mapping_summary.tsv
            first=0
        else
            tail -n +2 "\$f" >> subtractive_mapping_summary.tsv
        fi
    done
    """
}


process WRITE_SUBTRACTIVE_ASSEMBLY_SUMMARIES {
    tag "write_subtractive_assembly_summaries"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "subtractive_assembly_summary_manifest.tsv"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "subtractive_assembly_stats_summary.tsv"

    input:
    path manifest_records
    path stats_files

    output:
    path "subtractive_assembly_summary_manifest.tsv", emit: manifest
    path "subtractive_assembly_stats_summary.tsv", emit: stats_summary

    script:
    def manifest_file_list = manifest_records
        .collect { manifest_file ->
            manifest_file.name
        }
        .join(' ')

    def stats_file_list = stats_files
        .collect { stats_file ->
            stats_file.name
        }
        .join(' ')

    """
    set -euo pipefail

    printf 'sample_id\\tsafe_sample_id\\tassembly_sample_id\\tassembler\\tassembly_mode\\trenamed_fasta\\n' > subtractive_assembly_summary_manifest.tsv

    for f in ${manifest_file_list}; do
        cat "\$f" >> subtractive_assembly_summary_manifest.tsv
    done

    first=1
    : > subtractive_assembly_stats_summary.tsv

    for f in ${stats_file_list}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" >> subtractive_assembly_stats_summary.tsv
            first=0
        else
            tail -n +2 "\$f" >> subtractive_assembly_stats_summary.tsv
        fi
    done
    """
}


process WRITE_SUBTRACTIVE_MODULE3_INPUTS {
    tag "write_subtractive_module3_inputs"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "subtractive_trimmed_manifest.tsv"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "subtractive_module3_assembly_manifest.tsv"

    input:
    path trimmed_manifest_records
    path assembly_manifest_records

    output:
    path "subtractive_trimmed_manifest.tsv", emit: trimmed_manifest
    path "subtractive_module3_assembly_manifest.tsv", emit: assembly_manifest

    script:
    def trimmed_files = trimmed_manifest_records
        .collect { record_file ->
            record_file.name
        }
        .join(' ')

    def assembly_files = assembly_manifest_records
        .collect { record_file ->
            record_file.name
        }
        .join(' ')

    """
    set -euo pipefail

    printf 'sample_id\\tsafe_sample_id\\tlayout\\tread1\\tread2\\tinterleaved\\tmerged\\tfastp_html\\tfastp_json\\n' > subtractive_trimmed_manifest.tsv

    if [[ -n "${trimmed_files}" ]]; then
        for f in ${trimmed_files}; do
            if [[ -s "\$f" ]]; then
                cat "\$f" >> subtractive_trimmed_manifest.tsv
            fi
        done
    fi

    printf 'sample_id\\tsafe_sample_id\\tassembly_sample_id\\tassembler\\tassembly_mode\\trarefaction_label\\tassembly_strategy\\trenamed_fasta\\n' > subtractive_module3_assembly_manifest.tsv

    if [[ -n "${assembly_files}" ]]; then
        for f in ${assembly_files}; do
            if [[ -s "\$f" ]]; then
                cat "\$f" >> subtractive_module3_assembly_manifest.tsv
            fi
        done
    fi
    """
}


process RUN_SECOND_PASS_BINNING {
    tag "run_second_pass_binning"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "second_pass_binning_status.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "second_pass_binning.log"

    input:
    path subtractive_trimmed_manifest
    path subtractive_assembly_manifest

    output:
    path "second_pass_binning_status.tsv", emit: status
    path "second_pass_binning.log", emit: log_file

    script:
    def binner_args = []

    if (params.secondpass_metabat2.toString().toBoolean()) {
        binner_args << "--metabat2"
    }

    if (params.secondpass_quickbin.toString().toBoolean()) {
        binner_args << "--quickbin"
    }

    if (params.secondpass_maxbin2.toString().toBoolean()) {
        binner_args << "--maxbin2"
    }

    def binner_arg_string = binner_args.join(' ')
    def thread_arg = params.threads != null ? "--threads ${params.threads}" : ""

    """
    set -euo pipefail

    LOG_FILE="second_pass_binning.log"
    SECOND_PASS_MODULE3_DIR="${params.secondpass_dir}/module_3_binning"
    SECOND_PASS_BINNING_MANIFEST="\$SECOND_PASS_MODULE3_DIR/summary/binning_manifest.tsv"

    echo "Second-pass binning started: \$(date)" > "\$LOG_FILE"
    echo "Subtractive trimmed manifest: ${subtractive_trimmed_manifest}" >> "\$LOG_FILE"
    echo "Subtractive assembly manifest: ${subtractive_assembly_manifest}" >> "\$LOG_FILE"
    echo "Module 3 script: ${params.module3_script}" >> "\$LOG_FILE"
    echo "Second-pass dir: ${params.secondpass_dir}" >> "\$LOG_FILE"
    echo "Binner args: ${binner_arg_string}" >> "\$LOG_FILE"

    printf 'step\\tstatus\\texit_status\\tmessage\\n' > second_pass_binning_status.tsv

    ASSEMBLY_ROWS="\$(tail -n +2 "${subtractive_assembly_manifest}" | awk 'NF > 0' | wc -l | tr -d ' ')"

    echo "Subtractive assembly rows: \$ASSEMBLY_ROWS" >> "\$LOG_FILE"

    if [[ "\$ASSEMBLY_ROWS" -eq 0 ]]; then
        echo "No non-empty subtractive assemblies were available." >> "\$LOG_FILE"
        echo "No new MAGs were added with subtractive assembly." >> "\$LOG_FILE"
        printf 'second_pass_binning\\tskipped_no_subtractive_assemblies\\t0\\tNo new MAGs were added with subtractive assembly\\n' >> second_pass_binning_status.tsv
        exit 0
    fi

    if [[ -z "${binner_arg_string}" ]]; then
        echo "No second-pass binners selected." >> "\$LOG_FILE"
        printf 'second_pass_binning\\tskipped_no_binners\\t0\\tNo second-pass binners selected\\n' >> second_pass_binning_status.tsv
        exit 0
    fi

    if [[ ! -s "${params.module3_script}" ]]; then
        echo "ERROR: Module 3 script missing: ${params.module3_script}" >> "\$LOG_FILE"
        printf 'second_pass_binning\\tfailed\\t1\\tModule 3 script missing\\n' >> second_pass_binning_status.tsv
        exit 1
    fi

    SUB_TRIMMED="\$(pwd -P)/${subtractive_trimmed_manifest}"
    SUB_ASSEMBLY="\$(pwd -P)/${subtractive_assembly_manifest}"

    # Module 3 writes its results below a persistent second-pass directory. Remove
    # the prior Module 3 result tree so its manifest cannot be mistaken for output
    # from this invocation.
    if [[ -e "\$SECOND_PASS_MODULE3_DIR" ]]; then
        echo "Removing previous Module 3 second-pass output: \$SECOND_PASS_MODULE3_DIR" >> "\$LOG_FILE"
        rm -rf "\$SECOND_PASS_MODULE3_DIR"
    fi

    set +e
    ${params.nextflow_exe} run "${params.module3_script}" \\
        --working_dir "${params.secondpass_dir}" \\
        --input_trimmed_manifest "\$SUB_TRIMMED" \\
        --input_assembly_manifest "\$SUB_ASSEMBLY" \\
        ${thread_arg} \\
        ${binner_arg_string} \\
        >> "\$LOG_FILE" 2>&1

    STATUS_CODE="\$?"
    set -e

    BIN_ROWS=0

    if [[ -s "\$SECOND_PASS_BINNING_MANIFEST" ]]; then
        BIN_ROWS="\$(tail -n +2 "\$SECOND_PASS_BINNING_MANIFEST" | awk 'NF > 0' | wc -l | tr -d ' ')"
    fi

    echo "Second-pass Module 3 exit status: \$STATUS_CODE" >> "\$LOG_FILE"
    echo "Second-pass bin rows: \$BIN_ROWS" >> "\$LOG_FILE"

    if [[ "\$STATUS_CODE" -ne 0 ]]; then
        if [[ "\$BIN_ROWS" -eq 0 ]]; then
            echo "WARNING: Module 3 exited non-zero but produced no bins." >> "\$LOG_FILE"
            echo "Treating as non-fatal because no new MAGs were added." >> "\$LOG_FILE"
            printf 'second_pass_binning\\tskipped_no_new_bins\\t0\\tNo new MAGs were added with subtractive assembly\\n' >> second_pass_binning_status.tsv
            exit 0
        else
            printf 'second_pass_binning\\tfailed\\t%s\\tModule 3 second-pass binning failed\\n' "\$STATUS_CODE" >> second_pass_binning_status.tsv
            exit "\$STATUS_CODE"
        fi
    fi

    if [[ "\$BIN_ROWS" -eq 0 ]]; then
        echo "No bins were produced by second-pass binning." >> "\$LOG_FILE"
        echo "No new MAGs were added with subtractive assembly." >> "\$LOG_FILE"
        printf 'second_pass_binning\\tskipped_no_new_bins\\t0\\tNo new MAGs were added with subtractive assembly\\n' >> second_pass_binning_status.tsv
        exit 0
    fi

    printf 'second_pass_binning\\tcompleted\\t0\\tModule 3 second-pass binning completed and produced new bins\\n' >> second_pass_binning_status.tsv

    echo "Second-pass binning finished: \$(date)" >> "\$LOG_FILE"
    """
}


process WRITE_SECOND_PASS_DISABLED_STATUS {
    tag "second_pass_disabled"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "final_joint_refinement_status.tsv"

    output:
    path "final_joint_refinement_status.tsv", emit: status

    script:
    """
    set -euo pipefail

    printf 'step\\tstatus\\texit_status\\tmessage\\n' > final_joint_refinement_status.tsv
    printf 'final_joint_refinement\\tskipped_second_pass_disabled\\t0\\tSecond-pass binning and final joint refinement disabled; using original Module 4 refined MAGs\\n' >> final_joint_refinement_status.tsv
    """
}


process COMBINE_ORIGINAL_AND_SUBTRACTIVE_BINNING_MANIFESTS {
    tag "combine_original_and_subtractive_binning_manifests"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "combined_original_plus_subtractive_binning_manifest.tsv"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "combined_original_plus_subtractive_binning_manifest_stats.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "combine_binning_manifests.log"

    input:
    path original_binning_manifest
    path second_pass_binning_status

    output:
    path "combined_original_plus_subtractive_binning_manifest.tsv", emit: combined_manifest
    path "combined_original_plus_subtractive_binning_manifest_stats.tsv", emit: stats
    path "combine_binning_manifests.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG_FILE="combine_binning_manifests.log"

    ORIGINAL="${original_binning_manifest}"
    SECOND_PASS_STATUS_FILE="${second_pass_binning_status}"
    SUBTRACTIVE="${params.secondpass_dir}/module_3_binning/summary/binning_manifest.tsv"
    OUT="combined_original_plus_subtractive_binning_manifest.tsv"
    STATS="combined_original_plus_subtractive_binning_manifest_stats.tsv"

    echo "Combining original and subtractive binning manifests: \$(date)" > "\$LOG_FILE"
    echo "Original manifest: \$ORIGINAL" >> "\$LOG_FILE"
    echo "Subtractive manifest: \$SUBTRACTIVE" >> "\$LOG_FILE"
    echo "Second-pass binning status file: \$SECOND_PASS_STATUS_FILE" >> "\$LOG_FILE"

    if [[ ! -s "\$ORIGINAL" ]]; then
        echo "ERROR: Original binning manifest is missing or empty: \$ORIGINAL" >> "\$LOG_FILE"
        exit 1
    fi

    if [[ ! -s "\$SECOND_PASS_STATUS_FILE" ]]; then
        echo "ERROR: second-pass binning status is missing or empty: \$SECOND_PASS_STATUS_FILE" >> "\$LOG_FILE"
        exit 1
    fi

    SECOND_PASS_STATUS="\$(tail -n 1 "\$SECOND_PASS_STATUS_FILE" | cut -f2 | tr -d '[:space:]')"

    if [[ -z "\$SECOND_PASS_STATUS" ]]; then
        echo "ERROR: could not determine second-pass binning status." >> "\$LOG_FILE"
        exit 1
    fi

    echo "Second-pass status: \$SECOND_PASS_STATUS" >> "\$LOG_FILE"

    HEADER="\$(head -n 1 "\$ORIGINAL")"
    printf '%s\\n' "\$HEADER" > "\$OUT"

    ORIGINAL_ROWS="\$(tail -n +2 "\$ORIGINAL" | awk 'NF > 0' | wc -l | tr -d ' ')"
    SUBTRACTIVE_ROWS=0

    tail -n +2 "\$ORIGINAL" | awk 'NF > 0' >> "\$OUT"

    if [[ "\$SECOND_PASS_STATUS" == "completed" ]]; then
        if [[ ! -s "\$SUBTRACTIVE" ]]; then
            echo "ERROR: second-pass status is completed, but subtractive binning manifest is missing or empty: \$SUBTRACTIVE" >> "\$LOG_FILE"
            exit 1
        fi

        SUBTRACTIVE_ROWS="\$(tail -n +2 "\$SUBTRACTIVE" | awk 'NF > 0' | wc -l | tr -d ' ')"
        tail -n +2 "\$SUBTRACTIVE" | awk 'NF > 0' >> "\$OUT"
    else
        echo "Second-pass status is \$SECOND_PASS_STATUS; not using the subtractive manifest. Final joint refinement will use original bins only." >> "\$LOG_FILE"
    fi

    TOTAL_ROWS="\$(tail -n +2 "\$OUT" | awk 'NF > 0' | wc -l | tr -d ' ')"

    printf 'original_rows\\tsubtractive_rows\\ttotal_rows\\toriginal_manifest\\tsubtractive_manifest\\tcombined_manifest\\n' > "\$STATS"
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' "\$ORIGINAL_ROWS" "\$SUBTRACTIVE_ROWS" "\$TOTAL_ROWS" "\$ORIGINAL" "\$SUBTRACTIVE" "\$OUT" >> "\$STATS"

    echo "Original rows: \$ORIGINAL_ROWS" >> "\$LOG_FILE"
    echo "Subtractive rows: \$SUBTRACTIVE_ROWS" >> "\$LOG_FILE"
    echo "Total combined rows: \$TOTAL_ROWS" >> "\$LOG_FILE"
    echo "Combining manifests finished: \$(date)" >> "\$LOG_FILE"
    """
}


process RUN_FINAL_JOINT_REFINEMENT {
    tag "run_final_joint_refinement"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "final_joint_refinement_status.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "final_joint_refinement.log"

    input:
    path combined_binning_manifest
    path combined_manifest_stats

    output:
    path "final_joint_refinement_status.tsv", emit: status
    path "final_joint_refinement.log", emit: log_file

    script:
    def thread_arg = params.threads != null ? "--threads ${params.threads}" : ""

    def module4_extra_args = "--dependencies_dir \"${params.dependencies_dir}\""

    if (params.tigrfam_hmm) {
        module4_extra_args += " --tigrfam_hmm \"${params.tigrfam_hmm}\""
    }

    if (params.pfam_hmm) {
        module4_extra_args += " --pfam_hmm \"${params.pfam_hmm}\""
    }

    if (params.magscot_script) {
        module4_extra_args += " --magscot_script \"${params.magscot_script}\""
    }

    if (params.magscot_extra_args) {
        module4_extra_args += " --magscot_extra_args \"${params.magscot_extra_args}\""
    }

    """
    set -euo pipefail

    LOG_FILE="final_joint_refinement.log"
    FINAL_JOINT_MODULE4_DIR="${params.final_joint_dir}/module_4_binrefinement"

    echo "Final joint refinement started: \$(date)" > "\$LOG_FILE"
    echo "Combined binning manifest: ${combined_binning_manifest}" >> "\$LOG_FILE"
    echo "Combined manifest stats: ${combined_manifest_stats}" >> "\$LOG_FILE"
    echo "Module 4 script: ${params.module4_script}" >> "\$LOG_FILE"
    echo "Final joint working directory: ${params.final_joint_dir}" >> "\$LOG_FILE"

    printf 'step\\tstatus\\texit_status\\tmessage\\n' > final_joint_refinement_status.tsv

    TOTAL_ROWS="\$(tail -n +2 "${combined_binning_manifest}" | awk 'NF > 0' | wc -l | tr -d ' ')"

    SUBTRACTIVE_ROWS=0
    if [[ -s "${combined_manifest_stats}" ]]; then
        SUBTRACTIVE_ROWS="\$(tail -n 1 "${combined_manifest_stats}" | cut -f2 | tr -d ' ')"
    fi

    echo "Total combined rows: \$TOTAL_ROWS" >> "\$LOG_FILE"
    echo "Subtractive rows: \$SUBTRACTIVE_ROWS" >> "\$LOG_FILE"

    if [[ "\$TOTAL_ROWS" -eq 0 ]]; then
        echo "No bins were available for final joint refinement." >> "\$LOG_FILE"
        printf 'final_joint_refinement\\tskipped_no_bins\\t0\\tNo bins available for final joint refinement\\n' >> final_joint_refinement_status.tsv
        exit 0
    fi

    if [[ "\$SUBTRACTIVE_ROWS" -eq 0 ]]; then
        echo "No new MAGs were added with subtractive assembly." >> "\$LOG_FILE"
        echo "Skipping final joint MAGScoT refinement." >> "\$LOG_FILE"
        printf 'final_joint_refinement\\tskipped_no_new_mags\\t0\\tNo new MAGs were added with subtractive assembly; final joint refinement skipped\\n' >> final_joint_refinement_status.tsv
        exit 0
    fi

    if [[ ! -s "${params.module4_script}" ]]; then
        printf 'final_joint_refinement\\tfailed\\t1\\tModule 4 script missing\\n' >> final_joint_refinement_status.tsv
        exit 1
    fi

    COMBINED="\$(pwd -P)/${combined_binning_manifest}"

    # Module 4 writes its results below a persistent final-joint directory. Remove
    # the prior Module 4 result tree so its manifest cannot be mistaken for output
    # from this invocation.
    if [[ -e "\$FINAL_JOINT_MODULE4_DIR" ]]; then
        echo "Removing previous Module 4 final-joint output: \$FINAL_JOINT_MODULE4_DIR" >> "\$LOG_FILE"
        rm -rf "\$FINAL_JOINT_MODULE4_DIR"
    fi

    set +e
    ${params.nextflow_exe} run "${params.module4_script}" \\
        --working_dir "${params.final_joint_dir}" \\
        --input_binning_manifest "\$COMBINED" \\
        ${thread_arg} \\
        ${module4_extra_args} \\
        >> "\$LOG_FILE" 2>&1

    STATUS_CODE="\$?"
    set -e

    if [[ "\$STATUS_CODE" -ne 0 ]]; then
        printf 'final_joint_refinement\\tfailed\\t%s\\tFinal joint Module 4 refinement failed\\n' "\$STATUS_CODE" >> final_joint_refinement_status.tsv
        exit "\$STATUS_CODE"
    fi

    printf 'final_joint_refinement\\tcompleted\\t0\\tFinal joint Module 4 refinement completed\\n' >> final_joint_refinement_status.tsv

    echo "Final joint refinement finished: \$(date)" >> "\$LOG_FILE"
    """
}


process BUILD_FINAL_MAG_DATABASE_FROM_JOINT_REFINEMENT {
    tag "build_final_mag_database_from_joint_refinement"

    publishDir "${params.outdir}/final_mag_database", mode: params.publish_final_mags_mode, pattern: "final_mag_database/*.fa", saveAs: { filename -> filename.replaceFirst(/^final_mag_database\//, '') }
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "final_mag_database_*.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "build_final_mag_database.log"

    input:
    path final_joint_refinement_status

    output:
    path "final_mag_database", emit: final_mags
    path "final_mag_database_manifest.tsv", emit: manifest
    path "final_mag_database_stats.tsv", emit: stats
    path "build_final_mag_database.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG_FILE="build_final_mag_database.log"
    PUBLISHED_FINAL_MAG_DIR="${params.outdir}/final_mag_database"

    echo "Final MAG database construction started: \$(date)" > "\$LOG_FILE"

    if [[ -e "\$PUBLISHED_FINAL_MAG_DIR" ]]; then
        echo "Removing previous published final MAG database: \$PUBLISHED_FINAL_MAG_DIR" >> "\$LOG_FILE"
        rm -rf "\$PUBLISHED_FINAL_MAG_DIR"
    fi

    FINAL_STATUS="\$(tail -n 1 "${final_joint_refinement_status}" | cut -f2 || true)"
    FINAL_MESSAGE="\$(tail -n 1 "${final_joint_refinement_status}" | cut -f4- || true)"

    FINAL_JOINT_MANIFEST="${params.final_joint_dir}/module_4_binrefinement/summary/magscot_refined_bins_manifest.tsv"
    ORIGINAL_REFINED_MANIFEST="${params.module4_outdir}/summary/magscot_refined_bins_manifest.tsv"

    echo "Final joint refinement status file: ${final_joint_refinement_status}" >> "\$LOG_FILE"
    echo "Final status: \$FINAL_STATUS" >> "\$LOG_FILE"
    echo "Final message: \$FINAL_MESSAGE" >> "\$LOG_FILE"
    echo "Final joint refined MAG manifest: \$FINAL_JOINT_MANIFEST" >> "\$LOG_FILE"
    echo "Original refined MAG manifest: \$ORIGINAL_REFINED_MANIFEST" >> "\$LOG_FILE"
    mkdir -p final_mag_database

    if [[ "\$FINAL_STATUS" == "completed" ]]; then
        SELECTED_MANIFEST="\$FINAL_JOINT_MANIFEST"
        FINAL_DATABASE_MODE="final_joint_refinement"
    elif [[ "\$FINAL_STATUS" == "skipped_no_new_mags" || "\$FINAL_STATUS" == "skipped_no_bins" || "\$FINAL_STATUS" == "skipped_second_pass_disabled" ]]; then
        echo "Using original Module 4 refined MAGs as final MAG database." >> "\$LOG_FILE"
        SELECTED_MANIFEST="\$ORIGINAL_REFINED_MANIFEST"
        FINAL_DATABASE_MODE="original_refined_mags_only"
    else
        echo "ERROR: Final joint refinement did not complete or skip cleanly." >> "\$LOG_FILE"
        echo "Status: \$FINAL_STATUS" >> "\$LOG_FILE"
        exit 1
    fi

    echo "Building final MAG database from selected manifest: \$SELECTED_MANIFEST" >> "\$LOG_FILE"

    python3 - \\
        "\$SELECTED_MANIFEST" \\
        "\$FINAL_DATABASE_MODE" \\
        "final_mag_database" \\
        "final_mag_database_manifest.tsv" \\
        "final_mag_database_stats.tsv" \\
        "\$LOG_FILE" \\
        "${params.outdir}/final_mag_database" <<'PY'
import csv
import gzip
import re
import shutil
import sys
from pathlib import Path

(
    selected_manifest,
    final_database_mode,
    final_dir,
    final_manifest,
    final_stats,
    log_file,
    published_final_dir
) = sys.argv[1:]

selected_manifest = Path(selected_manifest).resolve()
selected_manifest_dir = selected_manifest.parent
final_dir = Path(final_dir)
final_manifest = Path(final_manifest)
final_stats = Path(final_stats)
log_file = Path(log_file)
published_final_dir = Path(published_final_dir)

final_dir.mkdir(parents=True, exist_ok=True)

def log(message):
    with log_file.open("a") as h:
        print(message, file=h)

def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("_")
    return value or "unnamed_mag"

def fasta_stats(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    contigs = 0
    bp = 0
    cur = 0
    seen = False
    with opener(path, "rt", errors="replace") as h:
        for line in h:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if seen:
                    bp += cur
                cur = 0
                seen = True
                contigs += 1
            else:
                cur += len(line.strip())
    if seen:
        bp += cur
    return contigs, bp

def copy_as_plain_fasta(source, destination):
    source = Path(source)
    destination = Path(destination)

    if str(source).endswith(".gz"):
        with gzip.open(source, "rt", errors="replace") as inp, destination.open("w") as out:
            shutil.copyfileobj(inp, out)
    else:
        shutil.copyfile(source, destination)

rows = []
manifest_rows = 0

if not selected_manifest.exists() or selected_manifest.stat().st_size == 0:
    log(f"WARNING: selected refined MAG manifest missing or empty: {selected_manifest}")
else:
    with selected_manifest.open() as h:
        reader = csv.DictReader(h, delimiter="\\t")
        for row in reader:
            manifest_rows += 1
            fasta = row.get("refined_bin_fasta", "").strip()
            bin_id = row.get("refined_bin_id", "").strip()
            if fasta:
                rows.append((bin_id, fasta))

seen_names = set()
copied = 0
missing = 0
total_contigs = 0
total_bp = 0

with final_manifest.open("w") as out:
    print(
        "final_mag_id",
        "source_refined_bin_id",
        "source_refined_bin_fasta",
        "final_mag_fasta",
        "contig_count",
        "total_bp",
        "final_database_mode",
        sep="\\t",
        file=out
    )

    for idx, (bin_id, fasta) in enumerate(rows, start=1):
        fasta = Path(fasta)

        if not fasta.is_absolute():
            fasta = selected_manifest_dir / fasta

        fasta = fasta.resolve()

        if not fasta.exists():
            log(f"WARNING: missing source MAG FASTA: {fasta}")
            missing += 1
            continue

        base = safe_id(bin_id or f"final_mag_{idx:06d}")
        name = f"{base}.fa"
        suffix = 1

        while name in seen_names:
            suffix += 1
            name = f"{base}_{suffix}.fa"

        seen_names.add(name)
        dest = final_dir / name

        copy_as_plain_fasta(fasta, dest)

        contigs, bp = fasta_stats(dest)

        copied += 1
        total_contigs += contigs
        total_bp += bp

        print(
            name.replace(".fa", ""),
            bin_id,
            str(fasta),
            str(published_final_dir / name),
            contigs,
            bp,
            final_database_mode,
            sep="\\t",
            file=out
        )

with final_stats.open("w") as out:
    print(
        "selected_manifest",
        "final_database_mode",
        "manifest_rows",
        "mags_copied",
        "missing_mags",
        "total_contigs",
        "total_bp",
        "final_mag_database_dir",
        sep="\\t",
        file=out
    )

    print(
        str(selected_manifest),
        final_database_mode,
        manifest_rows,
        copied,
        missing,
        total_contigs,
        total_bp,
        str(published_final_dir),
        sep="\\t",
        file=out
    )

log(f"Selected manifest: {selected_manifest}")
log(f"Final database mode: {final_database_mode}")
log(f"Manifest rows: {manifest_rows}")
log(f"MAGs copied: {copied}")
log(f"Missing MAGs: {missing}")
log(f"Total contigs: {total_contigs}")
log(f"Total bp: {total_bp}")
PY

    echo "Final MAG database construction finished: \$(date)" >> "\$LOG_FILE"
    """
}
