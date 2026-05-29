#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Module 3: MAG binning from Module 2 assemblies and Module 1 trimmed reads.
 *
 * Expected inputs:
 *
 *   <working_dir>/module_1_readtrimming/summary/trimmed_manifest.tsv
 *   <working_dir>/module_2_readassembly/summary/assembly_manifest.tsv
 *
 * Supported binners:
 *
 *   --metabat2
 *   --quickbin
 *   --maxbin2
 *
 * Example:
 *
 *   nextflow run module_3_binning.nf \
 *     --working_dir ./output_samwise \
 *     --threads 6 \
 *     --metabat2 \
 *     --quickbin \
 *     --maxbin2
 */


/*
 * Parameters
 */

params.working_dir = null
params.input_assembly_manifest = null
params.input_trimmed_manifest  = null
params.output_dir = null

params.metabat2 = false
params.quickbin = false
params.maxbin2  = false

params.auto_install = true
params.tool_env_dir = null

params.threads = null
params.mapping_threads = 4
params.binning_threads = 4

/*
 * Tool versions.
 */

params.seqkit_version   = "2.8.2"
params.bbmap_version    = "39.81"
params.samtools_version = "1.23.1"
params.metabat2_version = "2.18"
params.maxbin2_version  = "2.2.7"

/*
 * Assembly filtering.
 */
params.min_scaffold_length = 2500

/*
 * BBMap mapping parameters.
 *
 * bbmap_xmx is only used for the direct Java BBMap call. This avoids the
 * bbmap.sh wrapper classpath bug when the env is inside a path with spaces.
 */
params.bbmap_minid      = 0.90
params.bbmap_maxindel   = 10
params.bbmap_ambig      = "random"
params.bbmap_mateqtag   = true
params.bbmap_extra_args = ""
params.bbmap_xmx        = "4g"

/*
 * MetaBAT2 parameters.
 */
params.metabat2_min_contig = 2500
params.metabat2_extra_args = ""

/*
 * MaxBin2 parameters.
 */
params.maxbin2_extra_args = ""

/*
 * QuickBin parameters.
 */
params.quickbin_mincluster = "50k"
params.quickbin_mincontig  = 2500
params.quickbin_minseed    = 2500
params.quickbin_stringency = "normal"
params.quickbin_gzip       = false
params.quickbin_chaff      = false
params.quickbin_clade      = false
params.quickbin_sketch     = false
params.quickbin_server     = false
params.quickbin_xmx        = null
params.quickbin_extra_args = ""
params.quickbin_use_positional_bam = false

/*
 * Publish modes.
 */
params.publish_filtered_assemblies_mode = "symlink"
params.publish_bam_mode  = "symlink"
params.publish_bins_mode = "symlink"

params.results_dir    = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.module1_outdir = "${params.results_dir}/module_1_readtrimming"
params.module2_outdir = "${params.results_dir}/module_2_readassembly"
params.outdir         = "${params.results_dir}/module_3_binning"


def absOrEmpty(value) {
    def s = value == null ? "" : value.toString().trim()

    if( !s || s == "null" || s == "NA" ) {
        return ""
    }

    return java.nio.file.Paths.get(s).toAbsolutePath().normalize().toString()
}


workflow {

    def use_metabat2 = params.metabat2.toString().toBoolean()
    def use_quickbin = params.quickbin.toString().toBoolean()
    def use_maxbin2  = params.maxbin2.toString().toBoolean()

    if( !use_metabat2 && !use_quickbin && !use_maxbin2 ) {
        error """
        No binner selected.

        Please specify at least one of:
          --metabat2
          --quickbin
          --maxbin2

        Example:
          nextflow run module_3_binning.nf --working_dir ./output_samwise --threads 6 --metabat2 --quickbin --maxbin2
        """.stripIndent()
    }

    def assembly_manifest_file = params.input_assembly_manifest ?: "${params.module2_outdir}/summary/assembly_manifest.tsv"
    def trimmed_manifest_file  = params.input_trimmed_manifest  ?: "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    log.info "Module 3 results directory: ${params.results_dir}"
    log.info "Using Module 2 assembly manifest: ${assembly_manifest_file}"
    log.info "Using Module 1 trimmed manifest: ${trimmed_manifest_file}"
    log.info "Writing Module 3 outputs to: ${params.outdir}"
    log.info "Module 3 conda environments directory: ${params.tool_env_dir ?: "${params.outdir}/conda_envs"}"
    log.info "Minimum scaffold length for binning: ${params.min_scaffold_length}"
    log.info "Binners selected: MetaBAT2=${use_metabat2}, QuickBin=${use_quickbin}, MaxBin2=${use_maxbin2}"

    def assembly_manifest_ch = channel.fromPath(
        assembly_manifest_file,
        type: 'file',
        checkIfExists: true
    )

    def trimmed_manifest_ch = channel.fromPath(
        trimmed_manifest_file,
        type: 'file',
        checkIfExists: true
    )

    /*
     * Module 2 assembly manifest columns:
     *
     * sample_id
     * safe_sample_id
     * assembly_sample_id
     * assembler
     * assembly_mode
     * rarefaction_label
     * assembly_strategy
     * renamed_fasta
     */
    def assemblies_ch = assembly_manifest_ch
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample_id.toString(),
                row.safe_sample_id.toString(),
                row.assembly_sample_id.toString(),
                row.assembler.toString(),
                row.assembly_mode.toString(),
                row.rarefaction_label.toString(),
                row.assembly_strategy.toString(),
                file(absOrEmpty(row.renamed_fasta))
            )
        }

    /*
     * Module 1 trimmed manifest columns:
     *
     * sample_id
     * safe_sample_id
     * layout
     * read1
     * read2
     * interleaved
     * merged
     * fastp_html
     * fastp_json
     */
    def trimmed_reads_ch = trimmed_manifest_ch
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sample_id.toString(),
                row.layout.toString(),
                absOrEmpty(row.read1),
                absOrEmpty(row.read2),
                absOrEmpty(row.interleaved)
            )
        }

    /*
     * Join assemblies to their original sample reads by sample_id.
     */
    def binning_jobs_ch = assemblies_ch
        .join(trimmed_reads_ch)
        .map {
            sample_id,
            safe_sample_id,
            assembly_sample_id,
            assembly_assembler,
            assembly_mode,
            rarefaction_label,
            assembly_strategy,
            renamed_fasta,
            layout,
            read1,
            read2,
            interleaved ->

            def rare_part = rarefaction_label ? "_${rarefaction_label}" : ""
            def binning_id = "${safe_sample_id}${rare_part}_${assembly_assembler}_${assembly_mode}"
                .replaceAll('[^A-Za-z0-9._-]+', '_')

            if( layout != 'paired' && layout != 'interleaved' ) {
                error "Unsupported layout in trimmed manifest for sample '${sample_id}': ${layout}"
            }

            tuple(
                sample_id,
                safe_sample_id,
                assembly_sample_id,
                assembly_assembler,
                assembly_mode,
                rarefaction_label,
                assembly_strategy,
                binning_id,
                renamed_fasta,
                layout,
                read1,
                read2,
                interleaved
            )
        }

    SETUP_MODULE3_TOOLS()

    FILTER_ASSEMBLY_BY_LENGTH(
        binning_jobs_ch.combine(SETUP_MODULE3_TOOLS.out.status)
    )

    WRITE_FILTERED_ASSEMBLY_STATS_SUMMARY(
        FILTER_ASSEMBLY_BY_LENGTH.out.stats_file.collect()
    )

    MAP_READS_BBMAP(
        FILTER_ASSEMBLY_BY_LENGTH.out.filtered_jobs.combine(SETUP_MODULE3_TOOLS.out.status)
    )

    def manifest_records_ch = channel.empty()
    def stats_files_ch      = channel.empty()

    if( use_metabat2 || use_maxbin2 ) {
        GENERATE_COVERAGE_FILES(
            MAP_READS_BBMAP.out.mapped_bam.combine(SETUP_MODULE3_TOOLS.out.status)
        )
    }

    if( use_metabat2 ) {
        RUN_METABAT2(
            GENERATE_COVERAGE_FILES.out.coverage_files.combine(SETUP_MODULE3_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(RUN_METABAT2.out.manifest_record)
        stats_files_ch      = stats_files_ch.mix(RUN_METABAT2.out.stats_file)
    }

    if( use_maxbin2 ) {
        RUN_MAXBIN2(
            GENERATE_COVERAGE_FILES.out.coverage_files.combine(SETUP_MODULE3_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(RUN_MAXBIN2.out.manifest_record)
        stats_files_ch      = stats_files_ch.mix(RUN_MAXBIN2.out.stats_file)
    }

    if( use_quickbin ) {
        RUN_QUICKBIN(
            MAP_READS_BBMAP.out.mapped_bam.combine(SETUP_MODULE3_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(RUN_QUICKBIN.out.manifest_record)
        stats_files_ch      = stats_files_ch.mix(RUN_QUICKBIN.out.stats_file)
    }

    WRITE_BINNING_SUMMARIES(
        manifest_records_ch.collect(),
        stats_files_ch.collect()
    )
}

process SETUP_MODULE3_TOOLS {
    tag "setup_binning_tools"

    publishDir "${params.outdir}/setup",
        mode: 'copy',
        pattern: "module3_tools_status.env"

    output:
    path "module3_tools_status.env", emit: status

    script:
    def base_env_dir = params.tool_env_dir ?: "${params.outdir}/conda_envs"
    def mapping_env_dir  = "${base_env_dir}/module3_mapping_tools"
    def metabat2_env_dir = "${base_env_dir}/module3_metabat2_tools"
    def maxbin2_env_dir  = "${base_env_dir}/module3_maxbin2_tools"

    def want_metabat2 = params.metabat2.toString().toBoolean()
    def want_quickbin = params.quickbin.toString().toBoolean()
    def want_maxbin2  = params.maxbin2.toString().toBoolean()

    def need_jgi_depth = want_metabat2 || want_maxbin2

    """
    set -euo pipefail

    STATUS_FILE="module3_tools_status.env"

    MAPPING_ENV="${mapping_env_dir}"
    METABAT2_ENV="${metabat2_env_dir}"
    MAXBIN2_ENV="${maxbin2_env_dir}"

    WANT_METABAT2="${want_metabat2}"
    WANT_QUICKBIN="${want_quickbin}"
    WANT_MAXBIN2="${want_maxbin2}"
    NEED_JGI_DEPTH="${need_jgi_depth}"

    # Use bash arrays so package names are passed to mamba/conda reliably.
    # IMPORTANT FIX:
    # openjdk=17.* is explicitly included in MAPPING_PACKAGES.
    MAPPING_PACKAGES=(
        "python=3.*"
        "openjdk=17.*"
        "seqkit=${params.seqkit_version}"
        "bbmap=${params.bbmap_version}"
        "samtools>=${params.samtools_version},<2.0a0"
    )

    METABAT2_PACKAGES=(
        "python=3.*"
        "metabat2=${params.metabat2_version}"
    )

    MAXBIN2_PACKAGES=(
        "python=3.*"
        "perl"
        "maxbin2=${params.maxbin2_version}"
    )

    echo "Module 3 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested MetaBAT2: \$WANT_METABAT2" >> "\$STATUS_FILE"
    echo "Requested QuickBin: \$WANT_QUICKBIN" >> "\$STATUS_FILE"
    echo "Requested MaxBin2: \$WANT_MAXBIN2" >> "\$STATUS_FILE"
    echo "JGI depth utility needed: \$NEED_JGI_DEPTH" >> "\$STATUS_FILE"

    echo "Requested SeqKit version: ${params.seqkit_version}" >> "\$STATUS_FILE"
    echo "Requested BBMap/BBTools version: ${params.bbmap_version}" >> "\$STATUS_FILE"
    echo "Requested Samtools lower bound: >=${params.samtools_version},<2.0a0" >> "\$STATUS_FILE"
    echo "Requested OpenJDK version: 17.*" >> "\$STATUS_FILE"
    echo "Requested MetaBAT2 version: ${params.metabat2_version}" >> "\$STATUS_FILE"
    echo "Requested MaxBin2 version: ${params.maxbin2_version}" >> "\$STATUS_FILE"

    echo "MAPPING_ENV=\$MAPPING_ENV" >> "\$STATUS_FILE"
    echo "QUICKBIN_ENV=\$MAPPING_ENV" >> "\$STATUS_FILE"
    echo "METABAT2_ENV=\$METABAT2_ENV" >> "\$STATUS_FILE"
    echo "JGI_ENV=\$METABAT2_ENV" >> "\$STATUS_FILE"
    echo "MAXBIN2_ENV=\$MAXBIN2_ENV" >> "\$STATUS_FILE"

    echo "Mapping package request:" >> "\$STATUS_FILE"
    printf '  %s\\n' "\${MAPPING_PACKAGES[@]}" >> "\$STATUS_FILE"

    echo "MetaBAT2 package request:" >> "\$STATUS_FILE"
    printf '  %s\\n' "\${METABAT2_PACKAGES[@]}" >> "\$STATUS_FILE"

    echo "MaxBin2 package request:" >> "\$STATUS_FILE"
    printf '  %s\\n' "\${MAXBIN2_PACKAGES[@]}" >> "\$STATUS_FILE"

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

    INSTALLER="\$(find_installer)"

    if [[ -z "\$INSTALLER" ]]; then
        echo "ERROR: Neither mamba nor conda found in PATH." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "Using installer: \$INSTALLER" >> "\$STATUS_FILE"

    exists_tool() {
        local prefix="\$1"
        local exe="\$2"

        if [[ "\$prefix" == "SYSTEM" ]]; then
            command -v "\$exe" >/dev/null 2>&1
        else
            [[ -x "\$prefix/bin/\$exe" ]]
        fi
    }

    run_tool() {
        local prefix="\$1"
        local exe="\$2"
        shift 2

        if [[ "\$prefix" == "SYSTEM" ]]; then
            "\$exe" "\$@" >> "\$STATUS_FILE" 2>&1
        else
            "\$prefix/bin/\$exe" "\$@" >> "\$STATUS_FILE" 2>&1
        fi
    }

    check_mapping_env() {
    local prefix="\$1"
    local ok="true"

    echo "Checking mapping environment: \$prefix" >> "\$STATUS_FILE"

    if exists_tool "\$prefix" seqkit; then
        run_tool "\$prefix" seqkit version || ok="false"
    else
        echo "Missing seqkit" >> "\$STATUS_FILE"
        ok="false"
    fi

    if exists_tool "\$prefix" bbmap.sh; then
        run_tool "\$prefix" bbmap.sh -h >/dev/null 2>&1 || true
    else
        echo "Missing bbmap.sh" >> "\$STATUS_FILE"
        ok="false"
    fi

    if exists_tool "\$prefix" quickbin.sh; then
        run_tool "\$prefix" quickbin.sh -h >/dev/null 2>&1 || true
    elif [[ "\$WANT_QUICKBIN" == "true" ]]; then
        echo "Missing quickbin.sh" >> "\$STATUS_FILE"
        ok="false"
    fi

    if exists_tool "\$prefix" samtools; then
        run_tool "\$prefix" samtools --version || ok="false"
    else
        echo "Missing samtools" >> "\$STATUS_FILE"
        ok="false"
    fi

    # Java handling for conda-forge openjdk on macOS.
    #
    # On your system, openjdk installed java here:
    #
    #   \$prefix/lib/jvm/bin/java
    #
    # but the workflow expects:
    #
    #   \$prefix/bin/java
    #
    # So if the lib/jvm copy exists, create a compatibility symlink.
    if [[ -x "\$prefix/bin/java" ]]; then
        echo "Found java at: \$prefix/bin/java" >> "\$STATUS_FILE"
        run_tool "\$prefix" java -version || true

    elif [[ -x "\$prefix/lib/jvm/bin/java" ]]; then
        echo "Found java at: \$prefix/lib/jvm/bin/java" >> "\$STATUS_FILE"
        echo "Creating compatibility symlink: \$prefix/bin/java -> \$prefix/lib/jvm/bin/java" >> "\$STATUS_FILE"

        mkdir -p "\$prefix/bin"
        ln -sfn "\$prefix/lib/jvm/bin/java" "\$prefix/bin/java"

        if [[ -x "\$prefix/bin/java" ]]; then
            echo "Java compatibility symlink created successfully." >> "\$STATUS_FILE"
            run_tool "\$prefix" java -version || true
        else
            echo "ERROR: Java symlink was created but is not executable." >> "\$STATUS_FILE"
            ok="false"
        fi

    else
        echo "Missing java" >> "\$STATUS_FILE"
        echo "Checked these locations:" >> "\$STATUS_FILE"
        echo "  \$prefix/bin/java" >> "\$STATUS_FILE"
        echo "  \$prefix/lib/jvm/bin/java" >> "\$STATUS_FILE"

        echo "Searching for java under prefix:" >> "\$STATUS_FILE"
        find "\$prefix" -path '*/bin/java' -print >> "\$STATUS_FILE" 2>&1 || true

        ok="false"
    fi

    local jar_count=0
    if [[ "\$prefix" != "SYSTEM" ]]; then
        jar_count=\$(find "\$prefix" -type f -name 'bbtools.jar' 2>/dev/null | wc -l | tr -d ' ')
        echo "bbtools.jar files found under mapping env: \$jar_count" >> "\$STATUS_FILE"

        if [[ "\$jar_count" -lt 1 ]]; then
            echo "Missing bbtools.jar under mapping env." >> "\$STATUS_FILE"
            ok="false"
        fi
    fi

    [[ "\$ok" == "true" ]]
}

    check_metabat2_env() {
        local prefix="\$1"
        local ok="true"

        echo "Checking MetaBAT2/JGI environment: \$prefix" >> "\$STATUS_FILE"

        if [[ "\$WANT_METABAT2" == "true" ]]; then
            if exists_tool "\$prefix" metabat2; then
                run_tool "\$prefix" metabat2 --help >/dev/null 2>&1 || true
            else
                echo "Missing metabat2" >> "\$STATUS_FILE"
                ok="false"
            fi
        fi

        if [[ "\$NEED_JGI_DEPTH" == "true" ]]; then
            if exists_tool "\$prefix" jgi_summarize_bam_contig_depths; then
                echo "Found jgi_summarize_bam_contig_depths" >> "\$STATUS_FILE"
            else
                echo "Missing jgi_summarize_bam_contig_depths" >> "\$STATUS_FILE"
                ok="false"
            fi
        fi

        [[ "\$ok" == "true" ]]
    }

    check_maxbin2_env() {
        local prefix="\$1"
        local ok="true"

        echo "Checking MaxBin2 environment: \$prefix" >> "\$STATUS_FILE"

        if [[ "\$WANT_MAXBIN2" == "true" ]]; then
            if exists_tool "\$prefix" run_MaxBin.pl; then
                run_tool "\$prefix" run_MaxBin.pl -h >/dev/null 2>&1 || true
            else
                echo "Missing run_MaxBin.pl" >> "\$STATUS_FILE"
                ok="false"
            fi
        fi

        [[ "\$ok" == "true" ]]
    }

    create_env() {
        local env_path="\$1"
        shift

        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: auto_install is false and environment is missing or broken: \$env_path" >> "\$STATUS_FILE"
            exit 1
        fi

        mkdir -p "\$(dirname "\$env_path")"

        echo "Creating environment: \$env_path" >> "\$STATUS_FILE"
        echo "Packages:" >> "\$STATUS_FILE"
        printf '  %s\\n' "\$@" >> "\$STATUS_FILE"

        "\$INSTALLER" create -y \\
            -p "\$env_path" \\
            -c conda-forge \\
            -c bioconda \\
            "\$@" \\
            >> "\$STATUS_FILE" 2>&1
    }

    echo "" >> "\$STATUS_FILE"
    echo "=== Mapping / QuickBin environment ===" >> "\$STATUS_FILE"

    if [[ -d "\$MAPPING_ENV" ]]; then
        if check_mapping_env "\$MAPPING_ENV"; then
            echo "Existing mapping environment passed checks." >> "\$STATUS_FILE"
        else
            echo "Existing mapping environment failed checks. Removing." >> "\$STATUS_FILE"
            rm -rf "\$MAPPING_ENV"
        fi
    fi

    if [[ ! -d "\$MAPPING_ENV" ]]; then
        create_env "\$MAPPING_ENV" "\${MAPPING_PACKAGES[@]}"

        if ! check_mapping_env "\$MAPPING_ENV"; then
            echo "ERROR: Newly created mapping environment failed checks." >> "\$STATUS_FILE"
            echo "Contents of \$MAPPING_ENV/bin, if present:" >> "\$STATUS_FILE"
            ls -lah "\$MAPPING_ENV/bin" >> "\$STATUS_FILE" 2>&1 || true
            exit 1
        fi
    fi

    echo "" >> "\$STATUS_FILE"
    echo "=== MetaBAT2 / JGI depth environment ===" >> "\$STATUS_FILE"

    if [[ "\$NEED_JGI_DEPTH" == "true" ]]; then
        if [[ -d "\$METABAT2_ENV" ]]; then
            if check_metabat2_env "\$METABAT2_ENV"; then
                echo "Existing MetaBAT2 environment passed checks." >> "\$STATUS_FILE"
            else
                echo "Existing MetaBAT2 environment failed checks. Removing." >> "\$STATUS_FILE"
                rm -rf "\$METABAT2_ENV"
            fi
        fi

        if [[ ! -d "\$METABAT2_ENV" ]]; then
            create_env "\$METABAT2_ENV" "\${METABAT2_PACKAGES[@]}"

            if ! check_metabat2_env "\$METABAT2_ENV"; then
                echo "ERROR: Newly created MetaBAT2 environment failed checks." >> "\$STATUS_FILE"
                exit 1
            fi
        fi
    else
        echo "MetaBAT2/JGI environment not needed for selected binners." >> "\$STATUS_FILE"
        METABAT2_ENV="NOT_USED"
    fi

    echo "" >> "\$STATUS_FILE"
    echo "=== MaxBin2 environment ===" >> "\$STATUS_FILE"

    if [[ "\$WANT_MAXBIN2" == "true" ]]; then
        if [[ -d "\$MAXBIN2_ENV" ]]; then
            if check_maxbin2_env "\$MAXBIN2_ENV"; then
                echo "Existing MaxBin2 environment passed checks." >> "\$STATUS_FILE"
            else
                echo "Existing MaxBin2 environment failed checks. Removing." >> "\$STATUS_FILE"
                rm -rf "\$MAXBIN2_ENV"
            fi
        fi

        if [[ ! -d "\$MAXBIN2_ENV" ]]; then
            create_env "\$MAXBIN2_ENV" "\${MAXBIN2_PACKAGES[@]}"

            if ! check_maxbin2_env "\$MAXBIN2_ENV"; then
                echo "ERROR: Newly created MaxBin2 environment failed checks." >> "\$STATUS_FILE"
                exit 1
            fi
        fi
    else
        echo "MaxBin2 environment not needed for selected binners." >> "\$STATUS_FILE"
        MAXBIN2_ENV="NOT_USED"
    fi

    echo "" >> "\$STATUS_FILE"
    echo "Final environment assignments:" >> "\$STATUS_FILE"
    echo "MAPPING_ENV=\$MAPPING_ENV" >> "\$STATUS_FILE"
    echo "QUICKBIN_ENV=\$MAPPING_ENV" >> "\$STATUS_FILE"
    echo "METABAT2_ENV=\$METABAT2_ENV" >> "\$STATUS_FILE"
    echo "JGI_ENV=\$METABAT2_ENV" >> "\$STATUS_FILE"
    echo "MAXBIN2_ENV=\$MAXBIN2_ENV" >> "\$STATUS_FILE"
    echo "Module 3 tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}


process FILTER_ASSEMBLY_BY_LENGTH {

    tag { binning_id }

    stageInMode 'copy'

    publishDir "${params.outdir}/filtered_assemblies",
        mode: params.publish_filtered_assemblies_mode,
        pattern: "*.min${params.min_scaffold_length}.fa"

    publishDir "${params.outdir}/summary/filtered_assembly_stats",
        mode: 'copy',
        pattern: "*.filtered_assembly_stats.tsv"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.filter_assembly.log"

    input:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(original_assembly_fasta),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          path(tools_status)

    output:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path("${binning_id}.min${params.min_scaffold_length}.fa"),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          emit: filtered_jobs

    path "${binning_id}.filtered_assembly_stats.tsv", emit: stats_file
    path "${binning_id}.filter_assembly.log", emit: log_file

    script:
    """
    set -euo pipefail

    MAPPING_ENV="\$(grep '^MAPPING_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$MAPPING_ENV" && "\$MAPPING_ENV" != "SYSTEM" && "\$MAPPING_ENV" != "NOT_USED" ]]; then
        export PATH="\$MAPPING_ENV/bin:\$PATH"
    fi

    if ! command -v seqkit >/dev/null 2>&1; then
        echo "ERROR: seqkit is not available after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    LOG_FILE="${binning_id}.filter_assembly.log"
    OUT_FASTA="${binning_id}.min${params.min_scaffold_length}.fa"
    STATS_FILE="${binning_id}.filtered_assembly_stats.tsv"

    echo "Assembly length filtering started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: ${sample_id}" >> "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Original assembly: ${original_assembly_fasta}" >> "\$LOG_FILE"
    echo "Minimum scaffold length: ${params.min_scaffold_length}" >> "\$LOG_FILE"
    echo "Filtered assembly: \$OUT_FASTA" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    seqkit seq \\
        -m ${params.min_scaffold_length} \\
        "${original_assembly_fasta}" \\
        > "\$OUT_FASTA"

    python3 - \\
        "${original_assembly_fasta}" \\
        "\$OUT_FASTA" \\
        "\$STATS_FILE" \\
        "${sample_id}" \\
        "${safe_sample_id}" \\
        "${assembly_sample_id}" \\
        "${assembly_assembler}" \\
        "${assembly_mode}" \\
        "${rarefaction_label}" \\
        "${assembly_strategy}" \\
        "${binning_id}" \\
        "${params.min_scaffold_length}" \\
        "${params.outdir}/filtered_assemblies/\$OUT_FASTA" <<'PY'
import gzip
import sys
from pathlib import Path

(
    original_fasta,
    filtered_fasta,
    stats_file,
    sample_id,
    safe_sample_id,
    assembly_sample_id,
    assembly_assembler,
    assembly_mode,
    rarefaction_label,
    assembly_strategy,
    binning_id,
    min_len,
    published_filtered_fasta
) = sys.argv[1:]

def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def fasta_lengths(path):
    lengths = []
    current = 0
    seen = False

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\\n")
            if line.startswith(">"):
                if seen:
                    lengths.append(current)
                current = 0
                seen = True
            else:
                current += len(line.strip())

    if seen:
        lengths.append(current)

    return lengths

orig_lengths = fasta_lengths(original_fasta)
filt_lengths = fasta_lengths(filtered_fasta)

orig_count = len(orig_lengths)
filt_count = len(filt_lengths)

orig_bp = sum(orig_lengths)
filt_bp = sum(filt_lengths)

removed_count = orig_count - filt_count
removed_bp = orig_bp - filt_bp

with open(stats_file, "w") as out:
    print(
        "sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembly_assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "binning_id",
        "min_scaffold_length",
        "original_contigs",
        "original_bp",
        "filtered_contigs",
        "filtered_bp",
        "removed_contigs",
        "removed_bp",
        "filtered_fasta",
        sep="\\t",
        file=out
    )
    print(
        sample_id,
        safe_sample_id,
        assembly_sample_id,
        assembly_assembler,
        assembly_mode,
        rarefaction_label,
        assembly_strategy,
        binning_id,
        min_len,
        orig_count,
        orig_bp,
        filt_count,
        filt_bp,
        removed_count,
        removed_bp,
        published_filtered_fasta,
        sep="\\t",
        file=out
    )

if filt_count == 0:
    raise SystemExit(
        f"ERROR: No contigs >= {min_len} bp remained after filtering for {binning_id}. "
        f"Original contigs: {orig_count}, original bp: {orig_bp}."
    )
PY

    echo "Assembly length filtering finished: \$(date)" >> "\$LOG_FILE"
    """
}


process MAP_READS_BBMAP {

    tag { "${sample_id}:${binning_id}" }

    stageInMode 'symlink'

    publishDir "${params.outdir}/mapping",
        mode: params.publish_bam_mode,
        pattern: "*.sorted.bam*"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.bbmap.log"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.mapping_threads as int
    }

    input:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(filtered_assembly_fasta),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          path(tools_status)

    output:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(filtered_assembly_fasta),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          path("${binning_id}.sorted.bam"),
          path("${binning_id}.sorted.bam.bai"),
          emit: mapped_bam

    path "${binning_id}.bbmap.log", emit: log_file

    script:
    """
    set -euo pipefail

    MAPPING_ENV="\$(grep '^MAPPING_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$MAPPING_ENV" && "\$MAPPING_ENV" != "SYSTEM" && "\$MAPPING_ENV" != "NOT_USED" ]]; then
        export PATH="\$MAPPING_ENV/bin:\$PATH"
    fi

    if ! command -v samtools >/dev/null 2>&1; then
        echo "ERROR: samtools is not available after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    if ! command -v java >/dev/null 2>&1; then
        echo "ERROR: java is not available after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    LOG_FILE="${binning_id}.bbmap.log"

    echo "BBMap mapping started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: ${sample_id}" >> "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "Layout: ${layout}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "MAPPING_ENV: \$MAPPING_ENV" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    cp -L "${filtered_assembly_fasta}" assembly.fa

    if [[ ! -s assembly.fa ]]; then
        echo "ERROR: Filtered assembly could not be copied or is empty." >> "\$LOG_FILE"
        echo "filtered_assembly_fasta=${filtered_assembly_fasta}" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    BBTOOLS_JAR=""

    if [[ -n "\$MAPPING_ENV" && "\$MAPPING_ENV" != "SYSTEM" && "\$MAPPING_ENV" != "NOT_USED" ]]; then
        BBTOOLS_JAR="\$(find "\$MAPPING_ENV" -type f -name 'bbtools.jar' 2>/dev/null | head -n 1 || true)"
    fi

    if [[ -z "\$BBTOOLS_JAR" ]]; then
        echo "ERROR: Could not find bbtools.jar under MAPPING_ENV." >> "\$LOG_FILE"
        echo "MAPPING_ENV=\$MAPPING_ENV" >> "\$LOG_FILE"
        echo "This module uses direct Java calls to avoid bbmap.sh classpath issues with spaces." >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    echo "Using direct Java BBMap call." >> "\$LOG_FILE"
    echo "BBTOOLS_JAR: \$BBTOOLS_JAR" >> "\$LOG_FILE"

    BBMAP_COMMON_ARGS="ref=assembly.fa ambig=${params.bbmap_ambig} minid=${params.bbmap_minid} maxindel=${params.bbmap_maxindel} threads=${task.cpus} overwrite=t nodisk=t"

    if [[ "${params.bbmap_mateqtag.toString().toLowerCase()}" == "true" ]]; then
        BBMAP_COMMON_ARGS="\${BBMAP_COMMON_ARGS} mateqtag=t"
    fi

    if [[ -n "${params.bbmap_extra_args}" ]]; then
        BBMAP_COMMON_ARGS="\${BBMAP_COMMON_ARGS} ${params.bbmap_extra_args}"
    fi

    echo "BBMap common args: \$BBMAP_COMMON_ARGS" >> "\$LOG_FILE"

    run_bbmap_paired() {
        java \\
            -ea \\
            -Xmx${params.bbmap_xmx} \\
            -cp "\$BBTOOLS_JAR" \\
            align2.BBMap \\
            build=1 \\
            overwrite=true \\
            fastareadlen=500 \\
            \$BBMAP_COMMON_ARGS \\
            in1=input_R1.fastq.gz \\
            in2=input_R2.fastq.gz \\
            out=mapped.sam
    }

    run_bbmap_interleaved() {
        java \\
            -ea \\
            -Xmx${params.bbmap_xmx} \\
            -cp "\$BBTOOLS_JAR" \\
            align2.BBMap \\
            build=1 \\
            overwrite=true \\
            fastareadlen=500 \\
            \$BBMAP_COMMON_ARGS \\
            in=input_interleaved.fastq.gz \\
            interleaved=t \\
            out=mapped.sam
    }

    if [[ "${layout}" == "paired" ]]; then

        if [[ ! -s "${read1}" || ! -s "${read2}" ]]; then
            echo "ERROR: Missing paired read files for sample ${sample_id}" >> "\$LOG_FILE"
            echo "read1=${read1}" >> "\$LOG_FILE"
            echo "read2=${read2}" >> "\$LOG_FILE"
            cat "\$LOG_FILE" >&2
            exit 1
        fi

        ln -sfn "${read1}" input_R1.fastq.gz
        ln -sfn "${read2}" input_R2.fastq.gz

        if [[ ! -s input_R1.fastq.gz || ! -s input_R2.fastq.gz ]]; then
            echo "ERROR: Local paired-read symlinks are missing or empty." >> "\$LOG_FILE"
            echo "read1=${read1}" >> "\$LOG_FILE"
            echo "read2=${read2}" >> "\$LOG_FILE"
            ls -lah >> "\$LOG_FILE" 2>&1 || true
            cat "\$LOG_FILE" >&2
            exit 1
        fi

        echo "Running BBMap paired-end mapping..." >> "\$LOG_FILE"

        if ! run_bbmap_paired >> "\$LOG_FILE" 2>&1; then
            echo "ERROR: BBMap paired-end mapping failed." >&2
            echo "Last 100 lines of \$LOG_FILE:" >&2
            tail -n 100 "\$LOG_FILE" >&2 || true
            exit 1
        fi

    elif [[ "${layout}" == "interleaved" ]]; then

        if [[ ! -s "${interleaved}" ]]; then
            echo "ERROR: Missing interleaved read file for sample ${sample_id}" >> "\$LOG_FILE"
            echo "interleaved=${interleaved}" >> "\$LOG_FILE"
            cat "\$LOG_FILE" >&2
            exit 1
        fi

        ln -sfn "${interleaved}" input_interleaved.fastq.gz

        if [[ ! -s input_interleaved.fastq.gz ]]; then
            echo "ERROR: Local interleaved-read symlink is missing or empty." >> "\$LOG_FILE"
            echo "interleaved=${interleaved}" >> "\$LOG_FILE"
            ls -lah >> "\$LOG_FILE" 2>&1 || true
            cat "\$LOG_FILE" >&2
            exit 1
        fi

        echo "Running BBMap interleaved mapping..." >> "\$LOG_FILE"

        if ! run_bbmap_interleaved >> "\$LOG_FILE" 2>&1; then
            echo "ERROR: BBMap interleaved mapping failed." >&2
            echo "Last 100 lines of \$LOG_FILE:" >&2
            tail -n 100 "\$LOG_FILE" >&2 || true
            exit 1
        fi

    else
        echo "ERROR: Unsupported layout: ${layout}" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    if [[ ! -s mapped.sam ]]; then
        echo "ERROR: BBMap did not produce a non-empty SAM file." >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    echo "Converting SAM to sorted BAM..." >> "\$LOG_FILE"

    if ! samtools view -@ ${task.cpus} -bS mapped.sam \\
        | samtools sort -@ ${task.cpus} -o "${binning_id}.sorted.bam" -
    then
        echo "ERROR: samtools view/sort failed." >&2
        echo "Last 100 lines of \$LOG_FILE:" >&2
        tail -n 100 "\$LOG_FILE" >&2 || true
        exit 1
    fi

    if [[ ! -s "${binning_id}.sorted.bam" ]]; then
        echo "ERROR: sorted BAM was not created or is empty." >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    samtools index "${binning_id}.sorted.bam"

    rm -f mapped.sam assembly.fa input_R1.fastq.gz input_R2.fastq.gz input_interleaved.fastq.gz

    echo "BBMap mapping finished: \$(date)" >> "\$LOG_FILE"
    """
}

process GENERATE_COVERAGE_FILES {

    tag { binning_id }

    stageInMode 'symlink'

    publishDir "${params.outdir}/coverage/metabat2_depth",
        mode: 'copy',
        pattern: "*.depth.txt"

    publishDir "${params.outdir}/coverage/maxbin2_abundance",
        mode: 'copy',
        pattern: "*.maxbin2_abundance.tsv"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.coverage.log"

    input:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(filtered_assembly_fasta),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          path(sorted_bam),
          path(sorted_bam_bai),
          path(tools_status)

    output:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(filtered_assembly_fasta),
          path(sorted_bam),
          path(sorted_bam_bai),
          path("${binning_id}.depth.txt"),
          path("${binning_id}.maxbin2_abundance.tsv"),
          emit: coverage_files

    path "${binning_id}.coverage.log", emit: log_file

    script:
    """
    set -euo pipefail

    JGI_ENV="\$(grep '^JGI_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$JGI_ENV" && "\$JGI_ENV" != "SYSTEM" && "\$JGI_ENV" != "NOT_USED" ]]; then
        export PATH="\$JGI_ENV/bin:\$PATH"
    fi

    if ! command -v jgi_summarize_bam_contig_depths >/dev/null 2>&1; then
        echo "ERROR: jgi_summarize_bam_contig_depths is not available after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    LOG_FILE="${binning_id}.coverage.log"

    echo "Coverage generation started: \$(date)" > "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "BAM: ${sorted_bam}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    jgi_summarize_bam_contig_depths \\
        --outputDepth "${binning_id}.depth.txt" \\
        "${sorted_bam}" \\
        >> "\$LOG_FILE" 2>&1

    if [[ ! -s "${binning_id}.depth.txt" ]]; then
        echo "ERROR: jgi_summarize_bam_contig_depths did not produce a non-empty depth file." >> "\$LOG_FILE"
        exit 1
    fi

    python3 - "${binning_id}.depth.txt" "${binning_id}.maxbin2_abundance.tsv" <<'PY'
import sys
from pathlib import Path

depth_path = Path(sys.argv[1])
abund_path = Path(sys.argv[2])

with depth_path.open() as inp, abund_path.open("w") as out:
    for idx, line in enumerate(inp):
        line = line.rstrip("\\n")
        if not line:
            continue

        fields = line.split("\\t")

        if idx == 0 and fields[0] == "contigName":
            continue

        if len(fields) < 3:
            continue

        contig = fields[0]
        total_avg_depth = fields[2]

        print(contig, total_avg_depth, sep="\\t", file=out)
PY

    if [[ ! -s "${binning_id}.maxbin2_abundance.tsv" ]]; then
        echo "ERROR: MaxBin2 abundance file is empty." >> "\$LOG_FILE"
        exit 1
    fi

    echo "Coverage generation finished: \$(date)" >> "\$LOG_FILE"
    """
}

process RUN_METABAT2 {

    tag { binning_id }

    stageInMode 'symlink'

    publishDir { "${params.outdir}/bins/metabat2/${binning_id}" },
        mode: params.publish_bins_mode,
        pattern: "bins/**",
        saveAs: { filename -> filename.replaceFirst(/^bins\//, '') }

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.metabat2.log"

    publishDir "${params.outdir}/summary/per_binner_stats",
        mode: 'copy',
        pattern: "*.metabat2.binning_stats.tsv"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.binning_threads as int
    }

    input:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(filtered_assembly_fasta),
          path(sorted_bam),
          path(sorted_bam_bai),
          path(depth_file),
          path(maxbin_abundance),
          path(tools_status)

    output:
    path "bins", emit: bin_dir
    path "${binning_id}.metabat2.binning_manifest_record.tsv", emit: manifest_record
    path "${binning_id}.metabat2.binning_stats.tsv", emit: stats_file
    path "${binning_id}.metabat2.log", emit: log_file

    script:
    """
    set -euo pipefail

    METABAT2_ENV="\$(grep '^METABAT2_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$METABAT2_ENV" && "\$METABAT2_ENV" != "SYSTEM" && "\$METABAT2_ENV" != "NOT_USED" ]]; then
        export PATH="\$METABAT2_ENV/bin:\$PATH"
    fi

    if ! command -v metabat2 >/dev/null 2>&1; then
        echo "ERROR: metabat2 is not available after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    BINNER="metabat2"
    BIN_DIR="bins"
    LOG_FILE="${binning_id}.metabat2.log"

    mkdir -p "\$BIN_DIR"

    echo "MetaBAT2 binning started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: ${sample_id}" >> "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "Depth file: ${depth_file}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    metabat2 \\
        -i "${filtered_assembly_fasta}" \\
        -a "${depth_file}" \\
        -o "\$BIN_DIR/${binning_id}.metabat2_bin" \\
        -t ${task.cpus} \\
        -m ${params.metabat2_min_contig} \\
        ${params.metabat2_extra_args} \\
        >> "\$LOG_FILE" 2>&1

    python3 - \\
        "\$BIN_DIR" \\
        "${binning_id}.metabat2.binning_manifest_record.tsv" \\
        "${binning_id}.metabat2.binning_stats.tsv" \\
        "${sample_id}" \\
        "${safe_sample_id}" \\
        "${assembly_sample_id}" \\
        "${assembly_assembler}" \\
        "${assembly_mode}" \\
        "${rarefaction_label}" \\
        "${assembly_strategy}" \\
        "${binning_id}" \\
        "\$BINNER" \\
        "${params.outdir}/bins/metabat2/${binning_id}" \\
        "${params.outdir}/logs/${binning_id}.metabat2.log" \\
        "" <<'PY'
import gzip
import sys
from pathlib import Path

(
    bin_dir,
    manifest_path,
    stats_path,
    sample_id,
    safe_sample_id,
    assembly_sample_id,
    assembly_assembler,
    assembly_mode,
    rarefaction_label,
    assembly_strategy,
    binning_id,
    binner,
    published_bin_dir,
    published_log,
    published_report
) = sys.argv[1:]

bin_dir = Path(bin_dir)

def fasta_len(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total

fasta_exts = {".fa", ".fasta", ".fna", ".fas"}
bin_files = []

for p in sorted(bin_dir.iterdir()):
    if p.is_file():
        suffixes = p.suffixes
        if p.suffix in fasta_exts or (len(suffixes) >= 2 and suffixes[-1] == ".gz" and suffixes[-2] in fasta_exts):
            bin_files.append(p)

total_bp = 0
largest = 0

with open(manifest_path, "w") as manifest:
    for idx, p in enumerate(bin_files, start=1):
        bp = fasta_len(p)
        total_bp += bp
        largest = max(largest, bp)

        bin_id = f"{binning_id}_{binner}_bin_{idx:03d}"
        published_path = str(Path(published_bin_dir) / p.name)

        print(
            sample_id,
            safe_sample_id,
            assembly_sample_id,
            assembly_assembler,
            assembly_mode,
            rarefaction_label,
            assembly_strategy,
            binning_id,
            binner,
            bin_id,
            bp,
            published_path,
            sep="\\t",
            file=manifest
        )

with open(stats_path, "w") as stats:
    print(
        "sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembly_assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "binning_id",
        "binner",
        "bin_count",
        "total_bin_bp",
        "largest_bin_bp",
        "bins_dir",
        "log_file",
        "report_file",
        sep="\\t",
        file=stats
    )
    print(
        sample_id,
        safe_sample_id,
        assembly_sample_id,
        assembly_assembler,
        assembly_mode,
        rarefaction_label,
        assembly_strategy,
        binning_id,
        binner,
        len(bin_files),
        total_bp,
        largest,
        published_bin_dir,
        published_log,
        published_report,
        sep="\\t",
        file=stats
    )
PY

    echo "MetaBAT2 binning finished: \$(date)" >> "\$LOG_FILE"
    """
}


process RUN_MAXBIN2 {
    tag { binning_id }
    stageInMode 'symlink'

    publishDir { "${params.outdir}/bins/maxbin2/${binning_id}" },
        mode: params.publish_bins_mode,
        pattern: "bins/**",
        saveAs: { filename -> filename.replaceFirst(/^bins\//, '') }

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.maxbin2.log"

    publishDir "${params.outdir}/summary/per_binner_stats",
        mode: 'copy',
        pattern: "*.maxbin2.binning_stats.tsv"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.binning_threads as int
    }

    input:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(filtered_assembly_fasta),
          path(sorted_bam),
          path(sorted_bam_bai),
          path(depth_file),
          path(maxbin_abundance),
          path(tools_status)

    output:
    path "bins", emit: bin_dir
    path "${binning_id}.maxbin2.binning_manifest_record.tsv", emit: manifest_record
    path "${binning_id}.maxbin2.binning_stats.tsv", emit: stats_file
    path "${binning_id}.maxbin2.log", emit: log_file

    script:
    """
    set -euo pipefail

    TASK_DIR="\$(pwd -P)"
    ORIG_MAXBIN2_ENV="\$(grep '^MAXBIN2_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"
    MAXBIN2_ENV="\$ORIG_MAXBIN2_ENV"

    BINNER="maxbin2"
    BIN_DIR="bins"
    LOG_FILE="\$TASK_DIR/${binning_id}.maxbin2.log"

    echo "MaxBin2 binning started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: ${sample_id}" >> "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "Abundance file: ${maxbin_abundance}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "Task directory: \$TASK_DIR" >> "\$LOG_FILE"
    echo "Original MAXBIN2_ENV: \$ORIG_MAXBIN2_ENV" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    if [[ ! -s "${filtered_assembly_fasta}" ]]; then
        echo "ERROR: Filtered assembly FASTA is missing or empty: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    if [[ ! -s "${maxbin_abundance}" ]]; then
        echo "ERROR: MaxBin2 abundance file is missing or empty: ${maxbin_abundance}" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    TMP_PARENT="\${TMPDIR:-/tmp}"

    if [[ -z "\$TMP_PARENT" || "\$TMP_PARENT" == *" "* ]]; then
        TMP_PARENT="/tmp"
    fi

    mkdir -p "\$TMP_PARENT"

    SANDBOX_RAW="\$(mktemp -d "\$TMP_PARENT/samwise_maxbin2.XXXXXX")"
    SANDBOX="\$(cd "\$SANDBOX_RAW" && pwd -P)"

    if [[ "\$SANDBOX" == *" "* ]]; then
        echo "ERROR: MaxBin2 sandbox path still contains spaces: \$SANDBOX" >> "\$LOG_FILE"
        echo "Set TMPDIR to a no-space path, for example:" >> "\$LOG_FILE"
        echo "  export TMPDIR=/tmp" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    trap 'rm -rf "\$SANDBOX"' EXIT

    RUN_DIR="\$SANDBOX/run"
    mkdir -p "\$RUN_DIR"

    echo "Using MaxBin2 sandbox: \$SANDBOX" >> "\$LOG_FILE"
    echo "MaxBin2 run directory: \$RUN_DIR" >> "\$LOG_FILE"

    NEED_ENV_COPY="false"

    if [[ -n "\$MAXBIN2_ENV" && "\$MAXBIN2_ENV" != "SYSTEM" && "\$MAXBIN2_ENV" != "NOT_USED" ]]; then
        if [[ "\$MAXBIN2_ENV" == *" "* ]]; then
            NEED_ENV_COPY="true"
        fi
    fi

    if [[ "\$NEED_ENV_COPY" == "true" ]]; then
        echo "MAXBIN2_ENV contains spaces." >> "\$LOG_FILE"
        echo "Creating no-space sandbox copy of MaxBin2 environment." >> "\$LOG_FILE"

        SANDBOX_MAXBIN2_ENV="\$SANDBOX/maxbin2_env"
        mkdir -p "\$SANDBOX_MAXBIN2_ENV"

        echo "Copying MaxBin2 env while preserving symlinks:" >> "\$LOG_FILE"
        echo "  from: \$MAXBIN2_ENV" >> "\$LOG_FILE"
        echo "  to:   \$SANDBOX_MAXBIN2_ENV" >> "\$LOG_FILE"

        # IMPORTANT:
        # Do NOT use tar -h here.
        # Preserving symlinks is important for conda-installed MaxBin2,
        # because run_MaxBin.pl may be a symlink into a share/ directory
        # where helper scripts such as _getmarker.pl live.
        if ! (cd "\$MAXBIN2_ENV" && tar -cf - .) | (cd "\$SANDBOX_MAXBIN2_ENV" && tar -xpf -); then
            echo "ERROR: Failed to copy MaxBin2 environment into sandbox." >> "\$LOG_FILE"
            cat "\$LOG_FILE" >&2
            exit 1
        fi

        MAXBIN2_ENV="\$SANDBOX_MAXBIN2_ENV"
        echo "Sandbox MAXBIN2_ENV: \$MAXBIN2_ENV" >> "\$LOG_FILE"
    else
        echo "MAXBIN2_ENV does not require sandbox copy." >> "\$LOG_FILE"
    fi

    if [[ -n "\$MAXBIN2_ENV" && "\$MAXBIN2_ENV" != "SYSTEM" && "\$MAXBIN2_ENV" != "NOT_USED" ]]; then
        export PATH="\$MAXBIN2_ENV/bin:\$PATH"
    fi

    echo "Effective MAXBIN2_ENV: \$MAXBIN2_ENV" >> "\$LOG_FILE"
    echo "PATH: \$PATH" >> "\$LOG_FILE"

    repair_maxbin2_helpers() {
        local env_dir="\$1"
        local original_env_dir="\$2"
        local helper_list="\$SANDBOX/maxbin2_helper_files.list"

        if [[ "\$env_dir" == "SYSTEM" || "\$env_dir" == "NOT_USED" || ! -d "\$env_dir" ]]; then
            return 0
        fi

        mkdir -p "\$env_dir/bin"
        : > "\$helper_list"

        echo "Repairing/checking MaxBin2 helper scripts." >> "\$LOG_FILE"

        if [[ -d "\$env_dir" ]]; then
            find -L "\$env_dir" -type f -name '_*.pl' -print >> "\$helper_list" 2>/dev/null || true
        fi

        if [[ -n "\$original_env_dir" && "\$original_env_dir" != "SYSTEM" && "\$original_env_dir" != "NOT_USED" && -d "\$original_env_dir" ]]; then
            find -L "\$original_env_dir" -type f -name '_*.pl' -print >> "\$helper_list" 2>/dev/null || true
        fi

        echo "Discovered MaxBin2 helper candidates:" >> "\$LOG_FILE"
        cat "\$helper_list" >> "\$LOG_FILE" 2>/dev/null || true

        while IFS= read -r helper; do
            [[ -z "\$helper" ]] && continue
            base="\$(basename "\$helper")"

            if [[ ! -e "\$env_dir/bin/\$base" ]]; then
                echo "Installing helper into sandbox env bin: \$base" >> "\$LOG_FILE"
                cp -L "\$helper" "\$env_dir/bin/\$base"
                chmod +x "\$env_dir/bin/\$base" 2>/dev/null || true
            fi
        done < "\$helper_list"

        if [[ ! -s "\$env_dir/bin/_getmarker.pl" ]]; then
            echo "WARNING: _getmarker.pl is not present in \$env_dir/bin after repair." >> "\$LOG_FILE"
            echo "This may still be okay if run_MaxBin.pl is executed from a share/ directory containing helpers." >> "\$LOG_FILE"
        else
            echo "Confirmed helper exists: \$env_dir/bin/_getmarker.pl" >> "\$LOG_FILE"
        fi

        echo "MaxBin2 helper repair/check complete." >> "\$LOG_FILE"
        ls -lah "\$env_dir/bin"/_*.pl >> "\$LOG_FILE" 2>&1 || true
    }

    repair_maxbin2_helpers "\$MAXBIN2_ENV" "\$ORIG_MAXBIN2_ENV"

    RUN_MAXBIN=""
    RUN_CANDIDATES="\$SANDBOX/run_maxbin_candidates.list"
    : > "\$RUN_CANDIDATES"

    if [[ "\$MAXBIN2_ENV" != "SYSTEM" && "\$MAXBIN2_ENV" != "NOT_USED" && -d "\$MAXBIN2_ENV" ]]; then
        find -L "\$MAXBIN2_ENV" -type f -name 'run_MaxBin.pl' -print >> "\$RUN_CANDIDATES" 2>/dev/null || true
    fi

    if command -v run_MaxBin.pl >/dev/null 2>&1; then
        command -v run_MaxBin.pl >> "\$RUN_CANDIDATES" 2>/dev/null || true
    fi

    echo "run_MaxBin.pl candidates:" >> "\$LOG_FILE"
    cat "\$RUN_CANDIDATES" >> "\$LOG_FILE" 2>/dev/null || true

    # Prefer a run_MaxBin.pl whose directory also contains _getmarker.pl.
    while IFS= read -r candidate; do
        [[ -z "\$candidate" ]] && continue
        candidate_dir="\$(dirname "\$candidate")"

        if [[ -s "\$candidate_dir/_getmarker.pl" ]]; then
            RUN_MAXBIN="\$candidate"
            echo "Selected run_MaxBin.pl with local helpers: \$RUN_MAXBIN" >> "\$LOG_FILE"
            break
        fi
    done < "\$RUN_CANDIDATES"

    # Fallback to env/bin/run_MaxBin.pl if needed.
    if [[ -z "\$RUN_MAXBIN" ]]; then
        if [[ "\$MAXBIN2_ENV" != "SYSTEM" && "\$MAXBIN2_ENV" != "NOT_USED" && -s "\$MAXBIN2_ENV/bin/run_MaxBin.pl" ]]; then
            RUN_MAXBIN="\$MAXBIN2_ENV/bin/run_MaxBin.pl"
            echo "Selected fallback env/bin/run_MaxBin.pl: \$RUN_MAXBIN" >> "\$LOG_FILE"
        else
            RUN_MAXBIN="\$(command -v run_MaxBin.pl || true)"
            echo "Selected fallback PATH run_MaxBin.pl: \$RUN_MAXBIN" >> "\$LOG_FILE"
        fi
    fi

    if [[ -z "\$RUN_MAXBIN" || ! -s "\$RUN_MAXBIN" ]]; then
        echo "ERROR: run_MaxBin.pl is not available after tool setup." >> "\$LOG_FILE"
        echo "MAXBIN2_ENV=\$MAXBIN2_ENV" >> "\$LOG_FILE"
        cat "${tools_status}" >> "\$LOG_FILE" 2>&1 || true
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    if [[ "\$MAXBIN2_ENV" != "SYSTEM" && "\$MAXBIN2_ENV" != "NOT_USED" && -x "\$MAXBIN2_ENV/bin/perl" ]]; then
        PERL_EXE="\$MAXBIN2_ENV/bin/perl"
    else
        PERL_EXE="\$(command -v perl || true)"
    fi

    if [[ -z "\$PERL_EXE" || ! -x "\$PERL_EXE" ]]; then
        echo "ERROR: perl is not available for MaxBin2." >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    echo "Final run_MaxBin.pl path: \$RUN_MAXBIN" >> "\$LOG_FILE"
    echo "Final run_MaxBin.pl directory: \$(dirname "\$RUN_MAXBIN")" >> "\$LOG_FILE"
    echo "Perl executable: \$PERL_EXE" >> "\$LOG_FILE"

    SOURCE_CONTIG="\$TASK_DIR/${filtered_assembly_fasta}"
    SOURCE_ABUNDANCE="\$TASK_DIR/${maxbin_abundance}"

    echo "Creating sandbox symlinks instead of copying large input files." >> "\$LOG_FILE"
    echo "Source contigs: \$SOURCE_CONTIG" >> "\$LOG_FILE"
    echo "Source abundance: \$SOURCE_ABUNDANCE" >> "\$LOG_FILE"

    ln -s "\$SOURCE_CONTIG" "\$RUN_DIR/maxbin2_contigs.fa"
    ln -s "\$SOURCE_ABUNDANCE" "\$RUN_DIR/maxbin2_abundance.tsv"

    if [[ ! -s "\$RUN_DIR/maxbin2_contigs.fa" ]]; then
        echo "ERROR: Sandbox contig symlink is missing, broken, or points to an empty file." >> "\$LOG_FILE"
        echo "Symlink: \$RUN_DIR/maxbin2_contigs.fa" >> "\$LOG_FILE"
        ls -lah "\$RUN_DIR" >> "\$LOG_FILE" 2>&1 || true
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    if [[ ! -s "\$RUN_DIR/maxbin2_abundance.tsv" ]]; then
        echo "ERROR: Sandbox abundance symlink is missing, broken, or points to an empty file." >> "\$LOG_FILE"
        echo "Symlink: \$RUN_DIR/maxbin2_abundance.tsv" >> "\$LOG_FILE"
        ls -lah "\$RUN_DIR" >> "\$LOG_FILE" 2>&1 || true
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    echo "Sandbox input symlinks:" >> "\$LOG_FILE"
    ls -lah "\$RUN_DIR" >> "\$LOG_FILE" 2>&1 || true

    rm -rf "\$TASK_DIR/bins"
    mkdir -p "\$TASK_DIR/bins"

    ln -s "\$TASK_DIR/bins" "\$RUN_DIR/bins"

    if [[ ! -d "\$RUN_DIR/bins" ]]; then
        echo "ERROR: Sandbox bins symlink is missing or broken." >> "\$LOG_FILE"
        echo "Symlink: \$RUN_DIR/bins" >> "\$LOG_FILE"
        ls -lah "\$RUN_DIR" >> "\$LOG_FILE" 2>&1 || true
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    echo "Sandbox bins symlink:" >> "\$LOG_FILE"
    ls -lah "\$RUN_DIR/bins" >> "\$LOG_FILE" 2>&1 || true

    cd "\$RUN_DIR"

    echo "Running MaxBin2 inside sandbox." >> "\$LOG_FILE"
    echo "Sandbox pwd: \$(pwd -P)" >> "\$LOG_FILE"

    MAXBIN_CMD=( "\$PERL_EXE" "\$RUN_MAXBIN" )
    MAXBIN_CMD+=( -contig maxbin2_contigs.fa )
    MAXBIN_CMD+=( -abund maxbin2_abundance.tsv )
    MAXBIN_CMD+=( -out "bins/${binning_id}.maxbin2_bin" )
    MAXBIN_CMD+=( -thread ${task.cpus} )

    if [[ -n "${params.maxbin2_extra_args}" ]]; then
        EXTRA_ARGS=( ${params.maxbin2_extra_args} )
        MAXBIN_CMD+=( "\${EXTRA_ARGS[@]}" )
    fi

    printf 'MaxBin2 command:' >> "\$LOG_FILE"
    printf ' %q' "\${MAXBIN_CMD[@]}" >> "\$LOG_FILE"
    printf '\\n' >> "\$LOG_FILE"

    set +e
    "\${MAXBIN_CMD[@]}" >> "\$LOG_FILE" 2>&1
    STATUS="\$?"
    set -e

    if [[ "\$STATUS" -ne 0 ]]; then
        echo "ERROR: MaxBin2 failed with exit status \$STATUS" >> "\$LOG_FILE"
        echo "Last 180 lines of MaxBin2 log:" >&2
        tail -n 180 "\$LOG_FILE" >&2 || true
        exit "\$STATUS"
    fi

    echo "MaxBin2 command completed." >> "\$LOG_FILE"
    echo "Sandbox bins directory contents:" >> "\$LOG_FILE"
    find -L bins -type f -print >> "\$LOG_FILE" 2>&1 || true

    cd "\$TASK_DIR"

    echo "MaxBin2 wrote bins directly into task directory via sandbox symlink." >> "\$LOG_FILE"
    echo "Task bins directory contents:" >> "\$LOG_FILE"
    find bins -type f -print >> "\$LOG_FILE" 2>&1 || true

    python3 - \\
        "bins" \\
        "${binning_id}.maxbin2.binning_manifest_record.tsv" \\
        "${binning_id}.maxbin2.binning_stats.tsv" \\
        "${sample_id}" \\
        "${safe_sample_id}" \\
        "${assembly_sample_id}" \\
        "${assembly_assembler}" \\
        "${assembly_mode}" \\
        "${rarefaction_label}" \\
        "${assembly_strategy}" \\
        "${binning_id}" \\
        "\$BINNER" \\
        "${params.outdir}/bins/maxbin2/${binning_id}" \\
        "${params.outdir}/logs/${binning_id}.maxbin2.log" \\
        "" <<'PY'
import gzip
import sys
from pathlib import Path

(
    bin_dir,
    manifest_path,
    stats_path,
    sample_id,
    safe_sample_id,
    assembly_sample_id,
    assembly_assembler,
    assembly_mode,
    rarefaction_label,
    assembly_strategy,
    binning_id,
    binner,
    published_bin_dir,
    published_log,
    published_report
) = sys.argv[1:]

bin_dir = Path(bin_dir)

def fasta_len(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total

fasta_exts = {".fa", ".fasta", ".fna", ".fas"}

bin_files = []
for p in sorted(bin_dir.iterdir()):
    if p.is_file():
        suffixes = p.suffixes
        if p.suffix in fasta_exts or (
            len(suffixes) >= 2 and suffixes[-1] == ".gz" and suffixes[-2] in fasta_exts
        ):
            bin_files.append(p)

total_bp = 0
largest = 0

with open(manifest_path, "w") as manifest:
    for idx, p in enumerate(bin_files, start=1):
        bp = fasta_len(p)
        total_bp += bp
        largest = max(largest, bp)
        bin_id = f"{binning_id}_{binner}_bin_{idx:03d}"
        published_path = str(Path(published_bin_dir) / p.name)

        print(
            sample_id,
            safe_sample_id,
            assembly_sample_id,
            assembly_assembler,
            assembly_mode,
            rarefaction_label,
            assembly_strategy,
            binning_id,
            binner,
            bin_id,
            bp,
            published_path,
            sep="\\t",
            file=manifest
        )

with open(stats_path, "w") as stats:
    print(
        "sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembly_assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "binning_id",
        "binner",
        "bin_count",
        "total_bin_bp",
        "largest_bin_bp",
        "bins_dir",
        "log_file",
        "report_file",
        sep="\\t",
        file=stats
    )

    print(
        sample_id,
        safe_sample_id,
        assembly_sample_id,
        assembly_assembler,
        assembly_mode,
        rarefaction_label,
        assembly_strategy,
        binning_id,
        binner,
        len(bin_files),
        total_bp,
        largest,
        published_bin_dir,
        published_log,
        published_report,
        sep="\\t",
        file=stats
    )
PY

    echo "MaxBin2 binning finished: \$(date)" >> "\$LOG_FILE"
    """
}

process RUN_QUICKBIN {
    tag { binning_id }
    stageInMode 'symlink'

    publishDir { "${params.outdir}/bins/quickbin/${binning_id}" },
        mode: params.publish_bins_mode,
        pattern: "bins/**",
        saveAs: { filename -> filename.replaceFirst(/^bins\//, '') }

    publishDir "${params.outdir}/coverage/quickbin_cov",
        mode: 'copy',
        pattern: "*.quickbin_cov.txt"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.quickbin.*"

    publishDir "${params.outdir}/summary/per_binner_stats",
        mode: 'copy',
        pattern: "*.quickbin.binning_stats.tsv"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.binning_threads as int
    }

    input:
    tuple val(sample_id),
          val(safe_sample_id),
          val(assembly_sample_id),
          val(assembly_assembler),
          val(assembly_mode),
          val(rarefaction_label),
          val(assembly_strategy),
          val(binning_id),
          path(filtered_assembly_fasta),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          path(sorted_bam),
          path(sorted_bam_bai),
          path(tools_status)

    output:
    path "bins", emit: bin_dir
    path "${binning_id}.quickbin_cov.txt", emit: cov_file
    path "${binning_id}.quickbin.report.tsv", emit: report_file
    path "${binning_id}.quickbin.binning_manifest_record.tsv", emit: manifest_record
    path "${binning_id}.quickbin.binning_stats.tsv", emit: stats_file
    path "${binning_id}.quickbin.log", emit: log_file

    script:
    def quickbin_xmx_arg = ""
    if( params.quickbin_xmx ) {
        def x = params.quickbin_xmx.toString()
        quickbin_xmx_arg = x.startsWith("-Xmx") ? x : "-Xmx${x}"
    }

    def quickbin_xmx_log = quickbin_xmx_arg ? quickbin_xmx_arg : "NONE"

    """
    set -euo pipefail

    QUICKBIN_ENV="\$(grep '^QUICKBIN_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$QUICKBIN_ENV" && "\$QUICKBIN_ENV" != "SYSTEM" && "\$QUICKBIN_ENV" != "NOT_USED" ]]; then
        export PATH="\$QUICKBIN_ENV/bin:\$PATH"
    fi

    if ! command -v java >/dev/null 2>&1; then
        echo "ERROR: java is not available after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    BBTOOLS_JAR=""
    if [[ -n "\$QUICKBIN_ENV" && "\$QUICKBIN_ENV" != "SYSTEM" && "\$QUICKBIN_ENV" != "NOT_USED" ]]; then
        BBTOOLS_JAR="\$(find "\$QUICKBIN_ENV" -type f -name 'bbtools.jar' 2>/dev/null | head -n 1 || true)"
    fi

    if [[ -z "\$BBTOOLS_JAR" ]]; then
        echo "ERROR: Could not find bbtools.jar under QUICKBIN_ENV." >&2
        echo "QUICKBIN_ENV=\$QUICKBIN_ENV" >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    BINNER="quickbin"
    BIN_DIR="bins"
    LOG_FILE="${binning_id}.quickbin.log"
    REPORT_FILE="${binning_id}.quickbin.report.tsv"
    COV_FILE="${binning_id}.quickbin_cov.txt"

    mkdir -p "\$BIN_DIR"

    echo "QuickBin binning started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: ${sample_id}" >> "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "BAM: ${sorted_bam}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "QUICKBIN_ENV: \$QUICKBIN_ENV" >> "\$LOG_FILE"
    echo "BBTOOLS_JAR: \$BBTOOLS_JAR" >> "\$LOG_FILE"
    echo "Requested Java memory argument: ${quickbin_xmx_log}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    if [[ ! -s "${filtered_assembly_fasta}" ]]; then
        echo "ERROR: Filtered assembly FASTA is missing or empty: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    if [[ ! -s "${sorted_bam}" ]]; then
        echo "ERROR: Sorted BAM is missing or empty: ${sorted_bam}" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    QB_ARGS=()
    QB_ARGS+=( "mincluster=${params.quickbin_mincluster}" )
    QB_ARGS+=( "mincontig=${params.quickbin_mincontig}" )
    QB_ARGS+=( "minseed=${params.quickbin_minseed}" )
    QB_ARGS+=( "threads=${task.cpus}" )
    QB_ARGS+=( "gzip=${params.quickbin_gzip.toString().toLowerCase()}" )
    QB_ARGS+=( "clade=${params.quickbin_clade.toString().toLowerCase()}" )
    QB_ARGS+=( "sketch=${params.quickbin_sketch.toString().toLowerCase()}" )
    QB_ARGS+=( "server=${params.quickbin_server.toString().toLowerCase()}" )

    if [[ "${params.quickbin_chaff.toString().toLowerCase()}" == "true" ]]; then
        QB_ARGS+=( "chaff" )
    fi

    if [[ -n "${params.quickbin_stringency}" ]]; then
        QB_ARGS+=( "${params.quickbin_stringency}" )
    fi

    if [[ -n "${params.quickbin_extra_args}" ]]; then
        EXTRA_ARGS=( ${params.quickbin_extra_args} )
        QB_ARGS+=( "\${EXTRA_ARGS[@]}" )
    fi

    if [[ "${params.quickbin_use_positional_bam.toString().toLowerCase()}" == "true" ]]; then
        READS_ARG="${sorted_bam}"
    else
        READS_ARG="reads=${sorted_bam}"
    fi

    echo "Using direct Java QuickBin call instead of quickbin.sh." >> "\$LOG_FILE"
    echo "This avoids BBTools wrapper classpath breakage when paths contain spaces." >> "\$LOG_FILE"
    echo "READS_ARG: \$READS_ARG" >> "\$LOG_FILE"

    printf 'QB_ARGS:' >> "\$LOG_FILE"
    printf ' %q' "\${QB_ARGS[@]}" >> "\$LOG_FILE"
    printf '\\n' >> "\$LOG_FILE"

    JAVA_CMD=( java -ea )

    if [[ -n "${quickbin_xmx_arg}" ]]; then
        JAVA_CMD+=( "${quickbin_xmx_arg}" )
    fi

    JAVA_CMD+=( -cp "\$BBTOOLS_JAR" )
    JAVA_CMD+=( bin.QuickBin )
    JAVA_CMD+=( "in=${filtered_assembly_fasta}" )
    JAVA_CMD+=( "\$READS_ARG" )
    JAVA_CMD+=( "out=\$BIN_DIR" )
    JAVA_CMD+=( "covout=\$COV_FILE" )
    JAVA_CMD+=( "report=\$REPORT_FILE" )
    JAVA_CMD+=( "\${QB_ARGS[@]}" )

    printf 'Java command:' >> "\$LOG_FILE"
    printf ' %q' "\${JAVA_CMD[@]}" >> "\$LOG_FILE"
    printf '\\n' >> "\$LOG_FILE"

    "\${JAVA_CMD[@]}" >> "\$LOG_FILE" 2>&1

    if [[ ! -d "\$BIN_DIR" ]]; then
        echo "ERROR: QuickBin did not create bin directory: \$BIN_DIR" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    if [[ ! -s "\$COV_FILE" ]]; then
        echo "WARNING: QuickBin did not create a non-empty coverage file: \$COV_FILE" >> "\$LOG_FILE"
    fi

    if [[ ! -s "\$REPORT_FILE" ]]; then
        echo "WARNING: QuickBin did not create a non-empty report file: \$REPORT_FILE" >> "\$LOG_FILE"
    fi

    python3 - \\
        "\$BIN_DIR" \\
        "${binning_id}.quickbin.binning_manifest_record.tsv" \\
        "${binning_id}.quickbin.binning_stats.tsv" \\
        "${sample_id}" \\
        "${safe_sample_id}" \\
        "${assembly_sample_id}" \\
        "${assembly_assembler}" \\
        "${assembly_mode}" \\
        "${rarefaction_label}" \\
        "${assembly_strategy}" \\
        "${binning_id}" \\
        "\$BINNER" \\
        "${params.outdir}/bins/quickbin/${binning_id}" \\
        "${params.outdir}/logs/${binning_id}.quickbin.log" \\
        "${params.outdir}/logs/${binning_id}.quickbin.report.tsv" <<'PY'
import gzip
import sys
from pathlib import Path

(
    bin_dir,
    manifest_path,
    stats_path,
    sample_id,
    safe_sample_id,
    assembly_sample_id,
    assembly_assembler,
    assembly_mode,
    rarefaction_label,
    assembly_strategy,
    binning_id,
    binner,
    published_bin_dir,
    published_log,
    published_report
) = sys.argv[1:]

bin_dir = Path(bin_dir)

def fasta_len(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    total = 0
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total

fasta_exts = {".fa", ".fasta", ".fna", ".fas"}

bin_files = []
for p in sorted(bin_dir.iterdir()):
    if p.is_file():
        suffixes = p.suffixes
        if p.suffix in fasta_exts or (
            len(suffixes) >= 2 and suffixes[-1] == ".gz" and suffixes[-2] in fasta_exts
        ):
            bin_files.append(p)

total_bp = 0
largest = 0

with open(manifest_path, "w") as manifest:
    for idx, p in enumerate(bin_files, start=1):
        bp = fasta_len(p)
        total_bp += bp
        largest = max(largest, bp)
        bin_id = f"{binning_id}_{binner}_bin_{idx:03d}"
        published_path = str(Path(published_bin_dir) / p.name)

        print(
            sample_id,
            safe_sample_id,
            assembly_sample_id,
            assembly_assembler,
            assembly_mode,
            rarefaction_label,
            assembly_strategy,
            binning_id,
            binner,
            bin_id,
            bp,
            published_path,
            sep="\\t",
            file=manifest
        )

with open(stats_path, "w") as stats:
    print(
        "sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembly_assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "binning_id",
        "binner",
        "bin_count",
        "total_bin_bp",
        "largest_bin_bp",
        "bins_dir",
        "log_file",
        "report_file",
        sep="\\t",
        file=stats
    )

    print(
        sample_id,
        safe_sample_id,
        assembly_sample_id,
        assembly_assembler,
        assembly_mode,
        rarefaction_label,
        assembly_strategy,
        binning_id,
        binner,
        len(bin_files),
        total_bp,
        largest,
        published_bin_dir,
        published_log,
        published_report,
        sep="\\t",
        file=stats
    )
PY

    echo "QuickBin binning finished: \$(date)" >> "\$LOG_FILE"
    """
}

process WRITE_FILTERED_ASSEMBLY_STATS_SUMMARY {

    tag "write_filtered_assembly_stats_summary"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "filtered_assembly_stats_summary.tsv"

    input:
    path stats_files

    output:
    path "filtered_assembly_stats_summary.tsv", emit: summary

    script:
    def stats_file_list = stats_files.collect { stats -> stats.name }.join(' ')

    """
    set -euo pipefail

    if [[ -z "${stats_file_list}" ]]; then
        echo "ERROR: No filtered assembly stats files were received." >&2
        exit 1
    fi

    first=1
    : > filtered_assembly_stats_summary.tsv

    for f in ${stats_file_list}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" >> filtered_assembly_stats_summary.tsv
            first=0
        else
            tail -n +2 "\$f" >> filtered_assembly_stats_summary.tsv
        fi
    done
    """
}


process WRITE_BINNING_SUMMARIES {

    tag "write_binning_summaries"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "binning_manifest.tsv"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "binning_stats_summary.tsv"

    input:
    path manifest_records
    path stats_files

    output:
    path "binning_manifest.tsv", emit: manifest
    path "binning_stats_summary.tsv", emit: stats_summary

    script:
    def manifest_files = manifest_records.collect { record -> record.name }.join(' ')
    def stats_file_list = stats_files.collect { stats -> stats.name }.join(' ')

    """
    set -euo pipefail

    if [[ -z "${manifest_files}" ]]; then
        echo "ERROR: No binning manifest records were received." >&2
        exit 1
    fi

    if [[ -z "${stats_file_list}" ]]; then
        echo "ERROR: No binning stats files were received." >&2
        exit 1
    fi

    printf 'sample_id\\tsafe_sample_id\\tassembly_sample_id\\tassembly_assembler\\tassembly_mode\\trarefaction_label\\tassembly_strategy\\tbinning_id\\tbinner\\tbin_id\\tbin_bp\\tbin_fasta\\n' > binning_manifest.tsv

    for f in ${manifest_files}; do
        cat "\$f" >> binning_manifest.tsv
    done

    first=1
    : > binning_stats_summary.tsv

    for f in ${stats_file_list}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" >> binning_stats_summary.tsv
            first=0
        else
            tail -n +2 "\$f" >> binning_stats_summary.tsv
        fi
    done
    """
}