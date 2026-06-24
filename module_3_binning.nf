#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Module 3: MAG binning from Module 2 assemblies and Module 1 trimmed reads.
 */

params.working_dir = null
params.input_assembly_manifest = null
params.input_trimmed_manifest = null
params.input_coassembly_assembly_manifest = null
params.input_coassembly_trimmed_manifest = null
params.include_module2 = true
params.include_module2b = true
params.output_dir = null
params.metabat2 = false
params.quickbin = false
params.maxbin2 = false
params.auto_install = true
params.tool_env_dir = null
params.threads = null
params.mapping_threads = 4
params.binning_threads = 4
params.seqkit_version = "2.8.2"
params.bbmap_version = "39.81"
params.samtools_version = "1.23.1"
params.metabat2_version = "2.18"
params.maxbin2_version = "2.2.7"
params.min_scaffold_length = 2500
params.bbmap_minid = 0.90
params.bbmap_maxindel = 10
params.bbmap_ambig = "random"
params.bbmap_mateqtag = true
params.bbmap_extra_args = ""
params.bbmap_xmx = "4g"
params.metabat2_min_contig = 2500
params.metabat2_extra_args = ""
params.maxbin2_extra_args = ""
params.quickbin_mincluster = "50k"
params.quickbin_mincontig = 2500
params.quickbin_minseed = 2500
params.quickbin_stringency = "normal"
params.quickbin_gzip = false
params.quickbin_chaff = false
params.quickbin_clade = false
params.quickbin_sketch = false
params.quickbin_server = false
params.quickbin_xmx = null
params.quickbin_extra_args = ""
params.quickbin_use_positional_bam = false
params.publish_filtered_assemblies_mode = "symlink"
params.publish_bam_mode = "symlink"
params.publish_bins_mode = "symlink"
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.module1_outdir = "${params.results_dir}/module_1_readtrimming"
params.module2_outdir = "${params.results_dir}/module_2_readassembly"
params.module2b_outdir = "${params.results_dir}/module_2b_coassembly"
params.outdir = "${params.results_dir}/module_3_binning"

def absOrEmpty(value) {
    def s = value == null ? "" : value.toString().trim()
    if (!s || s == "null" || s == "NA") {
        return ""
    }
    return java.nio.file.Paths.get(s).toAbsolutePath().normalize().toString()
}

def resolveForDiscoveryPath(value, launchDir) {
    def p = java.nio.file.Paths.get(value.toString())

    if (!p.isAbsolute()) {
        p = java.nio.file.Paths.get(launchDir.toString()).resolve(p)
    }

    return p.toAbsolutePath().normalize().toString()
}

def existsForDiscoveryPath(value, launchDir) {
    return java.nio.file.Files.exists(
        java.nio.file.Paths.get(
            resolveForDiscoveryPath(value, launchDir)
        )
    )
}

workflow {
    def use_metabat2 = params.metabat2.toString().toBoolean()
    def use_quickbin = params.quickbin.toString().toBoolean()
    def use_maxbin2 = params.maxbin2.toString().toBoolean()

    if (!use_metabat2 && !use_quickbin && !use_maxbin2) {
        error(
            """
        No binner selected.
        Please specify at least one of:
          --metabat2
          --quickbin
          --maxbin2
        Example:
          nextflow run module_3_binning.nf --working_dir ./output_samwise --threads 6 --metabat2 --quickbin --maxbin2
        """.stripIndent()
        )
    }

    def include_module2 = params.include_module2.toString().toBoolean()
    def include_module2b = params.include_module2b.toString().toBoolean()
    def launch_dir = workflow.launchDir.toString()

    /*
     * Regular Module 2 inputs.
     */
    def module2_assembly_manifest_file = params.input_assembly_manifest ?: "${params.module2_outdir}/summary/assembly_manifest.tsv"
    def module2_trimmed_manifest_file = params.input_trimmed_manifest ?: "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    /*
     * Module 2b coassembly inputs.
     */
    def module2b_assembly_manifest_file = params.input_coassembly_assembly_manifest ?: "${params.module2b_outdir}/summary/assembly_manifest.tsv"
    def module2b_trimmed_manifest_file = params.input_coassembly_trimmed_manifest ?: "${params.module2b_outdir}/summary/coassembly_trimmed_manifest.tsv"

    def assembly_manifest_files = []
    def trimmed_manifest_files = []
    def module2_explicit = params.input_assembly_manifest != null
    def module2_exists = existsForDiscoveryPath(module2_assembly_manifest_file, launch_dir)

    if (module2_explicit && !module2_exists) {
        error(
            """
        Explicit Module 2 assembly manifest was supplied but does not exist:

          --input_assembly_manifest ${module2_assembly_manifest_file}
        """.stripIndent()
        )
    }

    if (module2_explicit || (include_module2 && module2_exists)) {
        if (!existsForDiscoveryPath(module2_trimmed_manifest_file, launch_dir)) {
            error(
                """
        Regular Module 2 assembly manifest was selected, but the corresponding
        trimmed-read manifest does not exist.

        Assembly manifest:
          ${module2_assembly_manifest_file}

        Trimmed manifest:
          ${module2_trimmed_manifest_file}

        If you only want to bin Module 2b coassemblies, run with:

          --include_module2 false
        """.stripIndent()
            )
        }

        assembly_manifest_files << resolveForDiscoveryPath(module2_assembly_manifest_file, launch_dir)
        trimmed_manifest_files << resolveForDiscoveryPath(module2_trimmed_manifest_file, launch_dir)
    }

    def module2b_explicit = (params.input_coassembly_assembly_manifest != null || params.input_coassembly_trimmed_manifest != null)

    def module2b_exists = existsForDiscoveryPath(module2b_assembly_manifest_file, launch_dir)

    if (params.input_coassembly_assembly_manifest != null && !module2b_exists) {
        error(
            """
        Explicit Module 2b coassembly assembly manifest was supplied but does not exist:

          --input_coassembly_assembly_manifest ${module2b_assembly_manifest_file}
        """.stripIndent()
        )
    }

    if (module2b_explicit || (include_module2b && module2b_exists)) {
        if (!existsForDiscoveryPath(module2b_assembly_manifest_file, launch_dir)) {
            error(
                """
        Module 2b coassembly discovery was requested, but the coassembly
        assembly manifest does not exist:

          ${module2b_assembly_manifest_file}
        """.stripIndent()
            )
        }

        if (!existsForDiscoveryPath(module2b_trimmed_manifest_file, launch_dir)) {
            error(
                """
        Module 2b coassembly assembly manifest was selected, but the corresponding
        coassembly trimmed-read manifest does not exist.

        Coassembly assembly manifest:
          ${module2b_assembly_manifest_file}

        Coassembly trimmed manifest:
          ${module2b_trimmed_manifest_file}

        Expected Module 2b output:
          ${params.module2b_outdir}/summary/coassembly_trimmed_manifest.tsv
        """.stripIndent()
            )
        }

        assembly_manifest_files << resolveForDiscoveryPath(module2b_assembly_manifest_file, launch_dir)
        trimmed_manifest_files << resolveForDiscoveryPath(module2b_trimmed_manifest_file, launch_dir)
    }

    assembly_manifest_files = assembly_manifest_files.unique()
    trimmed_manifest_files = trimmed_manifest_files.unique()

    if (assembly_manifest_files.isEmpty()) {
        error(
            """
        No assembly manifests were found for Module 3.

        Checked regular Module 2:
          ${module2_assembly_manifest_file}

        Checked Module 2b coassembly:
          ${module2b_assembly_manifest_file}

        At least one assembly manifest is required.

        Options:
          1. Run Module 2 first.
          2. Run Module 2b first.
          3. Provide explicit manifest paths.

        Examples:

          --input_assembly_manifest path/to/module2/assembly_manifest.tsv
          --input_trimmed_manifest path/to/module1/trimmed_manifest.tsv

        or:

          --input_coassembly_assembly_manifest path/to/module2b/assembly_manifest.tsv
          --input_coassembly_trimmed_manifest path/to/module2b/coassembly_trimmed_manifest.tsv
        """.stripIndent()
        )
    }

    if (trimmed_manifest_files.isEmpty()) {
        error(
            """
        No read manifests were found for Module 3.

        This should not happen if assembly manifests were detected.
        Please check Module 1 and/or Module 2b outputs.
        """.stripIndent()
        )
    }

    log.info("Module 3 results directory: ${params.results_dir}")
    log.info("Writing Module 3 outputs to: ${params.outdir}")
    log.info("Minimum scaffold length for binning: ${params.min_scaffold_length}")
    log.info("Binners selected: MetaBAT2=${use_metabat2}, QuickBin=${use_quickbin}, MaxBin2=${use_maxbin2}")
    log.info("Regular Module 2 discovery enabled: ${include_module2}")
    log.info("Module 2b coassembly discovery enabled: ${include_module2b}")

    log.info("Assembly manifests selected for Module 3:")
    assembly_manifest_files.each { manifest_path ->
        log.info("  ${manifest_path}")
    }

    log.info("Read manifests selected for Module 3:")
    trimmed_manifest_files.each { manifest_path ->
        log.info("  ${manifest_path}")
    }

    def assembly_manifest_ch = channel.fromList(assembly_manifest_files)
        .map { manifest_path ->
            file(manifest_path)
        }

    def trimmed_manifest_ch = channel.fromList(trimmed_manifest_files)
        .map { manifest_path ->
            file(manifest_path)
        }

    /*
     * Module 2 and Module 2b assembly_manifest.tsv columns:
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
            def fasta_path = absOrEmpty(row.renamed_fasta)

            if (!fasta_path) {
                error("Empty renamed_fasta path in assembly manifest for sample '${row.sample_id}'")
            }

            tuple(
                row.sample_id.toString(),
                row.safe_sample_id.toString(),
                row.assembly_sample_id.toString(),
                row.assembler.toString(),
                row.assembly_mode.toString(),
                row.rarefaction_label == null ? "" : row.rarefaction_label.toString(),
                row.assembly_strategy.toString(),
                file(fasta_path),
            )
        }

    /*
     * Module 1 and Module 2b coassembly trimmed_manifest.tsv columns:
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
                absOrEmpty(row.interleaved),
            )
        }

    /*
     * Important:
     *
     * Do not use join() here.
     *
     * There can be multiple assemblies per sample:
     *   sample + assembler + single
     *   sample + assembler + rarefied a
     *   sample + assembler + rarefied b
     *   coassembly group + megahit + coassembly
     *
     * join() is unsafe with duplicate keys. combine(..., by: 0)
     * gives the desired one-to-many match by sample_id.
     */
    def binning_jobs_ch = assemblies_ch
        .combine(trimmed_reads_ch, by: 0)
        .map { sample_id, safe_sample_id, assembly_sample_id, assembly_assembler, assembly_mode, rarefaction_label, assembly_strategy, renamed_fasta, layout, read1, read2, interleaved ->
            def rare_part = rarefaction_label ? "_${rarefaction_label}" : ""
            def binning_id = "${safe_sample_id}${rare_part}_${assembly_assembler}_${assembly_mode}".replaceAll('[^A-Za-z0-9._-]+', '_')

            if (layout != 'paired' && layout != 'interleaved') {
                error("Unsupported layout in trimmed manifest for sample '${sample_id}': ${layout}")
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
                interleaved,
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
    def stats_files_ch = channel.empty()

    if (use_metabat2 || use_maxbin2) {
        GENERATE_COVERAGE_FILES(
            MAP_READS_BBMAP.out.mapped_bam.combine(SETUP_MODULE3_TOOLS.out.status)
        )
    }

    if (use_metabat2) {
        RUN_METABAT2(
            GENERATE_COVERAGE_FILES.out.coverage_files.combine(SETUP_MODULE3_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(RUN_METABAT2.out.manifest_record)
        stats_files_ch = stats_files_ch.mix(RUN_METABAT2.out.stats_file)
    }

    if (use_maxbin2) {
        RUN_MAXBIN2(
            GENERATE_COVERAGE_FILES.out.coverage_files.combine(SETUP_MODULE3_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(RUN_MAXBIN2.out.manifest_record)
        stats_files_ch = stats_files_ch.mix(RUN_MAXBIN2.out.stats_file)
    }

    if (use_quickbin) {
        RUN_QUICKBIN(
            MAP_READS_BBMAP.out.mapped_bam.combine(SETUP_MODULE3_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(RUN_QUICKBIN.out.manifest_record)
        stats_files_ch = stats_files_ch.mix(RUN_QUICKBIN.out.stats_file)
    }

    WRITE_BINNING_SUMMARIES(
        manifest_records_ch.collect(),
        stats_files_ch.collect(),
    )
}

process SETUP_MODULE3_TOOLS {
    tag "setup_binning_tools"

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "module3_tools_status.env"

    output:
    path "module3_tools_status.env", emit: status

    script:
    def base_env_dir = params.tool_env_dir ?: "${params.outdir}/conda_envs"

    def mapping_env_dir = "${base_env_dir}/module3_mapping_tools"
    def metabat2_env_dir = "${base_env_dir}/module3_metabat2_tools"
    def maxbin2_env_dir = "${base_env_dir}/module3_maxbin2_tools"

    def want_metabat2 = params.metabat2.toString().toBoolean()
    def want_quickbin = params.quickbin.toString().toBoolean()
    def want_maxbin2 = params.maxbin2.toString().toBoolean()
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

    echo "Module 3 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested MetaBAT2: \$WANT_METABAT2" >> "\$STATUS_FILE"
    echo "Requested QuickBin: \$WANT_QUICKBIN" >> "\$STATUS_FILE"
    echo "Requested MaxBin2: \$WANT_MAXBIN2" >> "\$STATUS_FILE"
    echo "Need JGI depth utility: \$NEED_JGI_DEPTH" >> "\$STATUS_FILE"
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

    create_env() {
        local env_path="\$1"
        shift

        if [[ "${params.auto_install}" != "true" ]]; then
            echo "ERROR: Environment missing/broken and --auto_install false: \$env_path" >> "\$STATUS_FILE"
            exit 1
        fi

        if [[ -z "\$INSTALLER" ]]; then
            echo "ERROR: Neither mamba nor conda was found in PATH." >> "\$STATUS_FILE"
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

    check_mapping_env() {
        local env_path="\$1"
        local ok="true"

        echo "Checking mapping/QuickBin environment: \$env_path" >> "\$STATUS_FILE"

        [[ -x "\$env_path/bin/seqkit" ]] || ok="false"
        [[ -x "\$env_path/bin/samtools" ]] || ok="false"

        if [[ ! -x "\$env_path/bin/java" && -x "\$env_path/lib/jvm/bin/java" ]]; then
            mkdir -p "\$env_path/bin"
            ln -sfn "\$env_path/lib/jvm/bin/java" "\$env_path/bin/java"
        fi

        [[ -x "\$env_path/bin/java" ]] || ok="false"

        local jar_count
        jar_count="\$(find "\$env_path" -type f -name 'bbtools.jar' 2>/dev/null | wc -l | tr -d ' ')"
        [[ "\$jar_count" -ge 1 ]] || ok="false"

        if [[ "\$WANT_QUICKBIN" == "true" ]]; then
            if ! find "\$env_path" -type f -name 'bbtools.jar' 2>/dev/null | grep -q .; then
                ok="false"
            fi
        fi

        [[ "\$ok" == "true" ]]
    }

    check_metabat2_env() {
        local env_path="\$1"
        local ok="true"

        echo "Checking MetaBAT2/JGI environment: \$env_path" >> "\$STATUS_FILE"

        if [[ "\$WANT_METABAT2" == "true" ]]; then
            [[ -x "\$env_path/bin/metabat2" ]] || ok="false"
        fi

        if [[ "\$NEED_JGI_DEPTH" == "true" ]]; then
            [[ -x "\$env_path/bin/jgi_summarize_bam_contig_depths" ]] || ok="false"
        fi

        [[ "\$ok" == "true" ]]
    }

    check_maxbin2_env() {
        local env_path="\$1"
        local ok="true"

        echo "Checking MaxBin2 environment: \$env_path" >> "\$STATUS_FILE"

        if [[ "\$WANT_MAXBIN2" == "true" ]]; then
            [[ -x "\$env_path/bin/run_MaxBin.pl" || -s "\$env_path/bin/run_MaxBin.pl" ]] || ok="false"
            [[ -x "\$env_path/bin/perl" ]] || ok="false"
        fi

        [[ "\$ok" == "true" ]]
    }

    echo "=== Mapping / BBMap / QuickBin environment ===" >> "\$STATUS_FILE"

    if [[ -d "\$MAPPING_ENV" ]]; then
        if check_mapping_env "\$MAPPING_ENV"; then
            echo "Existing mapping environment passed checks." >> "\$STATUS_FILE"
        else
            echo "Existing mapping environment failed checks. Removing." >> "\$STATUS_FILE"
            rm -rf "\$MAPPING_ENV"
        fi
    fi

    if [[ ! -d "\$MAPPING_ENV" ]]; then
        create_env "\$MAPPING_ENV" \\
            "python=3.*" \\
            "openjdk=17.*" \\
            "seqkit=${params.seqkit_version}" \\
            "bbmap=${params.bbmap_version}" \\
            "samtools>=${params.samtools_version},<2.0a0"

        if ! check_mapping_env "\$MAPPING_ENV"; then
            echo "ERROR: Newly created mapping environment failed checks." >> "\$STATUS_FILE"
            exit 1
        fi
    fi

    echo "=== MetaBAT2 / JGI environment ===" >> "\$STATUS_FILE"

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
            create_env "\$METABAT2_ENV" \\
                "python=3.*" \\
                "metabat2=${params.metabat2_version}"

            if ! check_metabat2_env "\$METABAT2_ENV"; then
                echo "ERROR: Newly created MetaBAT2 environment failed checks." >> "\$STATUS_FILE"
                exit 1
            fi
        fi
    else
        METABAT2_ENV="NOT_USED"
    fi

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
            create_env "\$MAXBIN2_ENV" \\
                "python=3.*" \\
                "perl" \\
                "maxbin2=${params.maxbin2_version}"

            if ! check_maxbin2_env "\$MAXBIN2_ENV"; then
                echo "ERROR: Newly created MaxBin2 environment failed checks." >> "\$STATUS_FILE"
                exit 1
            fi
        fi
    else
        MAXBIN2_ENV="NOT_USED"
    fi

    echo "----------------------------------------" >> "\$STATUS_FILE"
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

    publishDir "${params.outdir}/filtered_assemblies", mode: params.publish_filtered_assemblies_mode, pattern: "*.min${params.min_scaffold_length}.fa"

    publishDir "${params.outdir}/summary/filtered_assembly_stats", mode: 'copy', pattern: "*.filtered_assembly_stats.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.filter_assembly.log"

    input:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(original_assembly_fasta), val(layout), val(read1), val(read2), val(interleaved), path(tools_status)

    output:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path("${binning_id}.min${params.min_scaffold_length}.fa"), val(layout), val(read1), val(read2), val(interleaved), emit: filtered_jobs, optional: true

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
        echo "ERROR: seqkit is not available after setup." >&2
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
        orig_count - filt_count,
        orig_bp - filt_bp,
        published_filtered_fasta,
        sep="\\t",
        file=out
    )

if filt_count == 0:
    try:
        Path(filtered_fasta).unlink()
    except FileNotFoundError:
        pass
    print(
        f"WARNING: No contigs >= {min_len} bp remained after filtering for {binning_id}; "
        "skipping mapping and binning for this assembly.",
        file=sys.stderr
    )
PY

    if [[ ! -e "\$OUT_FASTA" ]]; then
        echo "WARNING: No contigs >= ${params.min_scaffold_length} bp remained after filtering; skipping mapping and binning for ${binning_id}." >> "\$LOG_FILE"
    fi
    
    echo "Assembly length filtering finished: \$(date)" >> "\$LOG_FILE"
    """
}


process MAP_READS_BBMAP {
    tag { "${sample_id}:${binning_id}" }

    stageInMode 'symlink'

    publishDir "${params.outdir}/mapping", mode: params.publish_bam_mode, pattern: "*.sorted.bam*"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.bbmap.log"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.mapping_threads as int
    }

    input:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(filtered_assembly_fasta), val(layout), val(read1), val(read2), val(interleaved), path(tools_status)

    output:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(filtered_assembly_fasta), val(layout), val(read1), val(read2), val(interleaved), path("${binning_id}.sorted.bam"), path("${binning_id}.sorted.bam.bai"), emit: mapped_bam

    path "${binning_id}.bbmap.log", emit: log_file

    script:
    """
    set -euo pipefail

    MAPPING_ENV="\$(grep '^MAPPING_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$MAPPING_ENV" && "\$MAPPING_ENV" != "SYSTEM" && "\$MAPPING_ENV" != "NOT_USED" ]]; then
        export PATH="\$MAPPING_ENV/bin:\$PATH"
    fi

    if ! command -v samtools >/dev/null 2>&1; then
        echo "ERROR: samtools is not available after setup." >&2
        exit 1
    fi

    if ! command -v java >/dev/null 2>&1; then
        echo "ERROR: java is not available after setup." >&2
        exit 1
    fi

    LOG_FILE="${binning_id}.bbmap.log"

    echo "BBMap mapping started: \$(date)" > "\$LOG_FILE"
    echo "Sample ID: ${sample_id}" >> "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "Layout: ${layout}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    cp -L "${filtered_assembly_fasta}" assembly.fa

    BBTOOLS_JAR="\$(find "\$MAPPING_ENV" -type f -name 'bbtools.jar' 2>/dev/null | head -n 1 || true)"

    if [[ -z "\$BBTOOLS_JAR" ]]; then
        echo "ERROR: Could not find bbtools.jar under MAPPING_ENV=\$MAPPING_ENV" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    BBMAP_COMMON_ARGS="ref=assembly.fa ambig=${params.bbmap_ambig} minid=${params.bbmap_minid} maxindel=${params.bbmap_maxindel} threads=${task.cpus} overwrite=t nodisk=t"

    if [[ "${params.bbmap_mateqtag.toString().toLowerCase()}" == "true" ]]; then
        BBMAP_COMMON_ARGS="\${BBMAP_COMMON_ARGS} mateqtag=t"
    fi

    if [[ -n "${params.bbmap_extra_args}" ]]; then
        BBMAP_COMMON_ARGS="\${BBMAP_COMMON_ARGS} ${params.bbmap_extra_args}"
    fi

    if [[ "${layout}" == "paired" ]]; then
        ln -sfn "${read1}" input_R1.fastq.gz
        ln -sfn "${read2}" input_R2.fastq.gz

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
            out=mapped.sam \\
            >> "\$LOG_FILE" 2>&1

    elif [[ "${layout}" == "interleaved" ]]; then
        ln -sfn "${interleaved}" input_interleaved.fastq.gz

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
            out=mapped.sam \\
            >> "\$LOG_FILE" 2>&1
    else
        echo "ERROR: Unsupported layout: ${layout}" >> "\$LOG_FILE"
        exit 1
    fi

    if [[ ! -s mapped.sam ]]; then
        echo "ERROR: BBMap did not produce mapped.sam" >> "\$LOG_FILE"
        cat "\$LOG_FILE" >&2
        exit 1
    fi

    samtools view -@ ${task.cpus} -bS mapped.sam \\
        | samtools sort -@ ${task.cpus} -o "${binning_id}.sorted.bam" -

    samtools index "${binning_id}.sorted.bam"

    rm -f mapped.sam assembly.fa input_R1.fastq.gz input_R2.fastq.gz input_interleaved.fastq.gz

    echo "BBMap mapping finished: \$(date)" >> "\$LOG_FILE"
    """
}


process GENERATE_COVERAGE_FILES {
    tag { binning_id }

    stageInMode 'symlink'

    publishDir "${params.outdir}/coverage/metabat2_depth", mode: 'copy', pattern: "*.depth.txt"

    publishDir "${params.outdir}/coverage/maxbin2_abundance", mode: 'copy', pattern: "*.maxbin2_abundance.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.coverage.log"

    input:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(filtered_assembly_fasta), val(layout), val(read1), val(read2), val(interleaved), path(sorted_bam), path(sorted_bam_bai), path(tools_status)

    output:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(filtered_assembly_fasta), path(sorted_bam), path(sorted_bam_bai), path("${binning_id}.depth.txt"), path("${binning_id}.maxbin2_abundance.tsv"), emit: coverage_files

    path "${binning_id}.coverage.log", emit: log_file

    script:
    """
    set -euo pipefail

    JGI_ENV="\$(grep '^JGI_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$JGI_ENV" && "\$JGI_ENV" != "SYSTEM" && "\$JGI_ENV" != "NOT_USED" ]]; then
        export PATH="\$JGI_ENV/bin:\$PATH"
    fi

    if ! command -v jgi_summarize_bam_contig_depths >/dev/null 2>&1; then
        echo "ERROR: jgi_summarize_bam_contig_depths is not available after setup." >&2
        exit 1
    fi

    LOG_FILE="${binning_id}.coverage.log"

    echo "Coverage generation started: \$(date)" > "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "BAM: ${sorted_bam}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    jgi_summarize_bam_contig_depths \\
        --outputDepth "${binning_id}.depth.txt" \\
        "${sorted_bam}" \\
        >> "\$LOG_FILE" 2>&1

    if [[ ! -s "${binning_id}.depth.txt" ]]; then
        echo "ERROR: jgi_summarize_bam_contig_depths did not produce depth file." >> "\$LOG_FILE"
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
        echo "WARNING: MaxBin2 abundance file was empty; creating empty abundance file." >> "\$LOG_FILE"
        : > "${binning_id}.maxbin2_abundance.tsv"
    fi

    echo "Coverage generation finished: \$(date)" >> "\$LOG_FILE"
    """
}


process RUN_METABAT2 {
    tag { binning_id }
    stageInMode 'symlink'

    publishDir { "${params.outdir}/bins/metabat2/${binning_id}" }, mode: params.publish_bins_mode, pattern: "bins/*", saveAs: { filename -> filename.replaceFirst(/^bins\//, '') }

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.metabat2.log"

    publishDir "${params.outdir}/summary/per_binner_stats", mode: 'copy', pattern: "*.metabat2.binning_stats.tsv"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.binning_threads as int
    }

    input:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(filtered_assembly_fasta), path(sorted_bam), path(sorted_bam_bai), path(depth_file), path(maxbin_abundance), path(tools_status)

    output:
    path "bins/*", optional: true, emit: bin_files
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

    BINNER="metabat2"
    BIN_DIR="bins"
    LOG_FILE="${binning_id}.metabat2.log"
    MANIFEST_FILE="${binning_id}.metabat2.binning_manifest_record.tsv"
    STATS_FILE="${binning_id}.metabat2.binning_stats.tsv"

    mkdir -p "\$BIN_DIR"

    echo "MetaBAT2 binning started: \$(date)" > "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "Depth file: ${depth_file}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    BINNER_EXIT_STATUS=0
    BINNER_STATUS="completed"
    BINNER_MESSAGE="MetaBAT2 completed"

    if ! command -v metabat2 >/dev/null 2>&1; then
        BINNER_EXIT_STATUS=127
        BINNER_STATUS="tool_missing"
        BINNER_MESSAGE="metabat2 was not available after setup"
        echo "WARNING: \$BINNER_MESSAGE" >> "\$LOG_FILE"
    else
        set +e
        metabat2 \\
            -i "${filtered_assembly_fasta}" \\
            -a "${depth_file}" \\
            -o "\$BIN_DIR/${binning_id}.metabat2_bin" \\
            -t ${task.cpus} \\
            -m ${params.metabat2_min_contig} \\
            ${params.metabat2_extra_args} \\
            >> "\$LOG_FILE" 2>&1
        BINNER_EXIT_STATUS="\$?"
        set -e

        if [[ "\$BINNER_EXIT_STATUS" -ne 0 ]]; then
            BINNER_STATUS="failed_nonfatal"
            BINNER_MESSAGE="MetaBAT2 exited non-zero; treating as non-fatal"
            echo "WARNING: \$BINNER_MESSAGE. Exit status: \$BINNER_EXIT_STATUS" >> "\$LOG_FILE"
        fi
    fi

    python3 - \\
        "\$BIN_DIR" \\
        "\$MANIFEST_FILE" \\
        "\$STATS_FILE" \\
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
        "" \\
        "\$BINNER_EXIT_STATUS" \\
        "\$BINNER_STATUS" \\
        "\$BINNER_MESSAGE" <<'PY'
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
    published_report,
    binner_exit_status,
    binner_status,
    binner_message
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

fasta_exts = (".fa", ".fasta", ".fna", ".fas")
bin_files = []

if bin_dir.exists():
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

if len(bin_files) == 0 and binner_status == "completed":
    binner_status = "no_bins"
    binner_message = "Binner completed but produced zero bins"

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
        "binner_exit_status",
        "binner_status",
        "binner_message",
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
        binner_exit_status,
        binner_status,
        binner_message,
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

    publishDir { "${params.outdir}/bins/maxbin2/${binning_id}" }, mode: params.publish_bins_mode, pattern: "bins/*", saveAs: { filename -> filename.replaceFirst(/^bins\//, '') }

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.maxbin2.log"

    publishDir "${params.outdir}/summary/per_binner_stats", mode: 'copy', pattern: "*.maxbin2.binning_stats.tsv"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.binning_threads as int
    }

    input:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(filtered_assembly_fasta), path(sorted_bam), path(sorted_bam_bai), path(depth_file), path(maxbin_abundance), path(tools_status)

    output:
    path "bins/*", optional: true, emit: bin_files
    path "${binning_id}.maxbin2.binning_manifest_record.tsv", emit: manifest_record
    path "${binning_id}.maxbin2.binning_stats.tsv", emit: stats_file
    path "${binning_id}.maxbin2.log", emit: log_file

    script:
    """
    set -euo pipefail

    MAXBIN2_ENV="\$(grep '^MAXBIN2_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"
    if [[ -n "\$MAXBIN2_ENV" && "\$MAXBIN2_ENV" != "SYSTEM" && "\$MAXBIN2_ENV" != "NOT_USED" ]]; then
        export PATH="\$MAXBIN2_ENV/bin:\$PATH"
    fi

    BINNER="maxbin2"
    BIN_DIR="bins"
    LOG_FILE="${binning_id}.maxbin2.log"
    MANIFEST_FILE="${binning_id}.maxbin2.binning_manifest_record.tsv"
    STATS_FILE="${binning_id}.maxbin2.binning_stats.tsv"

    mkdir -p "\$BIN_DIR"

    echo "MaxBin2 binning started: \$(date)" > "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "Abundance file: ${maxbin_abundance}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    BINNER_EXIT_STATUS=0
    BINNER_STATUS="completed"
    BINNER_MESSAGE="MaxBin2 completed"

    if ! command -v run_MaxBin.pl >/dev/null 2>&1; then
        BINNER_EXIT_STATUS=127
        BINNER_STATUS="tool_missing"
        BINNER_MESSAGE="run_MaxBin.pl was not available after setup"
        echo "WARNING: \$BINNER_MESSAGE" >> "\$LOG_FILE"
    elif [[ ! -s "${filtered_assembly_fasta}" ]]; then
        BINNER_EXIT_STATUS=2
        BINNER_STATUS="input_missing"
        BINNER_MESSAGE="Filtered assembly FASTA missing or empty"
        echo "WARNING: \$BINNER_MESSAGE" >> "\$LOG_FILE"
    elif [[ ! -s "${maxbin_abundance}" ]]; then
        BINNER_EXIT_STATUS=3
        BINNER_STATUS="input_missing"
        BINNER_MESSAGE="MaxBin2 abundance file missing or empty"
        echo "WARNING: \$BINNER_MESSAGE" >> "\$LOG_FILE"
    else
        ln -sfn "${filtered_assembly_fasta}" maxbin2_contigs.fa
        ln -sfn "${maxbin_abundance}" maxbin2_abundance.tsv

        set +e
        run_MaxBin.pl \\
            -contig maxbin2_contigs.fa \\
            -abund maxbin2_abundance.tsv \\
            -out "\$BIN_DIR/${binning_id}.maxbin2_bin" \\
            -thread ${task.cpus} \\
            ${params.maxbin2_extra_args} \\
            >> "\$LOG_FILE" 2>&1
        BINNER_EXIT_STATUS="\$?"
        set -e

        if [[ "\$BINNER_EXIT_STATUS" -ne 0 ]]; then
            BINNER_STATUS="failed_nonfatal"
            BINNER_MESSAGE="MaxBin2 exited non-zero; treating as non-fatal"
            echo "WARNING: \$BINNER_MESSAGE. Exit status: \$BINNER_EXIT_STATUS" >> "\$LOG_FILE"
        fi
    fi

    find "\$BIN_DIR" -type f -print >> "\$LOG_FILE" 2>&1 || true

    python3 - \\
        "\$BIN_DIR" \\
        "\$MANIFEST_FILE" \\
        "\$STATS_FILE" \\
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
        "" \\
        "\$BINNER_EXIT_STATUS" \\
        "\$BINNER_STATUS" \\
        "\$BINNER_MESSAGE" <<'PY'
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
    published_report,
    binner_exit_status,
    binner_status,
    binner_message
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

fasta_exts = (".fa", ".fasta", ".fna", ".fas")
bin_files = []

if bin_dir.exists():
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

if len(bin_files) == 0 and binner_status == "completed":
    binner_status = "no_bins"
    binner_message = "Binner completed but produced zero bins"

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
        "binner_exit_status",
        "binner_status",
        "binner_message",
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
        binner_exit_status,
        binner_status,
        binner_message,
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

    publishDir { "${params.outdir}/bins/quickbin/${binning_id}" }, mode: params.publish_bins_mode, pattern: "bins/*", saveAs: { filename -> filename.replaceFirst(/^bins\//, '') }

    publishDir "${params.outdir}/coverage/quickbin_cov", mode: 'copy', pattern: "*.quickbin_cov.txt"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.quickbin.*"

    publishDir "${params.outdir}/summary/per_binner_stats", mode: 'copy', pattern: "*.quickbin.binning_stats.tsv"

    cpus {
        params.threads != null
            ? params.threads as int
            : params.binning_threads as int
    }

    input:
    tuple val(sample_id), val(safe_sample_id), val(assembly_sample_id), val(assembly_assembler), val(assembly_mode), val(rarefaction_label), val(assembly_strategy), val(binning_id), path(filtered_assembly_fasta), val(layout), val(read1), val(read2), val(interleaved), path(sorted_bam), path(sorted_bam_bai), path(tools_status)

    output:
    path "bins/*", optional: true, emit: bin_files
    path "${binning_id}.quickbin_cov.txt", emit: cov_file
    path "${binning_id}.quickbin.report.tsv", emit: report_file
    path "${binning_id}.quickbin.binning_manifest_record.tsv", emit: manifest_record
    path "${binning_id}.quickbin.binning_stats.tsv", emit: stats_file
    path "${binning_id}.quickbin.log", emit: log_file

    script:
    def quickbin_xmx_arg = ""
    if (params.quickbin_xmx) {
        def x = params.quickbin_xmx.toString()
        quickbin_xmx_arg = x.startsWith("-Xmx") ? x : "-Xmx${x}"
    }

    """
    set -euo pipefail

    QUICKBIN_ENV="\$(grep '^QUICKBIN_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"
    if [[ -n "\$QUICKBIN_ENV" && "\$QUICKBIN_ENV" != "SYSTEM" && "\$QUICKBIN_ENV" != "NOT_USED" ]]; then
        export PATH="\$QUICKBIN_ENV/bin:\$PATH"
    fi

    BINNER="quickbin"
    BIN_DIR="bins"
    LOG_FILE="${binning_id}.quickbin.log"
    REPORT_FILE="${binning_id}.quickbin.report.tsv"
    COV_FILE="${binning_id}.quickbin_cov.txt"
    MANIFEST_FILE="${binning_id}.quickbin.binning_manifest_record.tsv"
    STATS_FILE="${binning_id}.quickbin.binning_stats.tsv"

    mkdir -p "\$BIN_DIR"

    echo "QuickBin binning started: \$(date)" > "\$LOG_FILE"
    echo "Binning ID: ${binning_id}" >> "\$LOG_FILE"
    echo "Filtered assembly FASTA: ${filtered_assembly_fasta}" >> "\$LOG_FILE"
    echo "BAM: ${sorted_bam}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    BINNER_EXIT_STATUS=0
    BINNER_STATUS="completed"
    BINNER_MESSAGE="QuickBin completed"

    BBTOOLS_JAR="\$(find "\$QUICKBIN_ENV" -type f -name 'bbtools.jar' 2>/dev/null | head -n 1 || true)"

    if [[ -z "\$BBTOOLS_JAR" ]]; then
        BINNER_EXIT_STATUS=127
        BINNER_STATUS="tool_missing"
        BINNER_MESSAGE="bbtools.jar was not available for QuickBin"
        echo "WARNING: \$BINNER_MESSAGE" >> "\$LOG_FILE"
    else
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

        printf 'QuickBin Java command:' >> "\$LOG_FILE"
        printf ' %q' "\${JAVA_CMD[@]}" >> "\$LOG_FILE"
        printf '\\n' >> "\$LOG_FILE"

        set +e
        "\${JAVA_CMD[@]}" >> "\$LOG_FILE" 2>&1
        BINNER_EXIT_STATUS="\$?"
        set -e

        if [[ "\$BINNER_EXIT_STATUS" -ne 0 ]]; then
            BINNER_STATUS="failed_nonfatal"
            BINNER_MESSAGE="QuickBin exited non-zero; treating as non-fatal"
            echo "WARNING: \$BINNER_MESSAGE. Exit status: \$BINNER_EXIT_STATUS" >> "\$LOG_FILE"
        fi
    fi

    mkdir -p "\$BIN_DIR"

    if [[ ! -s "\$COV_FILE" ]]; then
        echo "WARNING: QuickBin did not create a non-empty coverage file: \$COV_FILE" >> "\$LOG_FILE"
        : > "\$COV_FILE"
    fi

    if [[ ! -s "\$REPORT_FILE" ]]; then
        echo "WARNING: QuickBin did not create a non-empty report file: \$REPORT_FILE" >> "\$LOG_FILE"
        : > "\$REPORT_FILE"
    fi

    python3 - \\
        "\$BIN_DIR" \\
        "\$MANIFEST_FILE" \\
        "\$STATS_FILE" \\
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
        "${params.outdir}/logs/${binning_id}.quickbin.report.tsv" \\
        "\$BINNER_EXIT_STATUS" \\
        "\$BINNER_STATUS" \\
        "\$BINNER_MESSAGE" <<'PY'
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
    published_report,
    binner_exit_status,
    binner_status,
    binner_message
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

fasta_exts = (".fa", ".fasta", ".fna", ".fas")
bin_files = []

if bin_dir.exists():
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

if len(bin_files) == 0 and binner_status == "completed":
    binner_status = "no_bins"
    binner_message = "Binner completed but produced zero bins"

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
        "binner_exit_status",
        "binner_status",
        "binner_message",
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
        binner_exit_status,
        binner_status,
        binner_message,
        sep="\\t",
        file=stats
    )
PY

    echo "QuickBin binning finished: \$(date)" >> "\$LOG_FILE"
    """
}

process WRITE_FILTERED_ASSEMBLY_STATS_SUMMARY {
    tag "write_filtered_assembly_stats_summary"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "filtered_assembly_stats_summary.tsv"

    input:
    path stats_files

    output:
    path "filtered_assembly_stats_summary.tsv", emit: summary

    script:
    def stats_file_list = stats_files.collect { stats -> stats.name }.join(' ')

    """
    set -euo pipefail

    if [[ -z "${stats_file_list}" ]]; then
        echo "WARNING: No binning stats files were received; no assemblies passed the minimum scaffold length filter." >&2
        printf 'sample_id\tsafe_sample_id\tassembly_sample_id\tassembly_assembler\tassembly_mode\trarefaction_label\tassembly_strategy\tbinning_id\tbinner\tbin_count\ttotal_bin_bp\tlargest_bin_bp\tbins_dir\tlog_file\treport_file\tbinner_exit_status\tbinner_status\tbinner_message\n' > binning_stats_summary.tsv
        exit 0
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

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "binning_manifest.tsv"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "binning_stats_summary.tsv"

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

    printf 'sample_id\\tsafe_sample_id\\tassembly_sample_id\\tassembly_assembler\\tassembly_mode\\trarefaction_label\\tassembly_strategy\\tbinning_id\\tbinner\\tbin_id\\tbin_bp\\tbin_fasta\\n' > binning_manifest.tsv

    if [[ -n "${manifest_files}" ]]; then
        for f in ${manifest_files}; do
            if [[ -s "\$f" ]]; then
                cat "\$f" >> binning_manifest.tsv
            fi
        done
    fi

    if [[ -z "${stats_file_list}" ]]; then
        echo "ERROR: No binning stats files were received." >&2
        exit 1
    fi

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
