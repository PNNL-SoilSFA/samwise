#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
* Module 2: Read assembly from Module 1 trimmed reads.
*
* The task helper scripts live in bin/:
*
*   bin/samwise_rarefy_reads.py     rarefied read subsetting
*   bin/samwise_run_assembler.py    MEGAHIT / metaSPAdes invocation
*   bin/samwise_rename_contigs.py   contig renaming, statistics, and manifest rows
*
* Run this workflow with bin/module_2_slurm.config to submit each single or
* rarefied assembly as an independent SLURM job. The config caps concurrent
* submissions with --max_parallel_assemblies; without it, the workflow uses the
* local executor.
*
* Output contract (assembly_manifest.tsv / assembly_stats_summary.tsv columns,
* contig header format, file naming) is unchanged, so Module 3 consumes the
* results identically in either execution mode.
*/

// samwise_dir identifies the installed SAMWISE source tree; working_dir is
// independently the results/environment root. Never infer source assets from
// working_dir, because a results-only invocation must still find this script's
// bundled helpers and dependencies.
params.samwise_dir = java.nio.file.Paths
    .get((params.samwise_dir ?: projectDir).toString())
    .toAbsolutePath()
    .normalize()
    .toString()
params.working_dir = java.nio.file.Paths
    .get((params.working_dir ?: params.samwise_dir).toString())
    .toAbsolutePath()
    .normalize()
    .toString()
params.input_manifest = null
params.megahit = false
params.metaspades = false
params.single_assembly = true
params.rarefied_assembly = false
params.rarefaction_splits = 2
params.megahit_version = "1.2.9"
params.spades_version = "4.2.0"
params.auto_install = true
params.tool_env_dir = null
params.threads = null
params.assembly_threads = 4
params.memory_gb = 0
params.clean_partial_assembler_outputs = true
params.megahit_threads = null
params.megahit_preset = "meta-large"
params.publish_assemblies_mode = "symlink"
params.results_dir = params.working_dir
params.module1_outdir = "${params.results_dir}/module_1_readtrimming"
params.outdir = "${params.results_dir}/module_2_readassembly"

def rareLabelFromIndex(index) {

    def alphabet = "abcdefghijklmnopqrstuvwxyz"

    if (index < 0) {
        error("Rarefaction label index cannot be negative: ${index}")
    }

    if (index < 26) {
        return alphabet.charAt(index).toString()
    }

    def prefix_index = index.intdiv(26) - 1
    def suffix_index = index % 26

    return rareLabelFromIndex(prefix_index) + alphabet.charAt(suffix_index).toString()
}

def rareLabels(count) {
    return (0..<count).collect { idx -> rareLabelFromIndex(idx as int) }
}

def absoluteReadPath(value) {

    def text = value?.toString()?.trim()
    return text ? file(text).toAbsolutePath().normalize().toString() : ''
}

workflow {

    def use_megahit = params.megahit.toString().toBoolean()
    def use_metaspades = params.metaspades.toString().toBoolean()

    def do_single = params.single_assembly.toString().toBoolean()
    def do_rarefied = params.rarefied_assembly.toString().toBoolean()

    if (!use_megahit && !use_metaspades) {
        error(
            """
        No assembler selected.

        Please specify at least one of:
          --megahit
          --metaspades
        """.stripIndent()
        )
    }

    if (!do_single && !do_rarefied) {
        error(
            """
        No assembly mode selected.

        Please enable at least one of:
          --single_assembly true
          --rarefied_assembly true
        """.stripIndent()
        )
    }

    def rare_split_count = params.rarefaction_splits as int

    if (rare_split_count < 2) {
        error(
            """
        Invalid rarefaction split count: ${rare_split_count}

        Rarefied assembly requires at least 2 splits.
        """.stripIndent()
        )
    }

    def manifest_file = params.input_manifest ?: "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    log.info("Module 2 results directory: ${params.results_dir}")
    log.info("Using Module 1 trimmed manifest: ${manifest_file}")
    log.info("Writing Module 2 outputs to: ${params.outdir}")
    log.info("Assembler selected: MEGAHIT=${use_megahit}, metaSPAdes=${use_metaspades}")
    log.info("Assembly modes: single=${do_single}, rarefied=${do_rarefied}")
    log.info("Global threads: ${params.threads ?: params.assembly_threads}")
    log.info("MEGAHIT thread override: ${params.megahit_threads ?: 'not supplied'}")
    log.info("Global memory: ${(params.memory_gb as int) > 0 ? params.memory_gb + ' GB' : 'not supplied'}")

    def assembler_list = []

    if (use_megahit) {
        assembler_list << "megahit"
    }

    if (use_metaspades) {
        assembler_list << "metaspades"
    }

    def manifest_ch = channel.fromPath(
        manifest_file,
        type: 'file',
        checkIfExists: true,
    )

    def reads_ch = manifest_ch
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

            if (layout != 'paired' && layout != 'interleaved') {
                error("Unsupported layout in trimmed manifest for sample '${sample_id}': ${layout}")
            }

            // Module 1 manifests created before absolute paths were recorded
            // contain paths relative to the Nextflow launch directory. Resolve
            // them here because assembly tasks execute in isolated work dirs.
            tuple(
                sample_id,
                safe_id,
                assembly_sample_id,
                layout,
                absoluteReadPath(row.read1),
                absoluteReadPath(row.read2),
                absoluteReadPath(row.interleaved),
            )
        }

    def helper_scripts_ch = channel.value(
        tuple(
            file("${params.samwise_dir}/bin/samwise_rarefy_reads.py", checkIfExists: true),
            file("${params.samwise_dir}/bin/samwise_run_assembler.py", checkIfExists: true),
            file("${params.samwise_dir}/bin/samwise_rename_contigs.py", checkIfExists: true),
        )
    )

    SETUP_MODULE2_TOOLS()

    def manifest_records_ch = channel.empty()
    def stats_files_ch = channel.empty()

    if (do_single) {

        def single_jobs_ch = reads_ch.flatMap { sample_id, safe_id, assembly_sample_id, layout, read1, read2, interleaved ->

            assembler_list.collect { assembler ->
                tuple(
                    sample_id,
                    safe_id,
                    assembly_sample_id,
                    layout,
                    read1,
                    read2,
                    interleaved,
                    assembler,
                )
            }
        }

        ASSEMBLE_SINGLE(
            single_jobs_ch
                .combine(SETUP_MODULE2_TOOLS.out.status)
                .combine(helper_scripts_ch)
        )

        manifest_records_ch = manifest_records_ch.mix(ASSEMBLE_SINGLE.out.manifest_record)
        stats_files_ch = stats_files_ch.mix(ASSEMBLE_SINGLE.out.stats_file)
    }

    if (do_rarefied) {

        def rare_letters = rareLabels(rare_split_count)

        log.info("Rarefied assembly enabled with ${rare_split_count} splits: ${rare_letters.join(', ')}")

        def rare_jobs_ch = reads_ch.flatMap { sample_id, safe_id, assembly_sample_id, layout, read1, read2, interleaved ->

            def jobs = []

            assembler_list.each { assembler ->
                (1..rare_split_count).each { idx ->
                    def rare_letter = rare_letters[idx - 1]

                    jobs << tuple(
                        sample_id,
                        safe_id,
                        assembly_sample_id,
                        layout,
                        read1,
                        read2,
                        interleaved,
                        assembler,
                        idx,
                        rare_letter,
                        rare_split_count,
                    )
                }
            }

            return jobs
        }

        ASSEMBLE_RAREFIED(
            rare_jobs_ch
                .combine(SETUP_MODULE2_TOOLS.out.status)
                .combine(helper_scripts_ch)
        )

        manifest_records_ch = manifest_records_ch.mix(ASSEMBLE_RAREFIED.out.manifest_record)
        stats_files_ch = stats_files_ch.mix(ASSEMBLE_RAREFIED.out.stats_file)
    }

    WRITE_ASSEMBLY_SUMMARIES(
        manifest_records_ch.collect(),
        stats_files_ch.collect(),
    )
}


process SETUP_MODULE2_TOOLS {

    tag "setup_assemblers"
    cache false

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "module2_tools_status.env"

    output:
    path "module2_tools_status.env", emit: status

    script:
    def env_dir = params.tool_env_dir
        ? file(params.tool_env_dir).toAbsolutePath().toString()
        : file("${params.outdir}/conda_envs/module2_tools").toAbsolutePath().toString()

    def want_megahit = params.megahit.toString().toBoolean()
    def want_metaspades = params.metaspades.toString().toBoolean()

    def packages = []
    packages << "python"

    if (want_megahit) {
        packages << "megahit=${params.megahit_version}"
    }

    if (want_metaspades) {
        def spades_version_for_conda = params.spades_version.toString().replaceFirst(/-\d+$/, '')
        packages << "spades=${spades_version_for_conda}"
    }

    def package_string = packages.collect { pkg -> "\"${pkg}\"" }.join(" \\\n        ")

    """
    set -euo pipefail

    STATUS_FILE="module2_tools_status.env"
    TOOL_ENV="${env_dir}"
    WANT_MEGAHIT="${want_megahit}"
    WANT_METASPADES="${want_metaspades}"

    echo "Module 2 tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested MEGAHIT: \$WANT_MEGAHIT" >> "\$STATUS_FILE"
    echo "Requested metaSPAdes: \$WANT_METASPADES" >> "\$STATUS_FILE"
    echo "Requested MEGAHIT version: ${params.megahit_version}" >> "\$STATUS_FILE"
    echo "Requested SPAdes/metaSPAdes version: ${params.spades_version}" >> "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    check_tool_set() {
        local prefix="\$1"
        local ok="true"

        if [[ "\$WANT_MEGAHIT" == "true" ]]; then
            if [[ "\$prefix" == "SYSTEM" ]]; then
                if ! command -v megahit >/dev/null 2>&1; then
                    ok="false"
                else
                    megahit --version >> "\$STATUS_FILE" 2>&1 || ok="false"
                fi
            else
                if [[ ! -x "\$prefix/bin/megahit" ]]; then
                    ok="false"
                else
                    "\$prefix/bin/megahit" --version >> "\$STATUS_FILE" 2>&1 || ok="false"
                fi
            fi
        fi

        if [[ "\$WANT_METASPADES" == "true" ]]; then
            if [[ "\$prefix" == "SYSTEM" ]]; then
                if ! command -v metaspades.py >/dev/null 2>&1; then
                    ok="false"
                else
                    metaspades.py --version >> "\$STATUS_FILE" 2>&1 || ok="false"
                fi
            else
                if [[ ! -x "\$prefix/bin/metaspades.py" ]]; then
                    ok="false"
                else
                    "\$prefix/bin/metaspades.py" --version >> "\$STATUS_FILE" 2>&1 || ok="false"
                fi
            fi
        fi

        [[ "\$ok" == "true" ]]
    }

    if [[ -d "\$TOOL_ENV" ]]; then
        echo "Existing Module 2 environment detected: \$TOOL_ENV" >> "\$STATUS_FILE"

        if check_tool_set "\$TOOL_ENV"; then
            echo "Existing Module 2 environment passed requested tool checks." >> "\$STATUS_FILE"
            echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
            echo "Module 2 tool setup finished: \$(date)" >> "\$STATUS_FILE"
            exit 0
        else
            echo "Existing Module 2 environment is incomplete or broken. Removing." >> "\$STATUS_FILE"
            rm -rf "\$TOOL_ENV"
        fi
    fi

    echo "Checking system/runtime tools..." >> "\$STATUS_FILE"

    if check_tool_set "SYSTEM"; then
        echo "All requested assemblers are available from system/runtime PATH." >> "\$STATUS_FILE"
        echo "TOOL_ENV=SYSTEM" >> "\$STATUS_FILE"
        echo "Module 2 tool setup finished: \$(date)" >> "\$STATUS_FILE"
        exit 0
    fi

    echo "One or more requested assemblers are unavailable from system/runtime PATH." >> "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: auto_install is false and requested assembler tools are missing." >> "\$STATUS_FILE"
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

    echo "Creating Module 2 environment: \$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Packages: ${packages.join(' ')}" >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \\
        -p "\$TOOL_ENV" \\
        -c conda-forge \\
        -c bioconda \\
        ${package_string} \\
        >> "\$STATUS_FILE" 2>&1

    if ! check_tool_set "\$TOOL_ENV"; then
        echo "ERROR: Newly created Module 2 environment failed requested tool checks." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Module 2 tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}


process ASSEMBLE_SINGLE {

    tag { "${sample_id}:${assembler}:single" }

    publishDir "${params.outdir}/assemblies", mode: params.publish_assemblies_mode, pattern: "*.renamed.fa"

    publishDir "${params.outdir}/header_maps", mode: 'copy', pattern: "*.header_map.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.log"

    publishDir "${params.outdir}/summary/per_assembly_stats", mode: 'copy', pattern: "*.assembly_stats.tsv"

    cpus {
        if (assembler == 'megahit' && params.megahit_threads != null) {
            return params.megahit_threads as int
        }

        return params.threads != null
            ? params.threads as int
            : params.assembly_threads as int
    }

    memory {
        def gb = params.memory_gb as int
        return gb > 0 ? "${gb} GB" : null
    }

    input:
    tuple val(sample_id), val(safe_id), val(assembly_sample_id), val(layout), val(read1), val(read2), val(interleaved), val(assembler), path(tools_status), path(rarefy_script), path(assembler_script), path(rename_script)

    output:
    path "*.renamed.fa", emit: renamed_contigs
    path "*.header_map.tsv", emit: header_map
    path "*.assembly_stats.tsv", emit: stats_file
    path "*.assembly_manifest_record.tsv", emit: manifest_record
    path "*.log", emit: log_file

    script:
    def clean_stale = params.clean_partial_assembler_outputs.toString().toBoolean()
    def assembly_strategy = assembler == 'megahit' ? 'A' : 'B'
    def raw_out = assembler == 'megahit' ? 'megahit_out' : 'metaspades_out'
    def published_outdir = file(params.outdir).toAbsolutePath().normalize().toString()

    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    SAMPLE_ID="${sample_id}"
    SAFE_ID="${safe_id}"
    ASSEMBLY_SAMPLE_ID="${assembly_sample_id}"
    LAYOUT="${layout}"
    ASSEMBLER="${assembler}"

    OUT_FASTA="\${SAFE_ID}_\${ASSEMBLER}_single.renamed.fa"
    HEADER_MAP="\${SAFE_ID}_\${ASSEMBLER}_single.header_map.tsv"
    STATS_FILE="\${SAFE_ID}_\${ASSEMBLER}_single.assembly_stats.tsv"
    MANIFEST_RECORD="\${SAFE_ID}_\${ASSEMBLER}_single.assembly_manifest_record.tsv"
    LOG_FILE="\${SAFE_ID}_\${ASSEMBLER}_single.log"

    if [[ "${clean_stale}" == "true" ]]; then
        rm -f "\$OUT_FASTA" "\$HEADER_MAP" "\$STATS_FILE" "\$MANIFEST_RECORD"
        rm -f input_R1.fastq.gz input_R2.fastq.gz input_interleaved.fastq.gz
    fi

    if [[ "\$LAYOUT" == "paired" ]]; then
        if [[ ! -s "${read1}" || ! -s "${read2}" ]]; then
            echo "ERROR: Missing paired reads for sample \${SAMPLE_ID}" >&2
            exit 1
        fi

        ln -sfn "${read1}" input_R1.fastq.gz
        ln -sfn "${read2}" input_R2.fastq.gz

        READ1_LOCAL="input_R1.fastq.gz"
        READ2_LOCAL="input_R2.fastq.gz"
        INTERLEAVED_LOCAL=""

    elif [[ "\$LAYOUT" == "interleaved" ]]; then
        if [[ ! -s "${interleaved}" ]]; then
            echo "ERROR: Missing interleaved reads for sample \${SAMPLE_ID}" >&2
            exit 1
        fi

        ln -sfn "${interleaved}" input_interleaved.fastq.gz

        READ1_LOCAL=""
        READ2_LOCAL=""
        INTERLEAVED_LOCAL="input_interleaved.fastq.gz"

    else
        echo "ERROR: Unsupported layout: \$LAYOUT" >&2
        exit 1
    fi

    echo "Running \${ASSEMBLER} single assembly for \${SAMPLE_ID}" >> "\$LOG_FILE"

    python3 "${assembler_script}" \\
        --assembler "\$ASSEMBLER" \\
        --layout "\$LAYOUT" \\
        --read1 "\$READ1_LOCAL" \\
        --read2 "\$READ2_LOCAL" \\
        --interleaved "\$INTERLEAVED_LOCAL" \\
        --threads ${task.cpus} \\
        --memory-gb ${params.memory_gb} \\
        --out-dir "${raw_out}" \\
        --log-file "\$LOG_FILE" \\
        --result-file assembler_result.env \\
        --megahit-preset ${params.megahit_preset} \\
        --clean-stale-output ${clean_stale}

    # Defines SRC_FASTA, ASSEMBLY_STATUS, ASSEMBLY_WARNING (shell-quoted by the script).
    source assembler_result.env

    python3 "${rename_script}" \\
        --src-fasta "\$SRC_FASTA" \\
        --out-fasta "\$OUT_FASTA" \\
        --header-map "\$HEADER_MAP" \\
        --stats-file "\$STATS_FILE" \\
        --manifest-record "\$MANIFEST_RECORD" \\
        --sample-id "\$SAMPLE_ID" \\
        --safe-id "\$SAFE_ID" \\
        --assembly-sample-id "\$ASSEMBLY_SAMPLE_ID" \\
        --assembler "\$ASSEMBLER" \\
        --mode single \\
        --rarefaction-label "" \\
        --assembly-strategy ${assembly_strategy} \\
        --assembly-status "\$ASSEMBLY_STATUS" \\
        --assembly-warning "\$ASSEMBLY_WARNING" \\
        --published-fasta "${published_outdir}/assemblies/\$OUT_FASTA"

    rm -f input_R1.fastq.gz input_R2.fastq.gz input_interleaved.fastq.gz
    rm -rf "${raw_out}"
    """
}


process ASSEMBLE_RAREFIED {

    tag { "${sample_id}:${assembler}:rarefied:${rare_letter}" }

    publishDir "${params.outdir}/assemblies", mode: params.publish_assemblies_mode, pattern: "*.renamed.fa"

    publishDir "${params.outdir}/header_maps", mode: 'copy', pattern: "*.header_map.tsv"

    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.log"

    publishDir "${params.outdir}/summary/per_assembly_stats", mode: 'copy', pattern: "*.assembly_stats.tsv"

    cpus {
        if (assembler == 'megahit' && params.megahit_threads != null) {
            return params.megahit_threads as int
        }

        return params.threads != null
            ? params.threads as int
            : params.assembly_threads as int
    }

    memory {
        def gb = params.memory_gb as int
        return gb > 0 ? "${gb} GB" : null
    }

    input:
    tuple val(sample_id), val(safe_id), val(assembly_sample_id), val(layout), val(read1), val(read2), val(interleaved), val(assembler), val(rare_index), val(rare_letter), val(rare_split_count), path(tools_status), path(rarefy_script), path(assembler_script), path(rename_script)

    output:
    path "*.renamed.fa", emit: renamed_contigs
    path "*.header_map.tsv", emit: header_map
    path "*.assembly_stats.tsv", emit: stats_file
    path "*.assembly_manifest_record.tsv", emit: manifest_record
    path "*.log", emit: log_file

    script:
    def clean_stale = params.clean_partial_assembler_outputs.toString().toBoolean()
    def assembly_strategy = assembler == 'megahit' ? 'C' : 'D'
    def raw_out = assembler == 'megahit' ? 'megahit_rarefied_out' : 'metaspades_rarefied_out'
    def rare_zero_index = (rare_index as int) - 1
    def published_outdir = file(params.outdir).toAbsolutePath().normalize().toString()

    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    SAMPLE_ID="${sample_id}"
    SAFE_ID="${safe_id}"
    ASSEMBLY_SAMPLE_ID="${assembly_sample_id}${rare_letter}"
    LAYOUT="${layout}"
    ASSEMBLER="${assembler}"

    RAREFACTION_LABEL="${rare_letter}"
    RARE_SPLIT_COUNT="${rare_split_count}"

    OUT_FASTA="\${SAFE_ID}\${RAREFACTION_LABEL}_\${ASSEMBLER}_rarefied.renamed.fa"
    HEADER_MAP="\${SAFE_ID}\${RAREFACTION_LABEL}_\${ASSEMBLER}_rarefied.header_map.tsv"
    STATS_FILE="\${SAFE_ID}\${RAREFACTION_LABEL}_\${ASSEMBLER}_rarefied.assembly_stats.tsv"
    MANIFEST_RECORD="\${SAFE_ID}\${RAREFACTION_LABEL}_\${ASSEMBLER}_rarefied.assembly_manifest_record.tsv"
    LOG_FILE="\${SAFE_ID}\${RAREFACTION_LABEL}_\${ASSEMBLER}_rarefied.log"

    SUB_R1="subset_\${RAREFACTION_LABEL}_R1.fastq.gz"
    SUB_R2="subset_\${RAREFACTION_LABEL}_R2.fastq.gz"
    SUB_12="subset_\${RAREFACTION_LABEL}_interleaved.fastq.gz"

    if [[ "${clean_stale}" == "true" ]]; then
        rm -f "\$OUT_FASTA" "\$HEADER_MAP" "\$STATS_FILE" "\$MANIFEST_RECORD"
        rm -f "\$SUB_R1" "\$SUB_R2" "\$SUB_12"
    fi

    echo "Creating rarefied subset \${RAREFACTION_LABEL} of \${RARE_SPLIT_COUNT} for \${SAMPLE_ID}" >> "\$LOG_FILE"

    python3 "${rarefy_script}" \\
        --layout "\$LAYOUT" \\
        --read1 "${read1}" \\
        --read2 "${read2}" \\
        --interleaved "${interleaved}" \\
        --split-index ${rare_zero_index} \\
        --split-count "\$RARE_SPLIT_COUNT" \\
        --out-r1 "\$SUB_R1" \\
        --out-r2 "\$SUB_R2" \\
        --out-interleaved "\$SUB_12" \\
        >> "\$LOG_FILE" 2>&1

    echo "Rarefied subset created successfully." >> "\$LOG_FILE"

    echo "Running \${ASSEMBLER} rarefied assembly for \${SAMPLE_ID}, subset \${RAREFACTION_LABEL}" >> "\$LOG_FILE"

    python3 "${assembler_script}" \\
        --assembler "\$ASSEMBLER" \\
        --layout "\$LAYOUT" \\
        --read1 "\$SUB_R1" \\
        --read2 "\$SUB_R2" \\
        --interleaved "\$SUB_12" \\
        --threads ${task.cpus} \\
        --memory-gb ${params.memory_gb} \\
        --out-dir "${raw_out}" \\
        --log-file "\$LOG_FILE" \\
        --result-file assembler_result.env \\
        --megahit-preset ${params.megahit_preset} \\
        --clean-stale-output ${clean_stale}

    # Defines SRC_FASTA, ASSEMBLY_STATUS, ASSEMBLY_WARNING (shell-quoted by the script).
    source assembler_result.env

    python3 "${rename_script}" \\
        --src-fasta "\$SRC_FASTA" \\
        --out-fasta "\$OUT_FASTA" \\
        --header-map "\$HEADER_MAP" \\
        --stats-file "\$STATS_FILE" \\
        --manifest-record "\$MANIFEST_RECORD" \\
        --sample-id "\$SAMPLE_ID" \\
        --safe-id "\$SAFE_ID" \\
        --assembly-sample-id "\$ASSEMBLY_SAMPLE_ID" \\
        --assembler "\$ASSEMBLER" \\
        --mode rarefied \\
        --rarefaction-label "\$RAREFACTION_LABEL" \\
        --assembly-strategy ${assembly_strategy} \\
        --assembly-status "\$ASSEMBLY_STATUS" \\
        --assembly-warning "\$ASSEMBLY_WARNING" \\
        --published-fasta "${published_outdir}/assemblies/\$OUT_FASTA"

    rm -f "\$SUB_R1" "\$SUB_R2" "\$SUB_12"
    rm -rf "${raw_out}"
    """
}


process WRITE_ASSEMBLY_SUMMARIES {

    tag "write_assembly_summaries"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "assembly_manifest.tsv"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "assembly_stats_summary.tsv"

    input:
    path manifest_records
    path stats_files

    output:
    path "assembly_manifest.tsv", emit: manifest
    path "assembly_stats_summary.tsv", emit: stats_summary

    script:
    def manifest_files = manifest_records.collect { record -> record.name }.join(' ')
    def stats_file_list = stats_files.collect { stats -> stats.name }.join(' ')

    """
    set -euo pipefail

    if [[ -z "${manifest_files}" ]]; then
        echo "ERROR: No assembly manifest records were received." >&2
        exit 1
    fi

    if [[ -z "${stats_file_list}" ]]; then
        echo "ERROR: No assembly stats files were received." >&2
        exit 1
    fi

    printf 'sample_id\\tsafe_sample_id\\tassembly_sample_id\\tassembler\\tassembly_mode\\trarefaction_label\\tassembly_strategy\\trenamed_fasta\\n' > assembly_manifest.tsv

    for f in ${manifest_files}; do
        cat "\$f" >> assembly_manifest.tsv
    done

    first=1
    : > assembly_stats_summary.tsv

    for f in ${stats_file_list}; do
        if [[ "\$first" -eq 1 ]]; then
            cat "\$f" >> assembly_stats_summary.tsv
            first=0
        else
            tail -n +2 "\$f" >> assembly_stats_summary.tsv
        fi
    done
    """
}
