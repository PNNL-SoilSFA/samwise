#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Module 2: Read assembly from Module 1 trimmed reads.
 *
 * Global memory flag:
 *
 *   --memory_gb 512
 *
 * This controls:
 *   1. Nextflow memory request for assembly processes
 *   2. MEGAHIT -m, converted from GB to bytes
 *   3. metaSPAdes -m, in GB
 */

/*
 * Parameters
 */
params.working_dir       = null
params.input_manifest    = null
params.output_dir        = null

params.megahit           = false
params.metaspades        = false

params.single_assembly   = true
params.rarefied_assembly = false
params.rarefaction_splits = 2

params.megahit_version   = "1.2.9"
params.spades_version    = "4.2.0"

params.auto_install      = true
params.tool_env_dir      = null

params.threads           = null
params.assembly_threads  = 4

/*
 * Global memory in GB.
 *
 * Example:
 *   --memory_gb 512
 *
 * Use 0 to leave memory unset.
 */
params.memory_gb         = 0

/*
 * Optional MEGAHIT-specific thread override.
 *
 * Example:
 *   --megahit_threads 1
 */
params.megahit_threads   = null
params.megahit_preset    = "meta-large"

params.publish_assemblies_mode = "symlink"

params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.module1_outdir = "${params.results_dir}/module_1_readtrimming"
params.outdir = "${params.results_dir}/module_2_readassembly"


def rareLabelFromIndex(int index) {
    def alphabet = "abcdefghijklmnopqrstuvwxyz"
    if( index < 0 ) {
        error "Rarefaction label index cannot be negative: ${index}"
    }
    if( index < 26 ) {
        return alphabet.charAt(index).toString()
    }
    def prefix_index = ((int)(index / 26)) - 1
    def suffix_index = index % 26
    return rareLabelFromIndex(prefix_index) + alphabet.charAt(suffix_index).toString()
}


def rareLabels(int count) {
    return (0..<count).collect { idx -> rareLabelFromIndex(idx as int) }
}


workflow {

    def use_megahit    = params.megahit.toString().toBoolean()
    def use_metaspades = params.metaspades.toString().toBoolean()

    def do_single   = params.single_assembly.toString().toBoolean()
    def do_rarefied = params.rarefied_assembly.toString().toBoolean()

    if( !use_megahit && !use_metaspades ) {
        error """
        No assembler selected.
        Please specify at least one of:
          --megahit
          --metaspades
        """.stripIndent()
    }

    if( !do_single && !do_rarefied ) {
        error """
        No assembly mode selected.
        Please enable at least one of:
          --single_assembly true
          --rarefied_assembly true
        """.stripIndent()
    }

    def rare_split_count = params.rarefaction_splits as int
    if( rare_split_count < 2 ) {
        error """
        Invalid rarefaction split count: ${rare_split_count}
        Rarefied assembly requires at least 2 splits.
        """.stripIndent()
    }

    def manifest_file = params.input_manifest ?: "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    log.info "Module 2 results directory: ${params.results_dir}"
    log.info "Using Module 1 trimmed manifest: ${manifest_file}"
    log.info "Writing Module 2 outputs to: ${params.outdir}"
    log.info "Assembler selected: MEGAHIT=${use_megahit}, metaSPAdes=${use_metaspades}"
    log.info "Assembly modes: single=${do_single}, rarefied=${do_rarefied}"
    log.info "Global threads: ${params.threads ?: params.assembly_threads}"
    log.info "MEGAHIT thread override: ${params.megahit_threads ?: 'not supplied'}"
    log.info "Global memory: ${(params.memory_gb as int) > 0 ? params.memory_gb + ' GB' : 'not supplied'}"

    def assembler_list = []
    if( use_megahit ) {
        assembler_list << "megahit"
    }
    if( use_metaspades ) {
        assembler_list << "metaspades"
    }

    def manifest_ch = channel.fromPath(
        manifest_file,
        type: 'file',
        checkIfExists: true
    )

    def reads_ch = manifest_ch
        .splitCsv(header: true, sep: '\t')
        .map { row ->

            def sample_id = row.sample_id.toString()
            def safe_id   = row.safe_sample_id.toString()
            def layout    = row.layout.toString()

            def assembly_sample_id = sample_id.replaceAll('[^A-Za-z0-9]+', '')
            if( !assembly_sample_id ) {
                assembly_sample_id = safe_id.replaceAll('[^A-Za-z0-9]+', '')
            }
            if( !assembly_sample_id ) {
                error "Could not derive non-empty assembly SampleID from sample '${sample_id}'"
            }

            if( layout != 'paired' && layout != 'interleaved' ) {
                error "Unsupported layout in trimmed manifest for sample '${sample_id}': ${layout}"
            }

            tuple(
                sample_id,
                safe_id,
                assembly_sample_id,
                layout,
                row.read1.toString(),
                row.read2.toString(),
                row.interleaved.toString()
            )
        }

    SETUP_MODULE2_TOOLS()

    def manifest_records_ch = channel.empty()
    def stats_files_ch      = channel.empty()

    if( do_single ) {

        def single_jobs_ch = reads_ch.flatMap {
            sample_id,
            safe_id,
            assembly_sample_id,
            layout,
            read1,
            read2,
            interleaved ->

            assembler_list.collect { assembler ->
                tuple(
                    sample_id,
                    safe_id,
                    assembly_sample_id,
                    layout,
                    read1,
                    read2,
                    interleaved,
                    assembler
                )
            }
        }

        ASSEMBLE_SINGLE(
            single_jobs_ch.combine(SETUP_MODULE2_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(ASSEMBLE_SINGLE.out.manifest_record)
        stats_files_ch      = stats_files_ch.mix(ASSEMBLE_SINGLE.out.stats_file)
    }

    if( do_rarefied ) {

        def rare_letters = rareLabels(rare_split_count)
        log.info "Rarefied assembly enabled with ${rare_split_count} splits: ${rare_letters.join(', ')}"

        def rare_jobs_ch = reads_ch.flatMap {
            sample_id,
            safe_id,
            assembly_sample_id,
            layout,
            read1,
            read2,
            interleaved ->

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
                        rare_split_count
                    )
                }
            }

            return jobs
        }

        ASSEMBLE_RAREFIED(
            rare_jobs_ch.combine(SETUP_MODULE2_TOOLS.out.status)
        )

        manifest_records_ch = manifest_records_ch.mix(ASSEMBLE_RAREFIED.out.manifest_record)
        stats_files_ch      = stats_files_ch.mix(ASSEMBLE_RAREFIED.out.stats_file)
    }

    WRITE_ASSEMBLY_SUMMARIES(
        manifest_records_ch.collect(),
        stats_files_ch.collect()
    )
}


process SETUP_MODULE2_TOOLS {

    tag "setup_assemblers"

    publishDir "${params.outdir}/setup",
        mode: 'copy',
        pattern: "module2_tools_status.env"

    output:
    path "module2_tools_status.env", emit: status

    script:

    def env_dir = params.tool_env_dir ?: "${params.outdir}/conda_envs/module2_tools"

    def want_megahit    = params.megahit.toString().toBoolean()
    def want_metaspades = params.metaspades.toString().toBoolean()

    def packages = []
    packages << "python"

    if( want_megahit ) {
        packages << "megahit=${params.megahit_version}"
    }

    if( want_metaspades ) {
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

    publishDir "${params.outdir}/assemblies",
        mode: params.publish_assemblies_mode,
        pattern: "*.renamed.fa"

    publishDir "${params.outdir}/header_maps",
        mode: 'copy',
        pattern: "*.header_map.tsv"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.log"

    publishDir "${params.outdir}/summary/per_assembly_stats",
        mode: 'copy',
        pattern: "*.assembly_stats.tsv"

    cpus {
        if( assembler == 'megahit' && params.megahit_threads != null ) {
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
    tuple val(sample_id),
          val(safe_id),
          val(assembly_sample_id),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          val(assembler),
          path(tools_status)

    output:
    path "*.renamed.fa", emit: renamed_contigs
    path "*.header_map.tsv", emit: header_map
    path "*.assembly_stats.tsv", emit: stats_file
    path "*.assembly_manifest_record.tsv", emit: manifest_record
    path "*.log", emit: log_file

    script:
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
    MODE="single"
    RAREFACTION_LABEL=""

    OUT_FASTA="\${SAFE_ID}_\${ASSEMBLER}_single.renamed.fa"
    HEADER_MAP="\${SAFE_ID}_\${ASSEMBLER}_single.header_map.tsv"
    STATS_FILE="\${SAFE_ID}_\${ASSEMBLER}_single.assembly_stats.tsv"
    MANIFEST_RECORD="\${SAFE_ID}_\${ASSEMBLER}_single.assembly_manifest_record.tsv"
    LOG_FILE="\${SAFE_ID}_\${ASSEMBLER}_single.log"

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

    if [[ "\$ASSEMBLER" == "megahit" ]]; then

        if ! command -v megahit >/dev/null 2>&1; then
            echo "ERROR: megahit is not available after tool setup." >&2
            cat "${tools_status}" >&2 || true
            exit 1
        fi

        ASSEMBLY_STRATEGY="A"
        RAW_OUT="megahit_out"

        MEGAHIT_MEM_ARG=""
        if [[ "${params.memory_gb}" != "0" ]]; then
            MEGAHIT_MEM_BYTES=\$(( ${params.memory_gb} * 1024 * 1024 * 1024 ))
            MEGAHIT_MEM_ARG="-m \$MEGAHIT_MEM_BYTES"
        fi

        echo "Running MEGAHIT single assembly for \${SAMPLE_ID}" > "\$LOG_FILE"
        echo "MEGAHIT threads: ${task.cpus}" >> "\$LOG_FILE"
        echo "Global memory GB: ${params.memory_gb}" >> "\$LOG_FILE"
        echo "MEGAHIT memory arg: \$MEGAHIT_MEM_ARG" >> "\$LOG_FILE"

        if [[ "\$LAYOUT" == "paired" ]]; then
            megahit \\
                -1 "\$READ1_LOCAL" \\
                -2 "\$READ2_LOCAL" \\
                -t ${task.cpus} \\
                \$MEGAHIT_MEM_ARG \\
                -o "\$RAW_OUT" \\
                --presets ${params.megahit_preset} \\
                >> "\$LOG_FILE" 2>&1
        else
            megahit \\
                --12 "\$INTERLEAVED_LOCAL" \\
                -t ${task.cpus} \\
                \$MEGAHIT_MEM_ARG \\
                -o "\$RAW_OUT" \\
                --presets ${params.megahit_preset} \\
                >> "\$LOG_FILE" 2>&1
        fi

        SRC_FASTA="\$RAW_OUT/final.contigs.fa"

        if [[ ! -s "\$SRC_FASTA" ]]; then
            echo "ERROR: MEGAHIT did not produce final.contigs.fa" >> "\$LOG_FILE"
            exit 1
        fi

    elif [[ "\$ASSEMBLER" == "metaspades" ]]; then

        if ! command -v metaspades.py >/dev/null 2>&1; then
            echo "ERROR: metaspades.py is not available after tool setup." >&2
            cat "${tools_status}" >&2 || true
            exit 1
        fi

        ASSEMBLY_STRATEGY="B"
        RAW_OUT="metaspades_out"

        SPADES_MEM_ARG=""
        if [[ "${params.memory_gb}" != "0" ]]; then
            SPADES_MEM_ARG="-m ${params.memory_gb}"
        fi

        echo "Running metaSPAdes single assembly for \${SAMPLE_ID}" > "\$LOG_FILE"
        echo "metaSPAdes threads: ${task.cpus}" >> "\$LOG_FILE"
        echo "Global memory GB: ${params.memory_gb}" >> "\$LOG_FILE"
        echo "metaSPAdes memory arg: \$SPADES_MEM_ARG" >> "\$LOG_FILE"

        if [[ "\$LAYOUT" == "paired" ]]; then
            metaspades.py \\
                -1 "\$READ1_LOCAL" \\
                -2 "\$READ2_LOCAL" \\
                -t ${task.cpus} \\
                \$SPADES_MEM_ARG \\
                -o "\$RAW_OUT" \\
                >> "\$LOG_FILE" 2>&1
        else
            metaspades.py \\
                --12 "\$INTERLEAVED_LOCAL" \\
                -t ${task.cpus} \\
                \$SPADES_MEM_ARG \\
                -o "\$RAW_OUT" \\
                >> "\$LOG_FILE" 2>&1
        fi

        if [[ -s "\$RAW_OUT/scaffolds.fasta" ]]; then
            SRC_FASTA="\$RAW_OUT/scaffolds.fasta"
        elif [[ -s "\$RAW_OUT/contigs.fasta" ]]; then
            SRC_FASTA="\$RAW_OUT/contigs.fasta"
        else
            echo "ERROR: metaSPAdes did not produce scaffolds.fasta or contigs.fasta" >> "\$LOG_FILE"
            exit 1
        fi

    else
        echo "ERROR: Unsupported assembler: \$ASSEMBLER" >&2
        exit 1
    fi

    python3 - \\
        "\$SRC_FASTA" \\
        "\$OUT_FASTA" \\
        "\$HEADER_MAP" \\
        "\$STATS_FILE" \\
        "\$MANIFEST_RECORD" \\
        "\$SAMPLE_ID" \\
        "\$SAFE_ID" \\
        "\$ASSEMBLY_SAMPLE_ID" \\
        "\$ASSEMBLER" \\
        "\$MODE" \\
        "\$RAREFACTION_LABEL" \\
        "\$ASSEMBLY_STRATEGY" \\
        "${params.outdir}/assemblies/\$OUT_FASTA" <<'PY'
import re
import sys
from pathlib import Path

(
    src_fasta,
    out_fasta,
    header_map,
    stats_file,
    manifest_record,
    sample_id,
    safe_id,
    assembly_sample_id,
    assembler,
    mode,
    rarefaction_label,
    assembly_strategy,
    published_fasta
) = sys.argv[1:]

src_fasta = Path(src_fasta)
out_fasta = Path(out_fasta)
header_map = Path(header_map)
stats_file = Path(stats_file)
manifest_record = Path(manifest_record)

lengths = []
contig_count = 0
current_len = 0

def n50(vals):
    if not vals:
        return 0
    vals = sorted(vals, reverse=True)
    half = sum(vals) / 2
    running = 0
    for v in vals:
        running += v
        if running >= half:
            return v
    return 0

with src_fasta.open() as inp, out_fasta.open("w") as out, header_map.open("w") as hmap:
    print("old_header", "new_header", sep="\\t", file=hmap)

    for line in inp:
        line = line.rstrip("\\n")

        if line.startswith(">"):
            if contig_count > 0:
                lengths.append(current_len)

            current_len = 0
            contig_count += 1

            old_header = line[1:].strip()
            first_token = old_header.split()[0] if old_header else f"contig_{contig_count}"

            if assembler == "megahit":
                m = re.search(r'(k\\d+)_(\\d+)', first_token)
                if m:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_{m.group(1)}_{m.group(2)}"
                else:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_k000_{contig_count}"

            elif assembler == "metaspades":
                m = re.search(r'NODE_(\\d+)', first_token)
                if m:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_NODE_{m.group(1)}"
                else:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_NODE_{contig_count}"

            else:
                raise RuntimeError(f"Unsupported assembler: {assembler}")

            print(f">{new_header}", file=out)
            print(old_header, new_header, sep="\\t", file=hmap)

        else:
            seq = line.strip()
            current_len += len(seq)
            print(seq, file=out)

if contig_count > 0:
    lengths.append(current_len)

total_bp = sum(lengths)
max_contig = max(lengths) if lengths else 0
n50_value = n50(lengths)

with stats_file.open("w") as stats:
    print(
        "sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "contigs",
        "total_bp",
        "max_contig_bp",
        "n50_bp",
        "renamed_fasta",
        sep="\\t",
        file=stats
    )

    print(
        sample_id,
        safe_id,
        assembly_sample_id,
        assembler,
        mode,
        rarefaction_label,
        assembly_strategy,
        contig_count,
        total_bp,
        max_contig,
        n50_value,
        published_fasta,
        sep="\\t",
        file=stats
    )

with manifest_record.open("w") as manifest:
    print(
        sample_id,
        safe_id,
        assembly_sample_id,
        assembler,
        mode,
        rarefaction_label,
        assembly_strategy,
        published_fasta,
        sep="\\t",
        file=manifest
    )
PY

    rm -f input_R1.fastq.gz input_R2.fastq.gz input_interleaved.fastq.gz
    rm -rf "\$RAW_OUT"
    """
}


process ASSEMBLE_RAREFIED {

    tag { "${sample_id}:${assembler}:rarefied:${rare_letter}" }

    publishDir "${params.outdir}/assemblies",
        mode: params.publish_assemblies_mode,
        pattern: "*.renamed.fa"

    publishDir "${params.outdir}/header_maps",
        mode: 'copy',
        pattern: "*.header_map.tsv"

    publishDir "${params.outdir}/logs",
        mode: 'copy',
        pattern: "*.log"

    publishDir "${params.outdir}/summary/per_assembly_stats",
        mode: 'copy',
        pattern: "*.assembly_stats.tsv"

    cpus {
        if( assembler == 'megahit' && params.megahit_threads != null ) {
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
    tuple val(sample_id),
          val(safe_id),
          val(assembly_sample_id),
          val(layout),
          val(read1),
          val(read2),
          val(interleaved),
          val(assembler),
          val(rare_index),
          val(rare_letter),
          val(rare_split_count),
          path(tools_status)

    output:
    path "*.renamed.fa", emit: renamed_contigs
    path "*.header_map.tsv", emit: header_map
    path "*.assembly_stats.tsv", emit: stats_file
    path "*.assembly_manifest_record.tsv", emit: manifest_record
    path "*.log", emit: log_file

    script:
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

    MODE="rarefied"
    RARE_INDEX="${rare_index}"
    RARE_ZERO_INDEX=\$((RARE_INDEX - 1))
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

    echo "Creating rarefied subset \${RAREFACTION_LABEL} of \${RARE_SPLIT_COUNT} for \${SAMPLE_ID}" > "\$LOG_FILE"

    python3 - \\
        "\$LAYOUT" \\
        "${read1}" \\
        "${read2}" \\
        "${interleaved}" \\
        "\$RARE_ZERO_INDEX" \\
        "\$RARE_SPLIT_COUNT" \\
        "\$SUB_R1" \\
        "\$SUB_R2" \\
        "\$SUB_12" <<'PY'
import gzip
import sys
from pathlib import Path

layout, read1, read2, interleaved, split_idx, split_count, out_r1, out_r2, out_12 = sys.argv[1:]

split_idx = int(split_idx)
split_count = int(split_count)

def open_fastq(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")

def read_fastq_records(path):
    with open_fastq(path) as handle:
        while True:
            h = handle.readline()
            if not h:
                break

            s = handle.readline()
            p = handle.readline()
            q = handle.readline()

            if not q:
                raise RuntimeError(f"Incomplete FASTQ record in {path}")

            yield h, s, p, q

records_written = 0

if layout == "paired":
    if not Path(read1).exists():
        raise RuntimeError(f"Read 1 file does not exist: {read1}")
    if not Path(read2).exists():
        raise RuntimeError(f"Read 2 file does not exist: {read2}")

    r1_iter = read_fastq_records(read1)
    r2_iter = read_fastq_records(read2)

    with gzip.open(out_r1, "wt") as o1, gzip.open(out_r2, "wt") as o2:
        for idx, (rec1, rec2) in enumerate(zip(r1_iter, r2_iter)):
            if idx % split_count == split_idx:
                o1.writelines(rec1)
                o2.writelines(rec2)
                records_written += 1

elif layout == "interleaved":
    if not Path(interleaved).exists():
        raise RuntimeError(f"Interleaved file does not exist: {interleaved}")

    rec_iter = read_fastq_records(interleaved)

    with gzip.open(out_12, "wt") as out:
        pair_idx = 0

        while True:
            try:
                rec1 = next(rec_iter)
            except StopIteration:
                break

            try:
                rec2 = next(rec_iter)
            except StopIteration:
                raise RuntimeError(f"Interleaved FASTQ has odd number of records: {interleaved}")

            if pair_idx % split_count == split_idx:
                out.writelines(rec1)
                out.writelines(rec2)
                records_written += 1

            pair_idx += 1

else:
    raise RuntimeError(f"Unsupported layout: {layout}")

if records_written == 0:
    raise RuntimeError(
        f"Rarefied subset {split_idx + 1} of {split_count} contains zero pairs/fragments."
    )
PY

    echo "Rarefied subset created successfully." >> "\$LOG_FILE"

    if [[ "\$ASSEMBLER" == "megahit" ]]; then

        if ! command -v megahit >/dev/null 2>&1; then
            echo "ERROR: megahit is not available after tool setup." >&2
            cat "${tools_status}" >&2 || true
            exit 1
        fi

        ASSEMBLY_STRATEGY="C"
        RAW_OUT="megahit_rarefied_out"

        MEGAHIT_MEM_ARG=""
        if [[ "${params.memory_gb}" != "0" ]]; then
            MEGAHIT_MEM_BYTES=\$(( ${params.memory_gb} * 1024 * 1024 * 1024 ))
            MEGAHIT_MEM_ARG="-m \$MEGAHIT_MEM_BYTES"
        fi

        echo "Running MEGAHIT rarefied assembly for \${SAMPLE_ID}, subset \${RAREFACTION_LABEL}" >> "\$LOG_FILE"
        echo "MEGAHIT threads: ${task.cpus}" >> "\$LOG_FILE"
        echo "Global memory GB: ${params.memory_gb}" >> "\$LOG_FILE"
        echo "MEGAHIT memory arg: \$MEGAHIT_MEM_ARG" >> "\$LOG_FILE"

        if [[ "\$LAYOUT" == "paired" ]]; then
            megahit \\
                -1 "\$SUB_R1" \\
                -2 "\$SUB_R2" \\
                -t ${task.cpus} \\
                \$MEGAHIT_MEM_ARG \\
                -o "\$RAW_OUT" \\
                --presets ${params.megahit_preset} \\
                >> "\$LOG_FILE" 2>&1
        else
            megahit \\
                --12 "\$SUB_12" \\
                -t ${task.cpus} \\
                \$MEGAHIT_MEM_ARG \\
                -o "\$RAW_OUT" \\
                --presets ${params.megahit_preset} \\
                >> "\$LOG_FILE" 2>&1
        fi

        SRC_FASTA="\$RAW_OUT/final.contigs.fa"

        if [[ ! -s "\$SRC_FASTA" ]]; then
            echo "ERROR: MEGAHIT did not produce final.contigs.fa" >> "\$LOG_FILE"
            exit 1
        fi

    elif [[ "\$ASSEMBLER" == "metaspades" ]]; then

        if ! command -v metaspades.py >/dev/null 2>&1; then
            echo "ERROR: metaspades.py is not available after tool setup." >&2
            cat "${tools_status}" >&2 || true
            exit 1
        fi

        ASSEMBLY_STRATEGY="D"
        RAW_OUT="metaspades_rarefied_out"

        SPADES_MEM_ARG=""
        if [[ "${params.memory_gb}" != "0" ]]; then
            SPADES_MEM_ARG="-m ${params.memory_gb}"
        fi

        echo "Running metaSPAdes rarefied assembly for \${SAMPLE_ID}, subset \${RAREFACTION_LABEL}" >> "\$LOG_FILE"
        echo "metaSPAdes threads: ${task.cpus}" >> "\$LOG_FILE"
        echo "Global memory GB: ${params.memory_gb}" >> "\$LOG_FILE"
        echo "metaSPAdes memory arg: \$SPADES_MEM_ARG" >> "\$LOG_FILE"

        if [[ "\$LAYOUT" == "paired" ]]; then
            metaspades.py \\
                -1 "\$SUB_R1" \\
                -2 "\$SUB_R2" \\
                -t ${task.cpus} \\
                \$SPADES_MEM_ARG \\
                -o "\$RAW_OUT" \\
                >> "\$LOG_FILE" 2>&1
        else
            metaspades.py \\
                --12 "\$SUB_12" \\
                -t ${task.cpus} \\
                \$SPADES_MEM_ARG \\
                -o "\$RAW_OUT" \\
                >> "\$LOG_FILE" 2>&1
        fi

        if [[ -s "\$RAW_OUT/scaffolds.fasta" ]]; then
            SRC_FASTA="\$RAW_OUT/scaffolds.fasta"
        elif [[ -s "\$RAW_OUT/contigs.fasta" ]]; then
            SRC_FASTA="\$RAW_OUT/contigs.fasta"
        else
            echo "ERROR: metaSPAdes did not produce scaffolds.fasta or contigs.fasta" >> "\$LOG_FILE"
            exit 1
        fi

    else
        echo "ERROR: Unsupported assembler: \$ASSEMBLER" >&2
        exit 1
    fi

    python3 - \\
        "\$SRC_FASTA" \\
        "\$OUT_FASTA" \\
        "\$HEADER_MAP" \\
        "\$STATS_FILE" \\
        "\$MANIFEST_RECORD" \\
        "\$SAMPLE_ID" \\
        "\$SAFE_ID" \\
        "\$ASSEMBLY_SAMPLE_ID" \\
        "\$ASSEMBLER" \\
        "\$MODE" \\
        "\$RAREFACTION_LABEL" \\
        "\$ASSEMBLY_STRATEGY" \\
        "${params.outdir}/assemblies/\$OUT_FASTA" <<'PY'
import re
import sys
from pathlib import Path

(
    src_fasta,
    out_fasta,
    header_map,
    stats_file,
    manifest_record,
    sample_id,
    safe_id,
    assembly_sample_id,
    assembler,
    mode,
    rarefaction_label,
    assembly_strategy,
    published_fasta
) = sys.argv[1:]

src_fasta = Path(src_fasta)
out_fasta = Path(out_fasta)
header_map = Path(header_map)
stats_file = Path(stats_file)
manifest_record = Path(manifest_record)

lengths = []
contig_count = 0
current_len = 0

def n50(vals):
    if not vals:
        return 0
    vals = sorted(vals, reverse=True)
    half = sum(vals) / 2
    running = 0
    for v in vals:
        running += v
        if running >= half:
            return v
    return 0

with src_fasta.open() as inp, out_fasta.open("w") as out, header_map.open("w") as hmap:
    print("old_header", "new_header", sep="\\t", file=hmap)

    for line in inp:
        line = line.rstrip("\\n")

        if line.startswith(">"):
            if contig_count > 0:
                lengths.append(current_len)

            current_len = 0
            contig_count += 1

            old_header = line[1:].strip()
            first_token = old_header.split()[0] if old_header else f"contig_{contig_count}"

            if assembler == "megahit":
                m = re.search(r'(k\\d+)_(\\d+)', first_token)
                if m:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_{m.group(1)}_{m.group(2)}"
                else:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_k000_{contig_count}"

            elif assembler == "metaspades":
                m = re.search(r'NODE_(\\d+)', first_token)
                if m:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_NODE_{m.group(1)}"
                else:
                    new_header = f"{assembly_sample_id}_{assembly_strategy}_NODE_{contig_count}"

            else:
                raise RuntimeError(f"Unsupported assembler: {assembler}")

            print(f">{new_header}", file=out)
            print(old_header, new_header, sep="\\t", file=hmap)

        else:
            seq = line.strip()
            current_len += len(seq)
            print(seq, file=out)

if contig_count > 0:
    lengths.append(current_len)

total_bp = sum(lengths)
max_contig = max(lengths) if lengths else 0
n50_value = n50(lengths)

with stats_file.open("w") as stats:
    print(
        "sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "contigs",
        "total_bp",
        "max_contig_bp",
        "n50_bp",
        "renamed_fasta",
        sep="\\t",
        file=stats
    )

    print(
        sample_id,
        safe_id,
        assembly_sample_id,
        assembler,
        mode,
        rarefaction_label,
        assembly_strategy,
        contig_count,
        total_bp,
        max_contig,
        n50_value,
        published_fasta,
        sep="\\t",
        file=stats
    )

with manifest_record.open("w") as manifest:
    print(
        sample_id,
        safe_id,
        assembly_sample_id,
        assembler,
        mode,
        rarefaction_label,
        assembly_strategy,
        published_fasta,
        sep="\\t",
        file=manifest
    )
PY

    rm -f "\$SUB_R1" "\$SUB_R2" "\$SUB_12"
    rm -rf "\$RAW_OUT"
    """
}


process WRITE_ASSEMBLY_SUMMARIES {

    tag "write_assembly_summaries"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "assembly_manifest.tsv"

    publishDir "${params.outdir}/summary",
        mode: 'copy',
        pattern: "assembly_stats_summary.tsv"

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
