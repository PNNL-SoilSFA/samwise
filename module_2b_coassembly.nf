#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Module 2b: MEGAHIT co-assembly from Module 1 trimmed reads.
 */

params.samwise_dir = params.samwise_dir ?: params.working_dir ?: projectDir
params.working_dir = params.working_dir ?: params.samwise_dir
params.output_dir = null
params.input_manifest = null
params.coassembly_groups = null
params.megahit_version = "1.2.9"
params.auto_install = true
params.tool_env_dir = null
params.threads = null
params.assembly_threads = 4
params.megahit_threads = null
params.memory_gb = 0
params.megahit_preset = "meta-large"
params.publish_assemblies_mode = "symlink"
params.publish_coassembly_reads_mode = "symlink"
params.results_dir = params.working_dir ? params.working_dir : (params.output_dir ? params.output_dir : ".")
params.module1_outdir = "${params.results_dir}/module_1_readtrimming"
params.outdir = "${params.results_dir}/module_2b_coassembly"

def absPath(value) {
    return java.nio.file.Paths
        .get(value.toString())
        .toAbsolutePath()
        .normalize()
        .toString()
}

workflow {
    if (!params.coassembly_groups) {
        error(
            """
        Missing required parameter: --coassembly_groups

        Expected two tab-separated columns:
          read_or_sample_id    group_id

        Example:
          sampleA    group_1
          sampleB    group_1
          sampleC    group_2
          sampleD    group_2

        Example command:
          nextflow run module_2b_coassembly.nf \\
            --working_dir ./results \\
            --coassembly_groups coassembly_groups.tsv
        """.stripIndent()
        )
    }

    def manifest_file = params.input_manifest ?: "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    log.info("Module 2b results directory: ${params.results_dir}")
    log.info("Using Module 1 trimmed manifest: ${manifest_file}")
    log.info("Using co-assembly groups table: ${params.coassembly_groups}")
    log.info("Writing Module 2b outputs to: ${params.outdir}")
    log.info("Assembler: MEGAHIT only")
    log.info("Assembly strategy letter: G")
    log.info("Threads: ${params.threads ?: params.assembly_threads}")
    log.info("MEGAHIT thread override: ${params.megahit_threads ?: 'not supplied'}")
    log.info("Global memory: ${(params.memory_gb as int) > 0 ? params.memory_gb + ' GB' : 'not supplied'}")

    def manifest_ch = channel.fromPath(
        manifest_file,
        type: 'file',
        checkIfExists: true,
    )

    def coassembly_groups_ch = channel.fromPath(
        params.coassembly_groups,
        type: 'file',
        checkIfExists: true,
    )

    SETUP_MODULE2B_TOOLS()

    PREPARE_COASSEMBLY_GROUPS(
        manifest_ch,
        coassembly_groups_ch,
    )

    def coassembly_jobs_ch = PREPARE_COASSEMBLY_GROUPS.out.group_files
        .flatten()
        .map { group_file ->
            def safe_group_id = group_file.name.replaceFirst(/\.coassembly_reads\.tsv$/, '')
            tuple(
                safe_group_id,
                safe_group_id,
                group_file,
            )
        }

    ASSEMBLE_COASSEMBLY(
        coassembly_jobs_ch.combine(SETUP_MODULE2B_TOOLS.out.status)
    )

    WRITE_ASSEMBLY_SUMMARY(
        ASSEMBLE_COASSEMBLY.out.manifest_record.collect(),
        ASSEMBLE_COASSEMBLY.out.stats_file.collect(),
    )

    WRITE_COASSEMBLY_TRIMMED_MANIFEST(
        ASSEMBLE_COASSEMBLY.out.coassembly_trimmed_manifest_record.collect()
    )
}

process SETUP_MODULE2B_TOOLS {
    tag "setup_megahit"
    cache false

    publishDir "${params.outdir}/setup", mode: 'copy', pattern: "module2b_tools_status.env"

    output:
    path "module2b_tools_status.env", emit: status

    script:
    def env_dir = params.tool_env_dir ?: "${params.outdir}/conda_envs/module2b_tools"

    """
    set -euo pipefail

    STATUS_FILE="module2b_tools_status.env"
    TOOL_ENV="${env_dir}"

    echo "Module 2b tool setup started: \$(date)" > "\$STATUS_FILE"
    echo "Requested MEGAHIT version: ${params.megahit_version}" >> "\$STATUS_FILE"
    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "----------------------------------------" >> "\$STATUS_FILE"

    check_megahit() {
        local prefix="\$1"

        if [[ "\$prefix" == "SYSTEM" ]]; then
            if ! command -v megahit >/dev/null 2>&1; then
                return 1
            fi
            megahit --version >> "\$STATUS_FILE" 2>&1 || return 1
            return 0
        fi

        if [[ ! -x "\$prefix/bin/megahit" ]]; then
            return 1
        fi

        "\$prefix/bin/megahit" --version >> "\$STATUS_FILE" 2>&1 || return 1
        return 0
    }

    if [[ -d "\$TOOL_ENV" ]]; then
        echo "Existing Module 2b environment detected: \$TOOL_ENV" >> "\$STATUS_FILE"

        if check_megahit "\$TOOL_ENV"; then
            echo "Existing Module 2b environment passed MEGAHIT check." >> "\$STATUS_FILE"
            echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
            echo "Module 2b tool setup finished: \$(date)" >> "\$STATUS_FILE"
            exit 0
        else
            echo "Existing Module 2b environment is incomplete or broken. Removing." >> "\$STATUS_FILE"
            rm -rf "\$TOOL_ENV"
        fi
    fi

    echo "Checking system/runtime MEGAHIT..." >> "\$STATUS_FILE"

    if check_megahit "SYSTEM"; then
        echo "MEGAHIT is available from system/runtime PATH." >> "\$STATUS_FILE"
        echo "TOOL_ENV=SYSTEM" >> "\$STATUS_FILE"
        echo "Module 2b tool setup finished: \$(date)" >> "\$STATUS_FILE"
        exit 0
    fi

    echo "MEGAHIT is unavailable from system/runtime PATH." >> "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then
        echo "ERROR: auto_install is false and MEGAHIT is missing." >> "\$STATUS_FILE"
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

    echo "Creating Module 2b environment: \$TOOL_ENV" >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \\
        -p "\$TOOL_ENV" \\
        -c conda-forge \\
        -c bioconda \\
        "python" \\
        "megahit=${params.megahit_version}" \\
        >> "\$STATUS_FILE" 2>&1

    if ! check_megahit "\$TOOL_ENV"; then
        echo "ERROR: Newly created Module 2b environment failed MEGAHIT check." >> "\$STATUS_FILE"
        exit 1
    fi

    echo "TOOL_ENV=\$TOOL_ENV" >> "\$STATUS_FILE"
    echo "Module 2b tool setup finished: \$(date)" >> "\$STATUS_FILE"
    """
}

process PREPARE_COASSEMBLY_GROUPS {
    tag "prepare_coassembly_groups"

    publishDir "${params.outdir}/coassembly/groups", mode: 'copy', pattern: "*.coassembly_reads.tsv"
    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "coassembly_group_summary.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "prepare_coassembly_groups.log"

    input:
    path trimmed_manifest
    path group_table

    output:
    path "*.coassembly_reads.tsv", emit: group_files
    path "coassembly_group_summary.tsv", emit: summary
    path "prepare_coassembly_groups.log", emit: log_file

    script:
    """
    set -euo pipefail

    LOG_FILE="prepare_coassembly_groups.log"

    echo "Preparing co-assembly groups: \$(date)" > "\$LOG_FILE"
    echo "Trimmed manifest: ${trimmed_manifest}" >> "\$LOG_FILE"
    echo "Group table: ${group_table}" >> "\$LOG_FILE"
    echo "Nextflow launchDir: ${workflow.launchDir}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    python3 - \\
        "${trimmed_manifest}" \\
        "${group_table}" \\
        "coassembly_group_summary.tsv" \\
        "\$LOG_FILE" \\
        "${workflow.launchDir}" <<'PY'
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

trimmed_manifest = Path(sys.argv[1])
group_table = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
log_path = Path(sys.argv[4])
launch_dir = Path(sys.argv[5]).resolve()

def log(msg):
    with log_path.open("a") as handle:
        print(msg, file=handle)

def safe_id(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("_")
    return value or "unnamed_group"

def strip_fastq_ext(name):
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name

def resolve_path(value):
    value = str(value or "").strip()
    if not value:
        return ""
    p = Path(value)
    if p.is_absolute():
        return str(p.resolve())
    return str((launch_dir / p).resolve())

def path_keys(value):
    value = str(value or "").strip()
    if not value:
        return set()

    resolved = resolve_path(value)
    p_original = Path(value)
    p_resolved = Path(resolved)

    keys = set()

    keys.add(value)
    keys.add(resolved)

    keys.add(p_original.name)
    keys.add(strip_fastq_ext(p_original.name))

    keys.add(p_resolved.name)
    keys.add(strip_fastq_ext(p_resolved.name))

    return keys

if not trimmed_manifest.exists():
    raise SystemExit(f"ERROR: trimmed manifest does not exist: {trimmed_manifest}")

if not group_table.exists():
    raise SystemExit(f"ERROR: coassembly group table does not exist: {group_table}")

sample_rows = []
lookup = {}

with trimmed_manifest.open() as handle:
    reader = csv.DictReader(handle, delimiter="\\t")

    required = {
        "sample_id",
        "safe_sample_id",
        "layout",
        "read1",
        "read2",
        "interleaved",
    }

    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(
            f"ERROR: trimmed manifest is missing required columns: {sorted(missing)}"
        )

    for row in reader:
        sample_id = row["sample_id"].strip()
        safe_sample_id = row["safe_sample_id"].strip()
        layout = row["layout"].strip()

        if layout not in {"paired", "interleaved"}:
            raise SystemExit(
                f"ERROR: unsupported layout in trimmed manifest for sample {sample_id}: {layout}"
            )

        row["read1"] = resolve_path(row.get("read1", ""))
        row["read2"] = resolve_path(row.get("read2", ""))
        row["interleaved"] = resolve_path(row.get("interleaved", ""))

        sample_rows.append(row)

        keys = set()
        keys.add(sample_id)
        keys.add(safe_sample_id)
        keys.update(path_keys(row.get("read1", "")))
        keys.update(path_keys(row.get("read2", "")))
        keys.update(path_keys(row.get("interleaved", "")))

        for key in keys:
            if not key:
                continue

            if key in lookup and lookup[key]["sample_id"] != sample_id:
                raise SystemExit(
                    f"ERROR: coassembly lookup key is ambiguous: {key}\\n"
                    f"  First sample: {lookup[key]['sample_id']}\\n"
                    f"  Second sample: {sample_id}"
                )

            lookup[key] = row

log(f"Samples in trimmed manifest: {len(sample_rows)}")
log(f"Lookup keys generated: {len(lookup)}")

assignments_by_sample = {}
unmatched = []

with group_table.open() as handle:
    for line_number, line in enumerate(handle, start=1):
        line = line.rstrip("\\n")

        if not line.strip():
            continue

        fields = line.split("\\t")

        if len(fields) < 2:
            raise SystemExit(
                f"ERROR: group table line {line_number} has fewer than 2 tab-separated columns: {line}"
            )

        read_or_sample_id = fields[0].strip()
        group_id = fields[1].strip()

        if line_number == 1:
            left = read_or_sample_id.lower()
            right = group_id.lower()

            if left in {
                "read",
                "read_id",
                "read_or_sample_id",
                "fastq",
                "fastq_id",
                "fastq_file",
                "fastq_file_id",
                "sample",
                "sample_id",
            } and right in {
                "group",
                "group_id",
                "coassembly_group",
            }:
                log("Detected and skipped header row in coassembly group table.")
                continue

        if not read_or_sample_id or not group_id:
            raise SystemExit(
                f"ERROR: empty read/sample ID or group ID at group table line {line_number}"
            )

        possible_keys = {
            read_or_sample_id,
            resolve_path(read_or_sample_id),
            Path(read_or_sample_id).name,
            strip_fastq_ext(Path(read_or_sample_id).name),
        }

        matched_row = None

        for key in possible_keys:
            if key in lookup:
                matched_row = lookup[key]
                break

        if matched_row is None:
            unmatched.append((line_number, read_or_sample_id, group_id))
            continue

        sample_id = matched_row["sample_id"].strip()

        if sample_id in assignments_by_sample:
            previous_group = assignments_by_sample[sample_id]

            if previous_group != group_id:
                raise SystemExit(
                    f"ERROR: sample {sample_id} was assigned to multiple coassembly groups:\\n"
                    f"  {previous_group}\\n"
                    f"  {group_id}"
                )

        assignments_by_sample[sample_id] = group_id

if unmatched:
    log("WARNING: Some coassembly group table rows did not match the trimmed manifest.")
    for line_number, read_id, group_id in unmatched[:50]:
        log(f"  unmatched line {line_number}: {read_id}\\t{group_id}")
    if len(unmatched) > 50:
        log(f"  ... {len(unmatched) - 50} more unmatched rows not shown")

if not assignments_by_sample:
    raise SystemExit(
        "ERROR: no valid coassembly assignments were found. "
        "Check that column 1 of the group table matches sample_id, safe_sample_id, "
        "or FASTQ filenames from trimmed_manifest.tsv."
    )

groups = defaultdict(list)

for row in sample_rows:
    sample_id = row["sample_id"].strip()

    if sample_id in assignments_by_sample:
        groups[assignments_by_sample[sample_id]].append(row)

if not groups:
    raise SystemExit("ERROR: no coassembly groups were created.")

used_safe_group_ids = set()

with summary_path.open("w") as summary:
    print(
        "group_id",
        "safe_group_id",
        "sample_count",
        "paired_sample_count",
        "interleaved_sample_count",
        "group_reads_tsv",
        sep="\\t",
        file=summary,
    )

    for group_id in sorted(groups):
        rows = groups[group_id]

        safe_group_id_base = safe_id(group_id)
        safe_group_id = safe_group_id_base
        suffix = 1

        while safe_group_id in used_safe_group_ids:
            suffix += 1
            safe_group_id = f"{safe_group_id_base}_{suffix}"

        used_safe_group_ids.add(safe_group_id)

        out_path = Path(f"{safe_group_id}.coassembly_reads.tsv")

        paired_count = sum(1 for r in rows if r["layout"].strip() == "paired")
        interleaved_count = sum(1 for r in rows if r["layout"].strip() == "interleaved")

        with out_path.open("w") as out:
            print(
                "group_id",
                "safe_group_id",
                "sample_id",
                "safe_sample_id",
                "layout",
                "read1",
                "read2",
                "interleaved",
                sep="\\t",
                file=out,
            )

            for row in rows:
                print(
                    group_id,
                    safe_group_id,
                    row["sample_id"].strip(),
                    row["safe_sample_id"].strip(),
                    row["layout"].strip(),
                    row.get("read1", "").strip(),
                    row.get("read2", "").strip(),
                    row.get("interleaved", "").strip(),
                    sep="\\t",
                    file=out,
                )

        print(
            group_id,
            safe_group_id,
            len(rows),
            paired_count,
            interleaved_count,
            str(out_path),
            sep="\\t",
            file=summary,
        )

        log(
            f"Prepared coassembly group {group_id} as {safe_group_id}: "
            f"{len(rows)} samples, {paired_count} paired, {interleaved_count} interleaved"
        )

log("Coassembly group preparation completed.")
PY

    echo "Co-assembly group preparation finished: \$(date)" >> "\$LOG_FILE"
    """
}

process ASSEMBLE_COASSEMBLY {
    tag { "${safe_group_id}:megahit:coassembly" }

    stageInMode 'symlink'

    publishDir "${params.outdir}/assemblies", mode: params.publish_assemblies_mode, pattern: "*.renamed.fa"
    publishDir "${params.outdir}/header_maps", mode: 'copy', pattern: "*.header_map.tsv"
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: "*.log"
    publishDir "${params.outdir}/summary/per_assembly_stats", mode: 'copy', pattern: "*.assembly_stats.tsv"
    publishDir "${params.outdir}/coassembly/interleaved_inputs", mode: params.publish_coassembly_reads_mode, pattern: "*.coassembly_interleaved.fastq.gz"

    cpus {
        if (params.megahit_threads != null) {
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
    tuple val(group_id), val(safe_group_id), path(group_reads_tsv), path(tools_status)

    output:
    path "*.renamed.fa", emit: renamed_contigs
    path "*.header_map.tsv", emit: header_map
    path "*.assembly_stats.tsv", emit: stats_file
    path "*.assembly_manifest_record.tsv", emit: manifest_record
    path "*_coassembly_trimmed_manifest_record.tsv", emit: coassembly_trimmed_manifest_record
    path "*.coassembly_interleaved.fastq.gz", emit: coassembly_interleaved_reads
    path "*.log", emit: log_file

    script:
    def outdir_abs = absPath(params.outdir)

    """
    set -euo pipefail

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" | tail -n 1 | cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    if ! command -v megahit >/dev/null 2>&1; then
        echo "ERROR: megahit is not available after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    GROUP_ID="${group_id}"
    SAFE_GROUP_ID="${safe_group_id}"
    ASSEMBLER="megahit"
    MODE="coassembly"
    ASSEMBLY_STRATEGY="G"
    ASSEMBLY_STATUS="ok"
    ASSEMBLY_WARNING=""

    ASSEMBLY_SAMPLE_ID="\$(python3 - <<'PY'
import re
value = "${group_id}"
clean = re.sub(r"[^A-Za-z0-9]+", "", value)
if not clean:
    clean = re.sub(r"[^A-Za-z0-9]+", "", "${safe_group_id}")
print(clean or "coassembly")
PY
)"

    OUT_FASTA="\${SAFE_GROUP_ID}_megahit_coassembly.renamed.fa"
    HEADER_MAP="\${SAFE_GROUP_ID}_megahit_coassembly.header_map.tsv"
    STATS_FILE="\${SAFE_GROUP_ID}_megahit_coassembly.assembly_stats.tsv"
    MANIFEST_RECORD="\${SAFE_GROUP_ID}_megahit_coassembly.assembly_manifest_record.tsv"
    COASSEMBLY_TRIMMED_RECORD="\${SAFE_GROUP_ID}_coassembly_trimmed_manifest_record.tsv"
    LOG_FILE="\${SAFE_GROUP_ID}_megahit_coassembly.log"
    COASSEMBLY_FASTQ="\${SAFE_GROUP_ID}.coassembly_interleaved.fastq.gz"

    echo "Running MEGAHIT co-assembly for group: \${GROUP_ID}" > "\$LOG_FILE"
    echo "Safe group ID: \${SAFE_GROUP_ID}" >> "\$LOG_FILE"
    echo "Assembly sample ID: \${ASSEMBLY_SAMPLE_ID}" >> "\$LOG_FILE"
    echo "Group reads TSV: ${group_reads_tsv}" >> "\$LOG_FILE"
    echo "Threads: ${task.cpus}" >> "\$LOG_FILE"
    echo "Global memory GB: ${params.memory_gb}" >> "\$LOG_FILE"
    echo "MEGAHIT preset: ${params.megahit_preset}" >> "\$LOG_FILE"
    echo "----------------------------------------" >> "\$LOG_FILE"

    python3 - \\
        "${group_reads_tsv}" \\
        "\$COASSEMBLY_FASTQ" \\
        "\$LOG_FILE" <<'PY'
import csv
import gzip
import sys
from itertools import zip_longest
from pathlib import Path

group_reads_tsv = Path(sys.argv[1])
out_fastq = Path(sys.argv[2])
log_file = Path(sys.argv[3])

def log(msg):
    with log_file.open("a") as handle:
        print(msg, file=handle)

def open_fastq(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")

def iter_fastq_records(path):
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

if not group_reads_tsv.exists():
    raise SystemExit(f"ERROR: group reads TSV does not exist: {group_reads_tsv}")

sample_count = 0
paired_count = 0
interleaved_count = 0
records_written = 0

with group_reads_tsv.open() as inp, gzip.open(out_fastq, "wt", compresslevel=4) as out:
    reader = csv.DictReader(inp, delimiter="\\t")

    required = {
        "group_id",
        "safe_group_id",
        "sample_id",
        "safe_sample_id",
        "layout",
        "read1",
        "read2",
        "interleaved",
    }

    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"ERROR: group reads TSV missing required columns: {sorted(missing)}")

    for row in reader:
        sample_count += 1
        sample_id = row["sample_id"].strip()
        layout = row["layout"].strip()

        if layout == "paired":
            paired_count += 1

            read1 = row["read1"].strip()
            read2 = row["read2"].strip()

            if not read1 or not read2:
                raise SystemExit(f"ERROR: paired sample {sample_id} is missing read1/read2 paths.")

            if not Path(read1).exists():
                raise SystemExit(f"ERROR: read1 does not exist for sample {sample_id}: {read1}")

            if not Path(read2).exists():
                raise SystemExit(f"ERROR: read2 does not exist for sample {sample_id}: {read2}")

            r1_iter = iter_fastq_records(read1)
            r2_iter = iter_fastq_records(read2)

            pair_count = 0

            for rec1, rec2 in zip_longest(r1_iter, r2_iter):
                if rec1 is None or rec2 is None:
                    raise SystemExit(
                        f"ERROR: paired FASTQ files have unequal record counts for sample {sample_id}"
                    )

                out.writelines(rec1)
                out.writelines(rec2)

                records_written += 2
                pair_count += 1

            log(f"Added paired sample {sample_id}: {pair_count} pairs")

        elif layout == "interleaved":
            interleaved_count += 1

            interleaved = row["interleaved"].strip()

            if not interleaved:
                raise SystemExit(f"ERROR: interleaved sample {sample_id} is missing interleaved path.")

            if not Path(interleaved).exists():
                raise SystemExit(f"ERROR: interleaved FASTQ does not exist for sample {sample_id}: {interleaved}")

            rec_count = 0

            for rec in iter_fastq_records(interleaved):
                out.writelines(rec)
                records_written += 1
                rec_count += 1

            if rec_count % 2 != 0:
                raise SystemExit(
                    f"ERROR: interleaved sample {sample_id} has odd number of records: {rec_count}"
                )

            log(f"Added interleaved sample {sample_id}: {rec_count} records")

        else:
            raise SystemExit(f"ERROR: unsupported layout for sample {sample_id}: {layout}")

if sample_count == 0:
    raise SystemExit("ERROR: no samples were present in the coassembly group TSV.")

if records_written == 0:
    raise SystemExit("ERROR: no FASTQ records were written for coassembly.")

if records_written % 2 != 0:
    raise SystemExit(
        f"ERROR: coassembly interleaved FASTQ has odd number of records: {records_written}"
    )

log("----------------------------------------")
log(f"Coassembly samples: {sample_count}")
log(f"Paired samples: {paired_count}")
log(f"Interleaved samples: {interleaved_count}")
log(f"Interleaved FASTQ records written: {records_written}")
log(f"Interleaved FASTQ pairs/fragments: {records_written // 2}")
PY

    echo "Co-assembly interleaved FASTQ created:" >> "\$LOG_FILE"
    ls -lh "\$COASSEMBLY_FASTQ" >> "\$LOG_FILE" 2>&1 || true
    echo "----------------------------------------" >> "\$LOG_FILE"

    MEGAHIT_MEM_ARG=""

    if [[ "${params.memory_gb}" != "0" ]]; then
        MEGAHIT_MEM_BYTES=\$(( ${params.memory_gb} * 1024 * 1024 * 1024 ))
        MEGAHIT_MEM_ARG="-m \$MEGAHIT_MEM_BYTES"
    fi

    echo "MEGAHIT memory arg: \$MEGAHIT_MEM_ARG" >> "\$LOG_FILE"

    megahit \\
        --12 "\$COASSEMBLY_FASTQ" \\
        -t ${task.cpus} \\
        \$MEGAHIT_MEM_ARG \\
        -o megahit_coassembly_out \\
        --presets ${params.megahit_preset} \\
        >> "\$LOG_FILE" 2>&1

    SRC_FASTA="megahit_coassembly_out/final.contigs.fa"

    if [[ ! -s "\$SRC_FASTA" ]]; then
        echo "ERROR: MEGAHIT did not produce final.contigs.fa for coassembly group \${GROUP_ID}" >> "\$LOG_FILE"
        exit 1
    fi

    python3 - \\
        "\$SRC_FASTA" \\
        "\$OUT_FASTA" \\
        "\$HEADER_MAP" \\
        "\$STATS_FILE" \\
        "\$MANIFEST_RECORD" \\
        "\$COASSEMBLY_TRIMMED_RECORD" \\
        "\$GROUP_ID" \\
        "\$SAFE_GROUP_ID" \\
        "\$ASSEMBLY_SAMPLE_ID" \\
        "\$ASSEMBLER" \\
        "\$MODE" \\
        "\$ASSEMBLY_STRATEGY" \\
        "\$ASSEMBLY_STATUS" \\
        "\$ASSEMBLY_WARNING" \\
        "${outdir_abs}/assemblies/\$OUT_FASTA" \\
        "${outdir_abs}/coassembly/interleaved_inputs/\$COASSEMBLY_FASTQ" <<'PY'
import re
import sys
from pathlib import Path

(
    src_fasta,
    out_fasta,
    header_map,
    stats_file,
    manifest_record,
    coassembly_trimmed_record,
    group_id,
    safe_group_id,
    assembly_sample_id,
    assembler,
    mode,
    assembly_strategy,
    assembly_status,
    assembly_warning,
    published_fasta,
    published_coassembly_fastq,
) = sys.argv[1:]

src_fasta = Path(src_fasta)
out_fasta = Path(out_fasta)
header_map = Path(header_map)
stats_file = Path(stats_file)
manifest_record = Path(manifest_record)
coassembly_trimmed_record = Path(coassembly_trimmed_record)

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

            m = re.search(r'(k\\d+)_(\\d+)', first_token)

            if m:
                new_header = f"{assembly_sample_id}_{assembly_strategy}_{m.group(1)}_{m.group(2)}"
            else:
                new_header = f"{assembly_sample_id}_{assembly_strategy}_k000_{contig_count}"

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
        "assembly_status",
        "assembly_warning",
        "contigs",
        "total_bp",
        "max_contig_bp",
        "n50_bp",
        "renamed_fasta",
        sep="\\t",
        file=stats,
    )

    print(
        group_id,
        safe_group_id,
        assembly_sample_id,
        assembler,
        mode,
        "",
        assembly_strategy,
        assembly_status,
        assembly_warning,
        contig_count,
        total_bp,
        max_contig,
        n50_value,
        published_fasta,
        sep="\\t",
        file=stats,
    )

with manifest_record.open("w") as manifest:
    print(
        group_id,
        safe_group_id,
        assembly_sample_id,
        assembler,
        mode,
        "",
        assembly_strategy,
        published_fasta,
        sep="\\t",
        file=manifest,
    )

with coassembly_trimmed_record.open("w") as out:
    print(
        group_id,
        safe_group_id,
        "interleaved",
        "",
        "",
        published_coassembly_fastq,
        "",
        "",
        "",
        sep="\\t",
        file=out,
    )
PY

    rm -rf megahit_coassembly_out

    echo "MEGAHIT co-assembly finished: \$(date)" >> "\$LOG_FILE"
    """
}

process WRITE_ASSEMBLY_SUMMARY {
    tag "write_coassembly_assembly_summary"

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
        echo "ERROR: No coassembly manifest records were received." >&2
        exit 1
    fi

    if [[ -z "${stats_file_list}" ]]; then
        echo "ERROR: No coassembly stats files were received." >&2
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

process WRITE_COASSEMBLY_TRIMMED_MANIFEST {
    tag "write_coassembly_trimmed_manifest"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: "coassembly_trimmed_manifest.tsv"

    input:
    path records

    output:
    path "coassembly_trimmed_manifest.tsv", emit: manifest

    script:
    def files = records.collect { record_file -> record_file.name }.join(' ')

    """
    set -euo pipefail

    printf 'sample_id\\tsafe_sample_id\\tlayout\\tread1\\tread2\\tinterleaved\\tmerged\\tfastp_html\\tfastp_json\\n' > coassembly_trimmed_manifest.tsv

    if [[ -z "${files}" ]]; then
        echo "ERROR: No coassembly trimmed-manifest records were received." >&2
        exit 1
    fi

    for f in ${files}; do
        if [[ -s "\$f" ]]; then
            cat "\$f" >> coassembly_trimmed_manifest.tsv
        fi
    done
    """
}
