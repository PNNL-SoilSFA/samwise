#!/usr/bin/env python3
"""Run MEGAHIT or metaSPAdes for one Module 2 assembly unit.

Standalone port of the assembler-invocation bash that was duplicated between
the ASSEMBLE_SINGLE and ASSEMBLE_RAREFIED processes of
module_2_readassembly.nf.

This is the unit of work Nextflow fans out: one invocation == one assembly ==
one SLURM job. It can be run by hand for debugging with no Nextflow involved.

On success it writes a shell-sourceable result file (--result-file, default
assembler_result.env) describing what happened:

    SRC_FASTA='megahit_out/final.contigs.fa'
    ASSEMBLY_STATUS='ok'
    ASSEMBLY_WARNING=''
    ASSEMBLER_EXIT_CODE='0'
    RAW_OUT='megahit_out'

metaSPAdes exit code 12 means "not enough memory". That is treated as
non-fatal so the run can still produce summaries, exactly as before: status
becomes 'failed_nonfatal' and an empty placeholder FASTA is emitted if
metaSPAdes wrote no contigs.

Run standalone:

  samwise_run_assembler.py \
      --assembler megahit --layout paired \
      --read1 R1.fq.gz --read2 R2.fq.gz \
      --threads 36 --memory-gb 0 \
      --out-dir megahit_out --log-file sample.log
"""

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

OOM_EXIT_CODE = 12
OOM_WARNING = "assembly failed - not enough memory"

# Fraction of detected node memory to hand to metaSPAdes, leaving headroom for
# the OS and for metaSPAdes' own accounting slop.
MEMORY_SAFETY_FRACTION = 0.9


def log(log_file, message):
    with open(log_file, "a") as handle:
        print(message, file=handle)


def detect_memory_gb():
    """Best-effort detection of memory available to this process, in whole GB.

    Returns None if nothing conclusive is found. Checked in order of
    specificity: SLURM's own accounting, then the cgroup limit, then physical
    RAM. Under `--exclusive --mem=0` SLURM reports 0, which we treat as unset
    and fall through to physical RAM -- the correct answer for a whole node.
    """
    mem_per_node = os.environ.get("SLURM_MEM_PER_NODE")

    if mem_per_node and mem_per_node.isdigit() and int(mem_per_node) > 0:
        return int(mem_per_node) // 1024

    mem_per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    cpus_on_node = os.environ.get("SLURM_CPUS_ON_NODE")

    if (
        mem_per_cpu and mem_per_cpu.isdigit() and int(mem_per_cpu) > 0
        and cpus_on_node and cpus_on_node.isdigit() and int(cpus_on_node) > 0
    ):
        return (int(mem_per_cpu) * int(cpus_on_node)) // 1024

    for cgroup_file in (
        "/sys/fs/cgroup/memory.max",
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
    ):
        try:
            raw = Path(cgroup_file).read_text().strip()
        except OSError:
            continue

        if raw.isdigit():
            limit = int(raw)

            # Unbounded cgroups report a sentinel close to 2**63.
            if 0 < limit < (1 << 62):
                return limit // (1024 ** 3)

    try:
        total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None

    if total > 0:
        return total // (1024 ** 3)

    return None


def resolve_memory_gb(requested_gb, log_file):
    """Return the memory budget in GB, or None to leave it to the assembler.

    A requested value of 0 means "auto": detect what the job actually has.
    Passing an explicit budget to metaSPAdes matters because its default cap is
    250 GB, which silently triggers exit 12 on smaller nodes.
    """
    if requested_gb > 0:
        return requested_gb

    detected = detect_memory_gb()

    if detected is None:
        log(log_file, "Memory auto-detection failed; leaving assembler default in place.")
        return None

    budget = int(detected * MEMORY_SAFETY_FRACTION)

    if budget < 1:
        log(log_file, f"Detected memory {detected} GB is too small to budget; using assembler default.")
        return None

    log(log_file, f"Auto-detected {detected} GB available; budgeting {budget} GB for the assembler.")

    return budget


def require_tool(tool, log_file):
    if shutil.which(tool) is None:
        message = f"ERROR: {tool} is not available on PATH after tool setup."
        log(log_file, message)
        print(message, file=sys.stderr)
        sys.exit(1)


def run_and_capture(cmd, log_file):
    """Run cmd, appending both stdout and stderr to log_file. Returns exit code."""
    log(log_file, f"Command: {' '.join(shlex.quote(c) for c in cmd)}")

    with open(log_file, "a") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)

    return proc.returncode


def build_read_args(layout, read1, read2, interleaved):
    if layout == "paired":
        return ["-1", read1, "-2", read2]

    return ["--12", interleaved]


def run_megahit(args, memory_gb, log_file):
    require_tool("megahit", log_file)

    cmd = ["megahit"]
    cmd += build_read_args(args.layout, args.read1, args.read2, args.interleaved)
    cmd += ["-t", str(args.threads)]

    if memory_gb is not None:
        cmd += ["-m", str(memory_gb * 1024 * 1024 * 1024)]

    cmd += ["-o", args.out_dir, "--presets", args.megahit_preset]

    log(log_file, f"MEGAHIT threads: {args.threads}")
    log(log_file, f"MEGAHIT memory budget GB: {memory_gb if memory_gb is not None else 'assembler default'}")

    exit_code = run_and_capture(cmd, log_file)

    if exit_code != 0:
        log(log_file, f"ERROR: MEGAHIT failed with exit code {exit_code}")
        sys.exit(exit_code)

    src_fasta = Path(args.out_dir) / "final.contigs.fa"

    if not src_fasta.exists() or src_fasta.stat().st_size == 0:
        log(log_file, "ERROR: MEGAHIT did not produce final.contigs.fa")
        sys.exit(1)

    return str(src_fasta), "ok", "", exit_code


def run_metaspades(args, memory_gb, log_file):
    require_tool("metaspades.py", log_file)

    cmd = ["metaspades.py"]
    cmd += build_read_args(args.layout, args.read1, args.read2, args.interleaved)
    cmd += ["-t", str(args.threads)]

    if memory_gb is not None:
        cmd += ["-m", str(memory_gb)]

    cmd += ["-o", args.out_dir]

    log(log_file, f"metaSPAdes threads: {args.threads}")
    log(log_file, f"metaSPAdes memory budget GB: {memory_gb if memory_gb is not None else 'assembler default'}")

    exit_code = run_and_capture(cmd, log_file)

    log(log_file, f"metaSPAdes exit code: {exit_code}")

    if exit_code not in (0, OOM_EXIT_CODE):
        log(log_file, f"ERROR: metaSPAdes failed with fatal exit code {exit_code}")
        sys.exit(exit_code)

    status = "ok"
    warning = ""

    if exit_code == OOM_EXIT_CODE:
        status = "failed_nonfatal"
        warning = OOM_WARNING
        log(log_file, "WARNING: metaSPAdes exited with code 12. Treating this as non-fatal so summaries can be written.")
        log(log_file, f"WARNING: {warning}")

    out_dir = Path(args.out_dir)

    for candidate in (out_dir / "scaffolds.fasta", out_dir / "contigs.fasta"):
        if candidate.exists() and candidate.stat().st_size > 0:
            return str(candidate), status, warning, exit_code

    if exit_code == OOM_EXIT_CODE:
        log(log_file, "WARNING: metaSPAdes exit code 12 produced no scaffolds.fasta or contigs.fasta. Creating empty placeholder FASTA.")
        placeholder = Path("metaspades_exit12_empty_contigs.fasta")
        placeholder.write_text("")
        return str(placeholder), status, warning, exit_code

    log(log_file, "ERROR: metaSPAdes did not produce scaffolds.fasta or contigs.fasta")
    sys.exit(1)


def write_result_file(path, src_fasta, status, warning, exit_code, raw_out):
    """Write a shell-sourceable result file. Values are quoted for `source`."""
    entries = (
        ("SRC_FASTA", src_fasta),
        ("ASSEMBLY_STATUS", status),
        ("ASSEMBLY_WARNING", warning),
        ("ASSEMBLER_EXIT_CODE", str(exit_code)),
        ("RAW_OUT", raw_out),
    )

    with open(path, "w") as handle:
        for key, value in entries:
            print(f"{key}={shlex.quote(value)}", file=handle)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run one MEGAHIT or metaSPAdes assembly for SAMWISE Module 2.",
    )

    parser.add_argument("--assembler", required=True, choices=["megahit", "metaspades"])
    parser.add_argument("--layout", required=True, choices=["paired", "interleaved"])

    parser.add_argument("--read1", default="", help="R1 FASTQ (paired layout).")
    parser.add_argument("--read2", default="", help="R2 FASTQ (paired layout).")
    parser.add_argument("--interleaved", default="", help="Interleaved FASTQ.")

    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--memory-gb", default=0, type=int,
                        help="Memory budget in GB. 0 means auto-detect from the "
                             "SLURM allocation / cgroup / physical RAM.")

    parser.add_argument("--out-dir", required=True,
                        help="Raw assembler output directory.")
    parser.add_argument("--log-file", required=True,
                        help="Log file; assembler stdout and stderr are appended here.")
    parser.add_argument("--result-file", default="assembler_result.env",
                        help="Shell-sourceable result file to write.")

    parser.add_argument("--megahit-preset", default="meta-large")
    parser.add_argument("--clean-stale-output", default="true",
                        help="Remove a pre-existing --out-dir before starting.")

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    log_file = args.log_file
    Path(log_file).touch()

    if args.layout == "paired":
        for label, path in (("read1", args.read1), ("read2", args.read2)):
            if not path or not Path(path).exists():
                message = f"ERROR: missing {label} for paired layout: {path!r}"
                log(log_file, message)
                print(message, file=sys.stderr)
                return 1
    else:
        if not args.interleaved or not Path(args.interleaved).exists():
            message = f"ERROR: missing interleaved reads: {args.interleaved!r}"
            log(log_file, message)
            print(message, file=sys.stderr)
            return 1

    if args.clean_stale_output.lower() == "true" and Path(args.out_dir).exists():
        log(log_file, f"Removing stale assembler output from previous attempt: {args.out_dir}")
        shutil.rmtree(args.out_dir, ignore_errors=True)

    memory_gb = resolve_memory_gb(args.memory_gb, log_file)

    if args.assembler == "megahit":
        src_fasta, status, warning, exit_code = run_megahit(args, memory_gb, log_file)
    else:
        src_fasta, status, warning, exit_code = run_metaspades(args, memory_gb, log_file)

    write_result_file(args.result_file, src_fasta, status, warning, exit_code, args.out_dir)

    log(log_file, f"Assembler finished. Source FASTA: {src_fasta} (status: {status})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
