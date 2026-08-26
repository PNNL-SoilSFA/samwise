#!/usr/bin/env python3
"""Rename assembly contigs and emit Module 2 stats / manifest records.

This is a standalone port of the Python that used to be embedded as a heredoc
inside the ASSEMBLE_SINGLE and ASSEMBLE_RAREFIED processes of
module_2_readassembly.nf. Both copies were byte-identical, so they are now a
single script.

Header rewriting rules (unchanged):

  megahit     >k141_1234 ...            -> >{assembly_sample_id}_{strategy}_k141_1234
              (no k<N>_<M> token)       -> >{assembly_sample_id}_{strategy}_k000_{n}
  metaspades  >NODE_12_length_...       -> >{assembly_sample_id}_{strategy}_NODE_12
              (no NODE_<N> token)       -> >{assembly_sample_id}_{strategy}_NODE_{n}

where {n} is the 1-based ordinal of the contig in the input file.

Run standalone:

  samwise_rename_contigs.py \
      --src-fasta megahit_out/final.contigs.fa \
      --out-fasta S1_megahit_single.renamed.fa \
      --header-map S1_megahit_single.header_map.tsv \
      --stats-file S1_megahit_single.assembly_stats.tsv \
      --manifest-record S1_megahit_single.assembly_manifest_record.tsv \
      --sample-id S1 --safe-id S1 --assembly-sample-id S1 \
      --assembler megahit --mode single --assembly-strategy A \
      --published-fasta /path/to/published/S1_megahit_single.renamed.fa
"""

import argparse
import re
import sys
from pathlib import Path

STATS_COLUMNS = (
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
)


def n50(vals):
    """Length-weighted median contig length."""
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


def new_header_for(assembler, first_token, assembly_sample_id, assembly_strategy, contig_count):
    if assembler == "megahit":
        m = re.search(r'(k\d+)_(\d+)', first_token)

        if m:
            return f"{assembly_sample_id}_{assembly_strategy}_{m.group(1)}_{m.group(2)}"

        return f"{assembly_sample_id}_{assembly_strategy}_k000_{contig_count}"

    if assembler == "metaspades":
        m = re.search(r'NODE_(\d+)', first_token)

        if m:
            return f"{assembly_sample_id}_{assembly_strategy}_NODE_{m.group(1)}"

        return f"{assembly_sample_id}_{assembly_strategy}_NODE_{contig_count}"

    raise RuntimeError(f"Unsupported assembler: {assembler}")


def rewrite_fasta(src_fasta, out_fasta, header_map, assembler, assembly_sample_id, assembly_strategy):
    """Stream src_fasta -> out_fasta, recording the header mapping.

    Returns (contig_count, list_of_contig_lengths).
    """
    lengths = []
    contig_count = 0
    current_len = 0

    with src_fasta.open() as inp, out_fasta.open("w") as out, header_map.open("w") as hmap:
        print("old_header", "new_header", sep="\t", file=hmap)

        for line in inp:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if contig_count > 0:
                    lengths.append(current_len)

                current_len = 0
                contig_count += 1

                old_header = line[1:].strip()
                first_token = old_header.split()[0] if old_header else f"contig_{contig_count}"

                new_header = new_header_for(
                    assembler,
                    first_token,
                    assembly_sample_id,
                    assembly_strategy,
                    contig_count,
                )

                print(f">{new_header}", file=out)
                print(old_header, new_header, sep="\t", file=hmap)

            else:
                seq = line.strip()
                current_len += len(seq)
                print(seq, file=out)

    if contig_count > 0:
        lengths.append(current_len)

    return contig_count, lengths


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Rename assembly contigs and write Module 2 stats/manifest records.",
    )

    parser.add_argument("--src-fasta", required=True, type=Path,
                        help="Raw assembler output FASTA.")
    parser.add_argument("--out-fasta", required=True, type=Path,
                        help="Renamed FASTA to write.")
    parser.add_argument("--header-map", required=True, type=Path,
                        help="TSV mapping old header -> new header.")
    parser.add_argument("--stats-file", required=True, type=Path,
                        help="Per-assembly stats TSV to write.")
    parser.add_argument("--manifest-record", required=True, type=Path,
                        help="Single-row manifest record TSV to write.")

    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--safe-id", required=True)
    parser.add_argument("--assembly-sample-id", required=True,
                        help="Sample ID embedded in contig headers. For rarefied "
                             "assemblies this already includes the rarefaction letter.")
    parser.add_argument("--assembler", required=True, choices=["megahit", "metaspades"])
    parser.add_argument("--mode", required=True,
                        help="Assembly mode, e.g. 'single' or 'rarefied'.")
    parser.add_argument("--rarefaction-label", default="",
                        help="Rarefaction letter, empty for single assemblies.")
    parser.add_argument("--assembly-strategy", required=True,
                        help="Strategy code: A/B for single, C/D for rarefied.")
    parser.add_argument("--assembly-status", default="ok")
    parser.add_argument("--assembly-warning", default="")
    parser.add_argument("--published-fasta", required=True,
                        help="Path the renamed FASTA will be published to. Recorded "
                             "in the stats and manifest rows for downstream modules.")

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if not args.src_fasta.exists():
        print(f"ERROR: source FASTA does not exist: {args.src_fasta}", file=sys.stderr)
        return 1

    contig_count, lengths = rewrite_fasta(
        args.src_fasta,
        args.out_fasta,
        args.header_map,
        args.assembler,
        args.assembly_sample_id,
        args.assembly_strategy,
    )

    total_bp = sum(lengths)
    max_contig = max(lengths) if lengths else 0
    n50_value = n50(lengths)

    row = (
        args.sample_id,
        args.safe_id,
        args.assembly_sample_id,
        args.assembler,
        args.mode,
        args.rarefaction_label,
        args.assembly_strategy,
        args.assembly_status,
        args.assembly_warning,
        contig_count,
        total_bp,
        max_contig,
        n50_value,
        args.published_fasta,
    )

    with args.stats_file.open("w") as stats:
        print(*STATS_COLUMNS, sep="\t", file=stats)
        print(*row, sep="\t", file=stats)

    # The manifest record drops the status/warning and size columns; Module 3
    # reads exactly these eight fields.
    with args.manifest_record.open("w") as manifest:
        print(
            args.sample_id,
            args.safe_id,
            args.assembly_sample_id,
            args.assembler,
            args.mode,
            args.rarefaction_label,
            args.assembly_strategy,
            args.published_fasta,
            sep="\t",
            file=manifest,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
