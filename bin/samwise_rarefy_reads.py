#!/usr/bin/env python3
"""Write one rarefied (strided) subset of a trimmed read set.

Standalone port of the Python that used to be embedded as a heredoc in the
ASSEMBLE_RAREFIED process of module_2_readassembly.nf.

Subsetting is deterministic striding over read *pairs*, not random sampling:
pair i goes to subset (i % split_count). So --split-count 2 produces two
disjoint, interleaved halves of the library.

Run standalone:

  samwise_rarefy_reads.py \
      --layout paired --read1 R1.fq.gz --read2 R2.fq.gz \
      --split-index 0 --split-count 2 \
      --out-r1 subset_a_R1.fastq.gz --out-r2 subset_a_R2.fastq.gz

Note: --split-index is 0-based, so subset 'a' of 2 is --split-index 0.

This currently writes a single subset per invocation, matching how the
Nextflow process calls it. Producing every subset in one pass over the reads
would avoid re-decompressing the input once per (assembler x split); the
write_*_subset functions below already stream, so that change is a loop over
output handles rather than a rewrite.
"""

import argparse
import gzip
import sys
from pathlib import Path
from itertools import zip_longest

def open_fastq(path):
    path = str(path)

    if path.endswith(".gz"):
        return gzip.open(path, "rt")

    return open(path, "rt")


def read_fastq_records(path):
    """Yield 4-line FASTQ records as (header, seq, plus, qual) tuples."""
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


def write_paired_subset(read1, read2, out_r1, out_r2, split_idx, split_count):
    if not Path(read1).exists():
        raise RuntimeError(f"Read 1 file does not exist: {read1}")

    if not Path(read2).exists():
        raise RuntimeError(f"Read 2 file does not exist: {read2}")

    r1_iter = read_fastq_records(read1)
    r2_iter = read_fastq_records(read2)

    records_written = 0

    with gzip.open(out_r1, "wt") as o1, gzip.open(out_r2, "wt") as o2:
        for idx, (rec1, rec2) in enumerate(zip_longest(r1_iter, r2_iter)):
            if rec1 is None or rec2 is None:
                raise RuntimeError(
                    "Paired FASTQ files have unequal record counts: "
                    f"{read1} and {read2}"
                )

            if idx % split_count == split_idx:
                o1.writelines(rec1)
                o2.writelines(rec2)
                records_written += 1

    return records_written


def write_interleaved_subset(interleaved, out_12, split_idx, split_count):
    if not Path(interleaved).exists():
        raise RuntimeError(f"Interleaved file does not exist: {interleaved}")

    rec_iter = read_fastq_records(interleaved)

    records_written = 0

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
                raise RuntimeError(
                    f"Interleaved FASTQ has odd number of records: {interleaved}"
                )

            if pair_idx % split_count == split_idx:
                out.writelines(rec1)
                out.writelines(rec2)
                records_written += 1

            pair_idx += 1

    return records_written


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Write one strided rarefied subset of a trimmed read set.",
    )

    parser.add_argument("--layout", required=True, choices=["paired", "interleaved"])
    parser.add_argument("--read1", default="", help="R1 FASTQ (paired layout).")
    parser.add_argument("--read2", default="", help="R2 FASTQ (paired layout).")
    parser.add_argument("--interleaved", default="", help="Interleaved FASTQ.")

    parser.add_argument("--split-index", required=True, type=int,
                        help="0-based index of the subset to emit.")
    parser.add_argument("--split-count", required=True, type=int,
                        help="Total number of subsets the library is divided into.")

    parser.add_argument("--out-r1", default="", help="Output R1 (paired layout).")
    parser.add_argument("--out-r2", default="", help="Output R2 (paired layout).")
    parser.add_argument("--out-interleaved", default="",
                        help="Output interleaved FASTQ (interleaved layout).")

    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.split_count < 2:
        print(f"ERROR: --split-count must be >= 2, got {args.split_count}", file=sys.stderr)
        return 1

    if not 0 <= args.split_index < args.split_count:
        print(
            f"ERROR: --split-index {args.split_index} out of range for "
            f"--split-count {args.split_count}",
            file=sys.stderr,
        )
        return 1

    if args.layout == "paired":
        if not args.out_r1 or not args.out_r2:
            print("ERROR: paired layout requires --out-r1 and --out-r2", file=sys.stderr)
            return 1

        records_written = write_paired_subset(
            args.read1, args.read2, args.out_r1, args.out_r2,
            args.split_index, args.split_count,
        )

    else:
        if not args.out_interleaved:
            print("ERROR: interleaved layout requires --out-interleaved", file=sys.stderr)
            return 1

        records_written = write_interleaved_subset(
            args.interleaved, args.out_interleaved,
            args.split_index, args.split_count,
        )

    if records_written == 0:
        raise RuntimeError(
            f"Rarefied subset {args.split_index + 1} of {args.split_count} "
            f"contains zero pairs/fragments."
        )

    print(
        f"Wrote {records_written} pairs to rarefied subset "
        f"{args.split_index + 1} of {args.split_count}"
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())