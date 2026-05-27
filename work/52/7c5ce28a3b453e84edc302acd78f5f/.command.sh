#!/bin/bash -ue
set -euo pipefail

    mkdir -p valid_reads

    cat > staged_files.list <<'EOF'
SM_0533_SRMG_trimmed_interleaved.fastq.gz
module_0_readprocess.nf
DN-Isolate20-Chitin-Soil_R2.fastq
DN-Isolate20-Chitin-Soil_R1.fastq
EOF

    python3 - <<'PY'
import os
import re
import sys
from pathlib import Path

report_path = Path("read_naming_report.txt")
manifest_path = Path("read_manifest.tsv")
valid_reads_dir = Path("valid_reads")
valid_reads_dir.mkdir(exist_ok=True)

r1_files = {}
r2_files = {}
interleaved_files = {}

r1_style = {}
r2_style = {}

seen_keys = set()
errors = 0
warnings = 0

read1_re = re.compile(r'^(.+)_(R1|1)\.(fastq|fq)(\.gz)?)
read2_re = re.compile(r'^(.+)_(R2|2)\.(fastq|fq)(\.gz)?)
interleaved_re = re.compile(r'^(.+)_interleaved\.(fastq|fq)(\.gz)?)
fastq_like_re = re.compile(r'.*\.(fastq|fq)(\.gz)?)


def normalize_read_name(filename):
    if filename.endswith(".fastq.gz"):
        return filename
    if filename.endswith(".fq.gz"):
        return filename[:-6] + ".fastq.gz"
    if filename.endswith(".fastq"):
        return filename
    if filename.endswith(".fq"):
        return filename[:-3] + ".fastq"
    return filename


def safe_symlink(src, dest):
    src_real = os.path.realpath(src)
    dest = Path(dest)

    if dest.exists() or dest.is_symlink():
        dest.unlink()

    os.symlink(src_real, dest)


with report_path.open("w") as report, manifest_path.open("w") as manifest:

    def log(message=""):
        print(message, file=report)

    print("sample_id", "layout", "read1", "read2", "interleaved", sep="\t", file=manifest)

    log("Read naming and pairing report")
    log("Started")
    log("----------------------------------------")
    log("")

    with open("staged_files.list") as handle:
        staged_files = [line.strip() for line in handle if line.strip()]

    if not staged_files:
        log("ERROR: No files were found in the input directory.")
        errors += 1
