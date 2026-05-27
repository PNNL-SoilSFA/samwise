#!/bin/bash -ue
set -euo pipefail

STATUS_FILE="fastqc_install_status.txt"

if command -v fastqc >/dev/null 2>&1; then
    echo "FastQC already installed: $(command -v fastqc)" > "$STATUS_FILE"
    fastqc --version >> "$STATUS_FILE" 2>&1 || true
    exit 0
fi

echo "FastQC not found in PATH." > "$STATUS_FILE"

if [[ "true" != "true" ]]; then
    echo "Auto-install disabled. Please install FastQC manually or run with --auto_install true" >> "$STATUS_FILE"
    exit 1
fi

if ! command -v mamba >/dev/null 2>&1; then
    echo "mamba not found. Please install mamba first." >> "$STATUS_FILE"
    exit 1
fi

echo "Attempting to install FastQC using mamba..." >> "$STATUS_FILE"

if mamba install -y -c bioconda -c conda-forge fastqc >> "$STATUS_FILE" 2>&1; then
    echo "FastQC installation successful." >> "$STATUS_FILE"
else
    echo "FastQC installation failed." >> "$STATUS_FILE"
    exit 1
fi

if command -v fastqc >/dev/null 2>&1; then
    echo "FastQC path: $(command -v fastqc)" >> "$STATUS_FILE"
    fastqc --version >> "$STATUS_FILE" 2>&1 || true
else
    echo "FastQC still not detected after installation." >> "$STATUS_FILE"
    exit 1
fi
