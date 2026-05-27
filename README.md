<img width="1400" height="2100" alt="SAMWISE_FULL-git" src="https://github.com/user-attachments/assets/45c8d15d-8008-4cf9-b226-cf00aee47ec7"/>

# Welcome to SAMWISE! The automated, end-to-end metagenomic read processing program


# module_0_readprocess.nf

`module_0_readprocess.nf` is a Nextflow DSL2 workflow for initial read preprocessing and validation. It checks read file names, detects paired-end or interleaved read layouts, normalizes `.fq` filenames to `.fastq`, validates FASTQ structure, and runs FastQC.

---

## What this workflow does

This module performs the following steps:

1. **Checks for FastQC**
   - Looks for `fastqc` in the current environment.
   - If missing, attempts to install FastQC using `mamba`.

2. **Scans a user-provided input directory**
   - Checks all files in the directory by default.

3. **Validates read naming conventions**
   - Detects paired-end reads using either:
     - `_R1` / `_R2`
     - `_1` / `_2`
   - Detects interleaved reads using:
     - `_interleaved`

4. **Normalizes `.fq` filenames**
   - Accepts `.fq` and `.fq.gz`.
   - Internally normalizes them to `.fastq` and `.fastq.gz` using symlinks.
   - Original input files are not modified.

5. **Checks read pairing**
   - Ensures every read 1 file has a matching read 2 file.
   - Ensures samples are not supplied as both paired-end and interleaved.

6. **Validates FASTQ structure**
   - Confirms FASTQ records are 4-line records.
   - Checks headers start with `@`.
   - Checks separator lines start with `+`.
   - Checks sequence and quality strings are the same length.

7. **Runs FastQC**
   - Runs FastQC on validated, normalized read files.

---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- Java, required by Nextflow
- `mamba`, if using automatic FastQC installation
- `fastqc`, if automatic installation is disabled

---

## Usage:

The workflow requires an input directory containing sequencing read files.

```bash
nextflow run module_0_readprocess.nf \
  --input_dir ./reads \
  --fastqc_threads 8 \
  --output_dir ./fastqc_out
