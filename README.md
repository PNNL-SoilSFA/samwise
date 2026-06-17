<img width="1406" height="240" alt="SAMWISE_FULL-git" src="https://github.com/user-attachments/assets/5f56bf43-fa98-4ce3-99b1-ea05313e85c5" />

# Welcome to SAMWISE!

SAMWISE is an automated, end-to-end metagenomic read processing program. Here is a quick conceptual rundown of what this software can enable you to do via Nextflow DSL2 workflows.

<img width="1406" height="1577" alt="SAMWISE_FULL-git" src="https://github.com/user-attachments/assets/469fe2de-8809-40af-8002-152213046439" />

---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- Java, required by Nextflow
- `mamba` (or `conda`), if using automatic package installation
  
---

# Step 0: module_0_readprocess.nf

`module_0_readprocess.nf` is a workflow for initial read preprocessing and validation. It checks read file names, detects paired-end or interleaved read layouts, validates FASTQ structure (optional), and runs FastQC.

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
   - Internally normalizes them to `.fastq` and `.fastq.gz` using symlinks for downstream.
   - Original input files are not modified.

```
*IMPORTANT*
Your reads MUST be in one of the naming formats shown in 4 and 5.
They are allowed to be gzipped or unzipped.
```

For example, this is a valid dir structure:
```
reads/
├── SampleA_R1.fastq.gz
├── SampleA_R2.fastq.gz
├── SampleB_1.fq
├── SampleB_2.fq
└── SampleC_interleaved.fastq
```
5. **Checks read pairing**
   - Ensures every read 1 file has a matching read 2 file.
   - Ensures samples are not supplied as both paired-end and interleaved.

6. **Validates FASTQ structure**
   - Confirms FASTQ records are 4-line records.
   - Checks headers start with `@`.
   - Checks separator lines start with `+`.
   - Checks sequence and quality strings are the same length.

```
*Note that this process can be very slow if files are gzipped.
*If you are certain that your file is a fastq file, you can safely skip this step with the flag --skip_validate
```

7. **Runs FastQC**
   - Runs FastQC on validated, normalized read files.

## Usage:

The workflow requires an input directory containing sequencing read files.

```bash
nextflow run module_0_readprocess.nf \
--input_dir ./reads_dir \
--threads 6 \
--working_dir ./output_samwise
```

# Step 1: module_1_readtrimming.nf

`module_1_readtrimming.nf` is a workflow for trimming of reads that have been validated in module 0.

This module performs the following steps:

1. **Trims reads using fastp**
   - Looks for `fastp` in the current environment.
   - If missing, attempts to install FastP using `mamba`.
   - Trims reads using fastp default trimming parameters - highly customizable with any and all fastp flags if needed.

2. **Provides trimming statistics**
   - Pre / Post trimming read quality statitsitcs in tabulated format

3. **Re-runs fastqc on trimmed reads**
   - Re-analysis of the fastqc outputs to confirm succesful trimming.

This workflow takes in the same working directory that was generated above and finds whatever files it needs.

## Usage:
```bash
nextflow run module_1_readtrimming.nf \
--working_dir ./output_samwise \
--threads 6
```

# Step 2: module_2_readassembly.nf

`module_2_readassembly.nf` is a workflow for assembly of reads that have been trimmed in module 1.

This module performs the following steps:

1. **Checks for assembly software and installs if necessary**
   - Looks for MEGAHIT and metaSPAdes and installs into a local conda environment if needed.
   
2. **Assembles using multiple assemblers and assembly methods**
   - Users can choose either assembler or both: with flags `--megahit` and `--metaspades`
   - This step can also perform rarefied assemblies by adding the flag `--rarefied_assembly TRUE` and specifying how many "fragments" you want the reads to be split into with `--rarefaction_splits #` (default is 2). Tl;dr - this will split the fastq files into # of split files and assemble them individually via a round robin by pair index approach.  
     
3. **Renames scaffold outputs and provides assembly statistics**
   - Implements naming scheme specifically:

```bash
   MEGAHIT single assembly:
   SampleID_A_k###_#
    
metaSPAdes single assembly:
   SampleID_B_NODE_#

MEGAHIT rarefied assembly:
   SampleIDa_C_k###_#
   SampleIDb_C_k###_#
   SampleIDc_C_k###_#

metaSPAdes rarefied assembly:
   SampleIDa_D_NODE_#
   SampleIDb_D_NODE_#
   SampleIDc_D_NODE_#
```   

This workflow takes in the same working directory that was generated above and finds whatever files it needs.

## Usage:
```bash
## metaspades only run:

nextflow run module_2_readassembly.nf \
--working_dir ./output_samwise \
--threads 6 \
--memory_gb 0 \
--metaspades

#--memory_gb 0 specifies use 90% of available memory

## metaspades and megahit with rarefied assemblies:
nextflow run module_2_readassembly.nf \
--working_dir ./output_samwise \
--threads 5 \
--megahit \
--metaspades \
--megahit_threads 1 \
--memory_gb 0 \
--rarefied_assembly TRUE \
--rarefaction_splits 2

#if on a mac, megahit running on more than 1 thread doesnt play nice, so there is an explicit --megahit_threads you can set separately from the global argument --threads which will set it for both.
```

# Step 2b: module_2b_coassembly.nf

`module_2b_coassembly.nf` performs **grouped co-assembly** from Module 1 trimmed reads using **MEGAHIT only**.

This module is designed to run alongside the normal Module 2 assembly workflow. It produces Module-3-compatible manifests so that Module 3 can automatically bin co-assemblies using the exact concatenated reads that were used to generate each co-assembly.

1. **Checks for assembly software and installs if necessary**
   - Looks for MEGAHIT and installs if needed
   
2. **Co-assembles reads as specified within reads manifest**
   - Users specify which read groupings are relevant for co-assembly and software automatically reads in from the previous Module 0 and Module 1.
   - SAMWISE concatenates the reads for each group into one interleaved FASTQ file (used later for binning).
   - Runs MEGAHIT co-assembly on each grouped interleaved FASTQ and renames scaffolds with letter E.
   
## Usage:
```bash
## 

nextflow run module_2b_coassembly.nf \
--working_dir ./output_samwise-main \
--coassembly_groups ./coassembly_manifest.txt \
--threads 5 \
--memory_gb 0

```

Coassembly_manifest.txt must be a tab-separated table and contain two columns:

| read_or_sample_id | group_id |
| sampleA	| group_1 |
| sampleB	| group_1 |
| sampleC	| group_2 |
| sampleD	| group_2 |
  
# Step 3: module_3_binning.nf

This module performs the following steps:

1. **Checks for binning software and installs if necessary**
   - Looks for Quickbin, Metabat2, and MaxBIN2 and installs if needed
   
2. **Assembles using multiple assemblers and assembly methods**
   - Users can choose either assembler or both: with flags `--quickbin`, `--metabat2`, `--maxbin2`
   - Binning will be run on all assemblies generated from the prior modules
   - Minimum scaffold length required for binning can be modified with `--min_scaffold_length` flag, default is 2500.

## Usage:
```bash
## All binners run on scaffolds >2500bp unless otherwise specified:

nextflow run module_3_binning.nf \
--working_dir ./output_samwise \
--threads 5 \
--metabat2 \
--quickbin \
--maxbin2 \
--min_scaffold_length 2500

# --publish_bins_mode copy tells code to copy the genomes instead of making a symlink

```

# Step 4: module_4_binrefinement.nf

This module performs the following steps:

1. **Checks for MAGScoT dependencies and installs if necessary**
   - Looks for r, r-base, r-optparse, r-dplyr, r-readr, r-funr, r-digest, hmmer, prodigal, parallel, pandas and installs if needed
   
2. **Runs the MAGScoT workflow and generates refined MAGs**
   - Runs MAGScoT which uses GTDB r207 to re-group / rebin genomes and outputs cleaned, refined MAG set.

## Usage:
```bash
nextflow run module_4_binrefinement.nf \
--working_dir ./output_samwise \
--magscot_threshold 0

#MAGScoT original code sets this threshold at 0.5, but since we are doing gtdb + checkm runs after on the latest databases, its better to just pass this as default 0 and retain all possible MAGs. Feel free to change that --magscot_threshold param to 0.5
```

# Step 5 (OPTIONAL): module_5_subassembly.nf

This module performs the following steps:

1. **Checks for dependencies and installs them if necessary**
2. **Runs subtractive assembly workflow**
3. **Runs binning**
   - Same as module 3: Runs the binning pipeline with all 3 binners (if needed).
4. **Runs consolidated MAG refinement**
   - Same as Module 4: Runs MAGScoT which uses GTDB r207 to re-group / rebin genomes and outputs cleaned, refined MAG set.

## Usage:
```bash
nextflow run module_5_subassembly.nf \
--working_dir ./output_samwise \
--threads 20 \
--megahit \
--metaspades \
--secondpass_metabat2 true \
--secondpass_quickbin true \
--secondpass_maxbin2 true \
--run_second_pass_binning_refinement true

# --run_second_pass_binning_refinement specifies whether or not you want it to re-bin after subassembly - some users may want to disable this if they want to make sure subassemblies are worth performing after looking at the assembly stats, but most should leave on. Default is true.
```
