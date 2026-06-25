<img width="1406" height="240" alt="SAMWISE_FULL-git" src="https://github.com/user-attachments/assets/5f56bf43-fa98-4ce3-99b1-ea05313e85c5" />

# Welcome to SAMWISE!

SAMWISE is an automated, end-to-end metagenomic read processing program. Here is a quick conceptual rundown of what this software can enable you to do via Nextflow DSL2 workflows.

<img width="1406" height="1577" alt="SAMWISE_FULL-git" src="https://github.com/user-attachments/assets/dfa4ee38-55ab-4b02-9837-7abfc0cb1528" />

---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- `mamba` (or `conda`) for automatic package installation
- `Java`, required by Nextflow (auto installed if installing with `mamba` / `conda`

To start with SAMWISE, you will want to make sure that you have `mamba` (or `conda`) installed. We recommend mamba, and you can follow the instructions here: https://conda-forge.org/download/

Then, you need to install NextFlow - this can be done via `mamba` / `conda`: https://anaconda.org/channels/bioconda/packages/nextflow/overview

Now, you are ready to proceed with SAMWISE!

---

# Quick recommended usage for the impatient
Alright alright - you want to run SAMWISE quickly and do not want to read through the full docs. All good. Here is how I would run this as an sbatch script on a server. 

NOTE: Your reads MUST be in one of the naming formats (_R1, _R2, _1, _2, _interleaved) and must
have extensions (.fq or .fastq - gzipped or not gzipped is fine).

```bash

# Pre-process your reads
nextflow run module_0_readprocess.nf \
--working_dir ./samwise-main \
--input_dir ./reads_dir \
--threads 36

# Trim your reads
nextflow run module_1_readtrimming.nf \
--working_dir ./samwise-main \
--threads 36

# Assemble your reads
nextflow run module_2_readassembly.nf \
--working_dir ./samwise-main \
--threads 36 \
--memory_gb 0 \
--megahit \
--metaspades \
--rarefied_assembly TRUE \
--rarefaction_splits 2

# Bin your assemblies
nextflow run module_3_binning.nf \
--working_dir ./samwise-main \
--threads 36 \
--metabat2 \
--quickbin \
--maxbin2 \
--min_scaffold_length 2500

# Refine the MAGs
nextflow run module_4_binrefinement.nf \
--working_dir ./samwise-main \
--threads 36

# Run a subtractive assembly
nextflow run module_5_subassembly.nf \
--working_dir ./samwise-main \
--threads 36 \
--megahit \
--metaspades \
--secondpass_metabat2 true \
--secondpass_quickbin true \
--secondpass_maxbin2 true

# Run final MAG annotation:
nextflow run module_6_magannotate.nf \
--working_dir ./samwise-main \
--run_checkm2 true \
--run_gtdbtk true \
--run_eggnog true \
--threads 36
```

```
Now that you got what you wanted, let's do a deep dive on the flags and modules that SAMWISE has to offer!
```

# Step 0: module_0_readprocess.nf

`module_0_readprocess.nf` is a workflow for initial read preprocessing and validation. It checks read file names, detects paired-end or interleaved read layouts, validates FASTQ structure (optional), and runs FastQC.

## Recommended Usage:

```bash
nextflow run module_0_readprocess.nf \
--input_dir ./reads_dir \
--threads 6 \
--working_dir ./output_samwise

# use `--` for any additional flags as well
```

## Module 0 Arguments:
| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Working directory where output files will be written. If provided, Module 0 outputs are written to `<working_dir>/module_0_readprocess`. |
| `input_dir` | `null` | Directory where the metagenomic reads are stored. Please see the required filename formats. |
| `threads` | `null` | Total number of threads to use. If provided, this overrides `fastqc_threads`. |
| `skip_validate` | `false` | Optional flag to skip validation of read files. This speeds up the process by skipping checks that confirm the reads are valid FASTQ files. If your reads are large, compressed, and you are certain they are valid FASTQ files, we recommend using this flag. |
| `file_pattern` | `*.{fastq.gz,fq.gz,fastq,fq}` | Glob/text pattern used to detect input files inside `input_dir`. The default detects `fastq.gz`, `fq.gz`, `fastq`, and `fq` files. We recommend leaving this unchanged. |
| `fastqc_threads` | `2` | Number of threads to use specifically for FastQC if the global `threads` argument is not passed. |
| `fastqc_version` | `0.12.1` | FastQC version to install if a different version is desired. |
| `auto_install` | `true` | Controls whether SAMWISE installs required packages, such as FastQC. If set to `false`, FastQC must already be available in your environment. |
| `outdir` | `<results_dir>/module_0_readprocess` | Only used if `working_dir` is not specified. |


```
*IMPORTANT*
Your reads MUST be in one of the naming formats (_R1, _R2, _1, _2, _interleaved) and must
have extensions (.fq or .fastq - gzipped or not gzipped is fine). We recommend naming
your reads something easy to detect that is all a single identifier, in other words,
remove "_", "-", ".", etc. and simply have files be like SampleA_R1.fastq.gz.
SAMWISE will trim the ids up to the first underscore and use it for downstream outputs.
```

For example, this is a valid reads_dir structure:
```
reads_dir/
├── SampleA_R1.fastq.gz
├── SampleA_R2.fastq.gz
├── SampleB_1.fq
├── SampleB_2.fq
└── SampleC_interleaved.fastq
```

---

# Step 1: module_1_readtrimming.nf

`module_1_readtrimming.nf` is a workflow for trimming of reads that have been validated in module 0. Module 1 will trim all reads with `fastp`, run `fastqc` on the trimmed reads, and provide a summary table of the trimming statistics.

Module 1 automatically detects folders and file inputs form Module 0, as such, not many flags are needed if you are running the full
SAMWISE workflow.

## Recommended Usage:

```bash
nextflow run module_1_readtrimming.nf \
--working_dir ./output_samwise \
--threads 6

# use `--` for any additional flags as well
```

## Module 1 Arguments

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Main working/results directory for the pipeline. If provided, Module 1 outputs are written to `<working_dir>/module_1_readtrimming`. |
| `input_manifest` | `null` | Input read manifest file. This is typically the `read_manifest.tsv` produced by Module 0. |
| `output_dir` | `null` | Alternative output directory used only if `--working_dir` is not provided. |
| `threads` | `null` | Global thread override. If provided, this can be used instead of module-specific thread settings. |
| `fastp_threads` | `4` | Number of threads to use for `fastp` read trimming. |
| `fastqc_threads` | `2` | Number of threads to use for FastQC after trimming. |
| `fastp_version` | `0.23.4` | Version of `fastp` to install/use. |
| `fastqc_version` | `0.12.1` | Version of FastQC to install/use. |
| `auto_install` | `true` | Whether to automatically install required tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `tool_env_dir` | `null` | Optional custom path for the conda environment containing Module 1 tools. |
| `publish_trimmed_mode` | `symlink` | How trimmed read files are published to the output directory. Options are `symlink`, `copy`, or `move`. |
| `compression` | `4` | Compression level used by `fastp` for output FASTQ files. |
| `detect_adapter_for_pe` | `true` | Enables adapter sequence detection for paired-end reads in `fastp`. |
| `enable_correction` | `false` | Enables base correction for paired-end data in `fastp`. |
| `cut_front` | `true` | Enables quality trimming from the front/start of reads. |
| `cut_tail` | `true` | Enables quality trimming from the tail/end of reads. |
| `cut_window_size` | `4` | Window size used for sliding-window quality trimming. |
| `cut_mean_quality` | `30` | Mean quality threshold required within the trimming window. |
| `qualified_quality_phred` | `30` | Phred quality score threshold for a base to be considered qualified. |
| `unqualified_percent` | `40` | Maximum allowed percentage of unqualified bases in a read. |
| `n_base_limit` | `5` | Maximum number of `N` bases allowed in a read before filtering. |
| `length_required` | `75` | Minimum read length required after trimming/filtering. Reads shorter than this are discarded. |
| `trim_poly_g` | `false` | Enables poly-G tail trimming. Often useful for reads generated on two-color Illumina platforms. |
| `trim_poly_x` | `false` | Enables poly-X tail trimming. |
| `results_dir` | null | Internal results directory. Uses `--working_dir` if provided, otherwise `--output_dir`, otherwise `.`. Usually does not need to be set directly. |
| `module0_outdir` | `<results_dir>/module_0_readprocess` | Expected Module 0 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_1_readtrimming` | Module 1 output directory. Usually derived automatically and does not need to be set directly. |

---

```
*IMPORTANT*
Currently, rarefied assemblies are set to run as paralell processes to single assemblies to speed things up.
In theory, they should play nice. However, if you run into issues with clobbering memory, we will be working
on adding a flag so that the rarefied assemblies run only after single assemblies are complete.
```

# Step 2: module_2_readassembly.nf

`module_2_readassembly.nf` is a workflow for assembly of reads that have been trimmed in module 1. 
This module will run single assemblies using either `megahit`, `metaspades` or both, and then also 
do rarefied assemblies (if specified) using `megahit`. It will rename and standardize all assembled
output scaffolds.

## Recommended Usage:

```bash
## metaspades only run:
nextflow run module_2_readassembly.nf \
--working_dir ./output_samwise \
--threads 6 \
--memory_gb 0 \
--metaspades

## metaspades and megahit with rarefied assemblies:
nextflow run module_2_readassembly.nf \
--working_dir ./output_samwise \
--threads 6 \
--megahit \
--metaspades \
--memory_gb 0 \
--rarefied_assembly TRUE \
--rarefaction_splits 2

#if on a mac, megahit running on more than 1 thread doesnt work, so there is an explicit arg:
#--megahit_threads that you can set separately from the global --threads (which will set it for both)

# use `--` for any additional flags as well
```

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Main working/results directory for the pipeline. If provided, Module 2 outputs are written to `<working_dir>/module_2_readassembly`. |
| `input_manifest` | `null` | Input manifest file containing reads for assembly. This is typically produced by Module 1 after read trimming. |
| `output_dir` | `null` | Alternative output directory used only if `--working_dir` is not provided. |
| `megahit` | `false` | Enables assembly with MEGAHIT. |
| `metaspades` | `false` | Enables assembly with metaSPAdes. |
| `single_assembly` | `true` | Performs a single assembly using the available input reads. |
| `rarefied_assembly` | `false` | Enables rarefied assembly mode. |
| `rarefaction_splits` | `2` | Number of rarefaction splits to generate when `--rarefied_assembly` is enabled. |
| `megahit_version` | `1.2.9` | Version of MEGAHIT to install/use. |
| `spades_version` | `4.2.0` | Version of SPAdes/metaSPAdes to install/use. |
| `auto_install` | `true` | Whether to automatically install required assembly tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `tool_env_dir` | `null` | Optional custom path for the conda environment containing Module 2 assembly tools. |
| `threads` | `null` | Global thread override. If provided, this can be used instead of module-specific thread settings. |
| `assembly_threads` | `4` | Number of threads to use for assembly if `--threads` is not provided. |
| `memory_gb` | `0` | Global memory limit in GB for assembly processes. Use `0` to leave memory unset / at max. |
| `megahit_threads` | `null` | Optional MEGAHIT-specific thread override. If provided, this overrides the general assembly thread setting for MEGAHIT. This is really only important for mac users that need to specify a single thread for it to work. |
| `megahit_preset` | `meta-large` | MEGAHIT preset to use for assembly. Default is `meta-large`. |
| `publish_assemblies_mode` | `symlink` | How final assembly files are published to the output directory. Common options are `symlink`, `copy`, or `move`. |
| `results_dir` | null | Internal results directory. Uses `--working_dir` if provided, otherwise `--output_dir`, otherwise `.`. Usually does not need to be set directly. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_2_readassembly` | Module 2 output directory. Usually derived automatically and does not need to be set directly. |

```
*IMPORTANT*
The default assembly outputs get written into the NextFlow work directories to save space. If you want it to
write out the output assemblies into a more accessible location, you can set publish_assemblies_mode to be
`copy`. Argument `move` here would also work but may cause issues with NextFlow not finding what it needs.
```

# Step 2b (optional): module_2b_coassembly.nf

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
#Run with max mem on 5 threads

nextflow run module_2b_coassembly.nf \
--working_dir ./output_samwise-main \
--coassembly_groups ./coassembly_manifest.txt \
--threads 6 \
--memory_gb 0

# use `--` for any additional flags as well
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
   - Automatically scans for all possible assemblies from both module 2 and 2b

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

# use `--` for any additional flags as well
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

# use `--` for any additional flags as well
```

# Step 5 (optional): module_5_subassembly.nf

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

# use `--` for any additional flags as well\
```
# Step 6: module_6_magannotate.nf

This module performs the following steps:

1. **Checks for dependencies and installs them if necessary**
   -Module 2 will attempt to find the MAG annotation tools and will annotate them if not found. 
   -Module 2 will download required databases for each tool if needed, but arguments can be passed to directly point to dbs.
2. **Runs CheckM2, GTDB-tk, and Eggnog (or DRAM2)**
   -We note that for right now, DRAM2 has been replaced with eggnog v2 since DRAM2 is undergoing significant development and is currently not fully installable. Once development is finished, we will update our module to incldue both DRAM2 and eggnog (and can updated to eggnog v3 once available as well).

## Usage:
```bash
nextflow run module_6_magannotate.nf \
--working_dir ./output_samwise \
--run_checkm2 true \
--run_gtdbtk true \
--run_eggnog true \
--threads 32

# If you already pre-downloaded the databases and have them installed elsewhere, you can directly pass arguments:
# --checkm2_db_path /path/to/uniref100.KO.1.dmnd
# --gtdbtk_data_path /path/to/gtdbtk/database_directory
# --eggnog_data_path /path/to/eggnog/database_directory

# You can also explicitly pass which directories you want the files downloaded into, for example:

# --checkm2_db_dir /path/to/download/checkm2_db_dir
# --gtdbtk_db_dir /path/to/download/gtdbtk_db_dir
# --eggnog_data_path /path/to/eggnog/database_directory

# use `--` for any additional flags as well
```


### Advanced workflow explanations:

## Module 0 Steps:
1. **Check for FastQC**

   The module first checks whether `fastqc` is available in the current environment.

   - If FastQC is found, the existing installation is used.
   - If FastQC is not found and `--auto_install true` is set, the module attempts to install FastQC using `mamba` or `conda`.
   - If `--auto_install false` is set, FastQC must already be available in your environment.

2. **Check Input Read Names**

   The module scans `--input_dir` for FASTQ files and checks that filenames follow supported naming conventions.

   Supported paired-end read patterns include:

   - `sample_R1.fastq.gz` and `sample_R2.fastq.gz`
   - `sample_R1.fastq` and `sample_R2.fastq`
   - `sample_1.fastq.gz` and `sample_2.fastq.gz`
   - `sample_1.fastq` and `sample_2.fastq`

   Supported interleaved read pattern:

   - `sample_interleaved.fastq.gz`
   - `sample_interleaved.fastq`

   The module also checks for common problems such as:

   - Missing R1 or R2 files
   - Duplicate read files for the same sample
   - Mixed naming styles, such as using both `_R1/_R2` and `_1/_2`
   - Samples with both paired-end and interleaved reads
   - Unsupported FASTQ filenames

   Original input files are not modified.

3. **Validate FASTQ Structure**

   Unless `--skip_validate` is used, the module validates the structure of each FASTQ file.

   The validation step checks that:

   - FASTQ records contain 4 lines
   - Header lines start with `@`
   - Separator lines start with `+`
   - Sequence and quality strings are the same length
   - The file is not empty
   - The total number of lines is divisible by 4

   This step can be slow for large compressed FASTQ files. If you are confident your reads are valid FASTQ files, you can skip this step using:
   `--skip_validate`


## Module 1 Steps:

1. **Trims reads using fastp**
   - Looks for `fastp` in the current environment.
   - If missing, attempts to install FastP using `mamba`.
   - Trims reads using fastp default trimming parameters - highly customizable with any and all fastp flags if needed.

2. **Provides trimming statistics**
   - Pre / Post trimming read quality statitsitcs in tabulated format

3. **Re-runs fastqc on trimmed reads**
   - Re-analysis of the fastqc outputs to confirm succesful trimming.


## Module 2 Steps:

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
