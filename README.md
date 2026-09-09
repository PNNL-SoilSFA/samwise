![SAMWISE title](images/SAMWISE_title_v2.png)

# Welcome to SAMWISE!

SAMWISE is an automated, end-to-end metagenomic read processing program. Here is a quick conceptual rundown of what this software can enable you to do via Nextflow DSL2 workflows.

![SAMWISE workflow](images/SAMWISE_FULL-manuscript_v2.png)


---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- `mamba` (or `conda`) for automatic package installation
- `Java`, required by Nextflow (auto installed if installing with `mamba` / `conda`

To start with SAMWISE, you will want to make sure that you have `mamba` (or `conda`) installed. We recommend mamba, and you can follow the instructions here: https://conda-forge.org/download/

Then, you need to install NextFlow - this can be done via `mamba` / `conda`: https://anaconda.org/channels/bioconda/packages/nextflow/overview
We recommend that you install NextFlow into its own, separate environment from your base environment. For example, with `mamba install -n nextflow -c bioconda nextflow` Then, when running SAMWISE, make sure that you activate your NextFlow environment with `mamba activate nextflow`!

Once NextFlow is installed, go ahead and clone this repo or download it / extract. You can click on `clone repo` in the top right on GitHub or just download the whole thing. Then, change directory into the directory of the cloned repo: `cd ./samwise-main`

A helpful note: the cloned `samwise` directory holds the Nextflow workflows, helper scripts, and bundled dependencies that are needed. If you would like, SAMWISE can keep source assets and results separate with two directory parameters:
- `--samwise_dir`: SAMWISE installation/clone directory containing the `.nf` workflows, `bin/`, and `dependencies/`.
- `--working_dir`: results directory. Module outputs, generated Conda environments, and downloaded tool databases are written below this directory.

For simplicity, if either `--working_dir` or `--samwise_dir` are omitted, they inherit the same value, meaning that existing in-place commands continue to work and all outputs are written to whatever directory is provided. If BOTH `--working_dir` and `--samwise_dir` are supplied on every module invocation, it will then split up the code directory and output directory.

Now, you are ready to proceed with SAMWISE!

---
![SAMWISE quickstart](images/quick_start.png)
---

Alright alright - you want to run SAMWISE quickly and do not want to read through the full docs. Here is how I would run this as an sbatch script on a server.

NOTE: Your reads MUST be in one of the naming formats (_R1, _R2, _1, _2, _interleaved) and must
have extensions (.fq or .fastq - gzipped or not gzipped is fine). See module_0 info below! The --input_dir flag just needs to point to any dir that has reads

```bash

#A quick note: Feel free to remove -c, --slurm_account, and --maxparallel_fastqc
#if you do not have a slurm manager or allocation. These are optional!

# Pre-process your reads
nextflow run module_0_readprocess.nf \
-c ./bin/module_0_slurm.config \
--working_dir ./samwise-main \
--slurm_account ChargeAccountID (or User Account ID) \
--max_parallel_fastqc 10 \
--working_dir ./samwise-main \
--input_dir ./reads_dir \
--threads 36

# Trim your reads
nextflow run module_1_readtrimming.nf \
-c ./bin/module_1_slurm.config \
--working_dir ./samwise-main \
--slurm_account ChargeAccountID (or User Account ID) \
--max_parallel_trimming 10 \
--working_dir ./samwise-main \
--threads 36

# Assemble reads in parallel across SLURM nodes
nextflow run module_2_readassembly.nf \
-c ./bin/module_2_slurm.config \
--working_dir ./samwise-main \
--slurm_account ChargeAccountID (or User Account ID) \
--max_parallel_assemblies 8 \
--threads 36 \
--memory_gb 0 \
--megahit \
--metaspades \
--rarefied_assembly TRUE \
--rarefaction_splits 2

# Keep --memory_gb 0 with this SLURM config: each assembly receives a whole
# node and auto-detects its available memory. Adjust --max_parallel_assemblies
# to the number of simultaneous assembly jobs permitted by your allocation.

# Coassemble each group in parallel across SLURM nodes
nextflow run module_2b_coassembly.nf \
-c ./bin/module_2b_slurm.config \
--working_dir ./samwise-main \
--coassembly_groups ./coassembly_manifest.txt \
--slurm_account ChargeAccountID (or User Account ID) \
--threads 36 \
--memory_gb 0 \
--max_parallel_coassemblies 4

# Each distinct group becomes one MEGAHIT job. Keep --memory_gb 0 with this
# SLURM profile; set --max_parallel_coassemblies to the allowed group-job count.

#check the module notes for coassembly_manifest.txt examples / format

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

#Note: you can only run this if you have run more than 1 binner.

# Run a subtractive assembly
nextflow run module_5_subassembly.nf \
-c ./bin/module_5_slurm.config \
--working_dir ./samwise-main \
--coassembly_groups ./coassembly_manifest.txt \
--slurm_account ChargeAccountID (or User Account ID) \
--max_parallel_subassembly 8 \
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
--run_microtrait true \
--threads 36

#Note: If running only 1 binner and not using MAGScoT, pass the MAG manifest from Module 3 directly to Module 6 with --input_mag_manifest.

# To save compute time, you can pre-download the databases before you run module 6.
# For this, see: CheckM2: https://zenodo.org/records/14897628, gtdbtk: https://ecogenomics.github.io/GTDBTk/installing/index.html,
# eggnog: https://github.com/eggnogdb/eggnog-mapper; command: download_eggnog_data.py --data_dir /path/to/eggnog-data
# If you already pre-downloaded the gtdb, checkm2, and eggnog databases, 
# you can directly pass the paths as arguments:
# --checkm2_db_path /path/to/uniref100.KO.1.dmnd
# --gtdbtk_data_path /path/to/gtdbtk/database_directory
# --eggnog_data_path /path/to/eggnog/database_directory

nextflow run module_7_gems.nf \
--working_dir ./output_samwise

nextflow run AuxModule_1_assemblyAnnotate.nf \
--working_dir ./output_samwise \
--threads 36 \
--min_scaffold_bp 1000

nextflow run AuxModule_2_mvp.nf \
--working_dir ./output_samwise \
--mvp_modules "0,1,2,3,6" \
--threads 36 \
--memory_gb 0
```

```
Now that you got what you wanted, let's do a deep dive on the flags and modules that SAMWISE has to offer! 
First, a quick note. If you ever have a module (for example, an assembly module) halt because of time or whatever issue, you can resume the assembly by simply passing "-resume" as an argument for that module.
```

![SAMWISE step0](images/step_0.png)

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
| `working_dir` | `null` | Results root. Module 0 outputs are written to `<working_dir>/module_0_readprocess`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_dir` | `null` | Directory where the metagenomic reads are stored. Please see the required filename formats. |
| `threads` | `null` | Total number of threads to use. If provided, this overrides `fastqc_threads`. |
| `skip_validate` | `false` | Skip FASTQ structure validation. Naming and pairing checks still run; every named read is passed to FastQC. |
| `ignore_invalid_fastq` | `false` | Continue after FASTQ structure validation failures. Each invalid read is excluded from the Module 0 manifest, its validation report includes `READ IGNORED - BAD FASTQ FILE`, and `validation_reports/ignored_bad_fastq_files.tsv` records the exclusion. The default (`false`) stops the workflow on the first invalid FASTQ. Has no effect with `--skip_validate true`. |
| `file_pattern` | `*.{fastq.gz,fq.gz,fastq,fq}` | Glob/text pattern used to detect input files inside `input_dir`. The default detects `fastq.gz`, `fq.gz`, `fastq`, and `fq` files. We recommend leaving this unchanged. |
| `fastqc_threads` | `2` | Number of threads to use specifically for FastQC if the global `threads` argument is not passed. |
| `fastqc_version` | `0.12.1` | FastQC version to install if a different version is desired. |
| `auto_install` | `true` | Controls whether SAMWISE installs required packages, such as FastQC. If set to `false`, FastQC must already be available in your environment. |
| `outdir` | `<working_dir>/module_0_readprocess` | Derived Module 0 output directory. |

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

![SAMWISE step1](images/step_1.png)

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
| `working_dir` | `null` | Results root. Module 1 outputs are written to `<working_dir>/module_1_readtrimming`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_manifest` | `null` | Input read manifest file. This is typically the `read_manifest.tsv` produced by Module 0. |
| `threads` | `null` | Global thread override. If provided, this can be used instead of module-specific thread settings. |
| `fastp_threads` | `4` | Number of threads to use for `fastp` read trimming. |
| `fastqc_threads` | `2` | Number of threads to use for FastQC after trimming. |
| `fastp_version` | `0.23.4` | Version of `fastp` to install/use. |
| `fastqc_version` | `0.12.1` | Version of FastQC to install/use. |
| `auto_install` | `true` | Whether to automatically install required tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `tool_env_dir` | `null` | Optional custom path for the conda environment containing Module 1 tools. |
| `publish_trimmed_mode` | `symlink` | How trimmed read files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
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
| `results_dir` | null | Internal results root, derived from the resolved `working_dir`. Usually does not need to be set directly. |
| `module0_outdir` | `<results_dir>/module_0_readprocess` | Expected Module 0 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_1_readtrimming` | Module 1 output directory. Usually derived automatically and does not need to be set directly. |


![SAMWISE step2](images/step_2.png)

`module_2_readassembly.nf` assembles reads trimmed by Module 1. It can run single assemblies with `megahit`, `metaspades`, or both; when `--rarefied_assembly` is enabled, it can also run rarefied assemblies with either selected assembler. Final contigs are renamed and standardized for downstream Module 3 binning.

Assembly strategies in `summary/assembly_manifest.tsv` are `A` (single MEGAHIT), `B` (single metaSPAdes), `C` (rarefied MEGAHIT), and `D` (rarefied metaSPAdes).

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

# For parallel SLURM execution, add these arguments immediately after the
# workflow filename. Every single and rarefied assembly becomes an sbatch job.
# Keep --memory_gb 0 because the profile requests the full node memory.
#
# -c ./bin/module_2_slurm.config \
# --slurm_account ChargeAccountID \
# --max_parallel_assemblies 8

# use `--` for any additional flags as well
```

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 2 outputs are written to `<working_dir>/module_2_readassembly`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_manifest` | `null` | Input manifest file containing reads for assembly. This is typically produced by Module 1 after read trimming. |
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
| `memory_gb` | `0` | Global assembly memory limit in GB. Use `0` to auto-detect memory available to the task and reserve 10% headroom. |
| `megahit_threads` | `null` | Optional MEGAHIT-specific thread override. If provided, this overrides the general assembly thread setting for MEGAHIT. This is really only important for mac users that need to specify a single thread for it to work. |
| `megahit_preset` | `meta-large` | MEGAHIT preset to use for assembly. Default is `meta-large`. Could also use `meta-sensitive` |
| `publish_assemblies_mode` | `symlink` | How final assembly files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results root, derived from the resolved `working_dir`. Usually does not need to be set directly. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_2_readassembly` | Module 2 output directory. Usually derived automatically and does not need to be set directly. |

```
*IMPORTANT*
The default assembly outputs get written into the NextFlow work directories to save space. If you want it to
write out the output assemblies into a more accessible location, you can set publish_assemblies_mode to be
`copy`. Argument `move` here would also work but may cause issues with NextFlow not finding what it needs.
```

![SAMWISE step2b](images/step_2b.png)

`module_2b_coassembly.nf` performs **grouped co-assembly** from Module 1 trimmed reads using **MEGAHIT only**. This module is designed to run alongside the normal Module 2 assembly workflow. It produces Module-3-compatible manifests so that Module 3 can automatically bin co-assemblies using the exact concatenated reads that were used to generate each co-assembly. Coassembly contigs use assembly strategy `G`.

This module needs a co-assembly reads manifest that shows which reads you want to co-assemble.

Coassembly_manifest.txt must be a tab-separated table and contain two columns:

| read_or_sample_id | group_id |
|---|---|
| sampleA	| group_1 |
| sampleB	| group_1 |
| sampleC	| group_2 |
| sampleD	| group_2 |

## Usage:
```bash 
nextflow run module_2b_coassembly.nf \
--working_dir ./output_samwise-main \
--coassembly_groups ./coassembly_manifest.txt \
--threads 6 \
--memory_gb 0

# use `--` for any additional flags as well
```

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 2B writes to `<working_dir>/module_2b_coassembly`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_manifest` | `null` | Module 1 trimmed-read manifest; defaults to Module 1's `summary/trimmed_manifest.tsv`. |
| `coassembly_groups` | `null` | Required two-column tab-separated read/sample-to-group manifest. |
| `megahit_version` | `1.2.9` | MEGAHIT version to install/use. |
| `auto_install` | `true` | Install MEGAHIT automatically when it is unavailable. |
| `tool_env_dir` | `null` | Optional custom environment directory for Module 2B tools. |
| `threads` | `null` | Global thread override. |
| `assembly_threads` | `4` | Assembly threads when `--threads` is unset. |
| `megahit_threads` | `null` | Optional MEGAHIT-specific thread override. |
| `memory_gb` | `0` | Assembly memory limit in GB; `0` auto-detects available task memory with safety headroom. |
| `megahit_preset` | `meta-large` | MEGAHIT preset. |
| `publish_assemblies_mode` | `symlink` | Publication mode for coassembly FASTAs: `symlink`, `copy`, or `move`. |
| `publish_coassembly_reads_mode` | `symlink` | Publication mode for concatenated coassembly reads. |
| `outdir` | `<results_dir>/module_2b_coassembly` | Module 2B output directory. |

```
*IMPORTANT*
-For the coassembly manifest, only include the reads that you want coassembled.
-Currently, only MEGAHIT is used for running co-assemblies because of memory constraints. If enough users want
metaSPAdes support for this we can certainly add it, just reach out to devs or open an issue.
```

![SAMWISE step3](images/step_3.png)

`module_3_binning.nf` performs binning of metagenome assembled genomes (MAGs). It produces a folder with all of the MAGs that can then be fed into refinement pipelines (Module 4). For this, we run 3 different binners: `quickbin`, `metabat2`, and `maxbin2`

## Usage:
```bash
nextflow run module_3_binning.nf \
--working_dir ./output_samwise \
--threads 6 \
--metabat2 \
--quickbin \
--maxbin2 \
--min_scaffold_length 2500

# use `--` for any additional flags as well
```

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 3 outputs are written to `<working_dir>/module_3_binning`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_assembly_manifest` | `null` | Input assembly manifest file. Used to provide assemblies for binning. |
| `input_trimmed_manifest` | `null` | Input trimmed-read manifest file. Used to provide trimmed reads for read mapping/coverage generation. |
| `input_coassembly_assembly_manifest` | `null` | Input coassembly assembly manifest file, typically produced by Module 2B/coassembly. |
| `input_coassembly_trimmed_manifest` | `null` | Input coassembly trimmed-read manifest file, typically used with Module 2B/coassembly outputs. |
| `include_module2` | `true` | Whether to include/use assemblies from Module 2. |
| `include_module2b` | `true` | Whether to include/use coassemblies from Module 2B. |
| `metabat2` | `false` | Enables binning with MetaBAT2. |
| `quickbin` | `false` | Enables binning with QuickBin. |
| `maxbin2` | `false` | Enables binning with MaxBin2. |
| `auto_install` | `true` | Whether to automatically install required tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `tool_env_dir` | `null` | Optional custom path for the conda environment containing Module 3 tools. |
| `threads` | `null` | Global thread override. If provided, this can be used instead of module-specific thread settings. |
| `mapping_threads` | `4` | Number of threads to use for read mapping steps, such as BBMap. |
| `binning_threads` | `4` | Number of threads to use for binning tools. |
| `seqkit_version` | `2.8.2` | Version of SeqKit to install/use. |
| `bbmap_version` | `39.81` | Version of BBMap to install/use. |
| `samtools_version` | `1.23.1` | Version of Samtools to install/use. |
| `metabat2_version` | `2.18` | Version of MetaBAT2 to install/use. |
| `maxbin2_version` | `2.2.7` | Version of MaxBin2 to install/use. |
| `min_scaffold_length` | `2500` | Minimum scaffold/contig length to keep before binning. |
| `bbmap_minid` | `0.90` | Minimum sequence identity for BBMap read mapping. |
| `bbmap_maxindel` | `10` | Maximum allowed indel length for BBMap alignments. |
| `bbmap_ambig` | `random` | How BBMap handles ambiguous read mappings. |
| `bbmap_mateqtag` | `true` | Enables BBMap mate quality tagging. |
| `bbmap_extra_args` | `""` | Additional custom arguments to pass directly to BBMap. |
| `bbmap_xmx` | `4g` | Java memory setting for BBMap, passed as an `-Xmx` value. |
| `metabat2_min_contig` | `2500` | Minimum contig length to use for MetaBAT2 binning. |
| `metabat2_extra_args` | `""` | Additional custom arguments to pass directly to MetaBAT2. |
| `maxbin2_extra_args` | `""` | Additional custom arguments to pass directly to MaxBin2. |
| `quickbin_mincluster` | `50k` | Minimum cluster size setting for QuickBin. |
| `quickbin_mincontig` | `2500` | Minimum contig length to use for QuickBin. |
| `quickbin_minseed` | `2500` | Minimum seed contig length for QuickBin. |
| `quickbin_stringency` | `normal` | QuickBin stringency setting. |
| `quickbin_gzip` | `false` | Enables gzip-compressed QuickBin output, if supported. |
| `quickbin_chaff` | `false` | Enables QuickBin chaff-related output/handling, if supported. |
| `quickbin_clade` | `false` | Enables QuickBin clade-related output/handling, if supported. |
| `quickbin_sketch` | `false` | Enables QuickBin sketch-related output/handling, if supported. |
| `quickbin_server` | `false` | Enables QuickBin server mode/options, if supported. |
| `quickbin_xmx` | `null` | Optional Java memory setting for QuickBin. If unset, no custom QuickBin memory value is used. |
| `quickbin_extra_args` | `""` | Additional custom arguments to pass directly to QuickBin. |
| `quickbin_use_positional_bam` | `false` | Whether to use positional BAM input behavior for QuickBin, if supported by the workflow. |
| `publish_filtered_assemblies_mode` | `symlink` | How filtered assembly files are published to the output directory. Options can be `symlink`, `copy`, or `move`.|
| `publish_bam_mode` | `symlink` | How BAM mapping files are published. Options can be `symlink`, `copy`, or `move`.| 
| `publish_bins_mode` | `symlink` | How final bin files are published. Options can be `symlink`, `copy`, or `move`.|
| `results_dir` | null | Internal results root, derived from the resolved `working_dir`. Usually does not need to be set directly.|
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `module2_outdir` | `<results_dir>/module_2_readassembly` | Expected Module 2 output directory. Usually derived automatically and does not need to be set directly. |
| `module2b_outdir` | `<results_dir>/module_2b_coassembly` | Expected Module 2B/coassembly output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_3_binning` | Module 3 output directory. Usually derived automatically and does not need to be set directly. |

![SAMWISE step4](images/step_4.png)

`module_4_binrefinement.nf` collects Module 3 bins, predicts genes with Prodigal, runs HMMER marker searches, runs MAGScoT, and reconstructs refined MAG FASTAs. Its primary downstream contract is `summary/magscot_refined_bins_manifest.tsv` together with `refined_bins/`.

MAGScoT is intentionally best-effort: missing usable marker/contig input, no refined assignments, or a nonzero MAGScoT exit produce a status/summary rather than an immediate workflow failure. Review the Module 4 summary and refined-bin manifest before continuing; Module 5 and Module 6 require usable refined MAGs.

## Usage:
```bash
nextflow run module_4_binrefinement.nf \
--working_dir ./output_samwise \
--magscot_threshold 0

#MAGScoT original code sets this threshold at 0.5, but since we are doing gtdb + checkm runs after on the latest databases, its better to just pass this as default 0 and retain all possible MAGs. Feel free to change that --magscot_threshold param to 0.5

# use `--` for any additional flags as well
```

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 4 outputs are written to `<working_dir>/module_4_binrefinement`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_binning_manifest` | `null` | Input binning manifest file, typically produced by Module 3. This should describe the bins to be refined. |
| `dependencies_dir` | `${samwise_dir}/dependencies` | Directory containing external dependency files used by Module 4. |
| `tigrfam_hmm` | `null` | Path to the TIGRFAM HMM database file. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `pfam_hmm` | `null` | Path to the Pfam HMM database file. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `magscot_script` | `null` | Path to the MAGSCOT script. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `magscot_profiles_dir` | `null` | Path to the MAGSCOT profiles directory. |
| `auto_install` | `true` | Whether to automatically install required tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `tool_env_dir` | `null` | Optional custom path for the conda environment containing Module 4 tools. |
| `threads` | `null` | Global thread override. If provided, this can be used instead of module-specific thread settings. |
| `hmm_threads` | `8` | Number of threads to use for HMM-related steps. Used if `threads` is not provided. |
| `r_base_version` | `null` | Version of `r-base` to install/use. If `null`, the environment/tool setup may use its default version. |
| `hmmer_version` | `null` | Version of HMMER to install/use. If `null`, the environment/tool setup may use its default version. |
| `prodigal_version` | `null` | Version of Prodigal to install/use. If `null`, the environment/tool setup may use its default version. |
| `parallel_version` | `null` | Version of GNU Parallel to install/use. If `null`, the environment/tool setup may use its default version. |
| `magscot_extra_args` | `""` | Additional custom arguments to pass directly to MAGSCOT. |
| `magscot_threshold` | `0` | MAGSCOT threshold value used during bin refinement/scoring. |
| `publish_gathered_bins_mode` | `copy` | How gathered bin files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_refined_bins_mode` | `copy` | How refined bin files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `results_dir` | Derived | Internal results root, derived from the resolved `working_dir`. Usually does not need to be set directly. |
| `module3_outdir` | `<results_dir>/module_3_binning` | Expected Module 3 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_4_binrefinement` | Module 4 output directory. Usually derived automatically and does not need to be set directly. |

![SAMWISE step5](images/step_5.png)

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
--megahit true \
--metaspades true \
--secondpass_metabat2 true \
--secondpass_quickbin true \
--secondpass_maxbin2 true \
--run_second_pass_binning_refinement true

# --run_second_pass_binning_refinement specifies whether or not you want it to re-bin after subassembly - some users may want to disable this if they want to make sure subassemblies are worth performing after looking at the assembly stats, but most should leave on. Default is true. When disabled, Module 5 still writes final_mag_database using the original Module 4 refined MAGs.

# use `--` for any additional flags as well\
```
| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 5 outputs are written to `<working_dir>/module_5_subassembly`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_trimmed_manifest` | `null` | Input trimmed-read manifest file, typically produced by Module 1. |
| `input_original_binning_manifest` | `null` | Input original binning manifest file, typically produced by Module 3. |
| `input_refined_manifest` | `null` | Input refined-bin manifest file, typically produced by Module 4. |
| `megahit` | `true` | Enables subtractive assembly with MEGAHIT. At least one subtractive assembler must be enabled. |
| `metaspades` | `false` | Enables subtractive assembly with metaSPAdes. It can be enabled together with MEGAHIT. |
| `auto_install` | `true` | Whether to automatically install required tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `tool_env_dir` | `null` | Optional custom path for the conda environment containing Module 5 tools. |
| `threads` | `null` | Global thread override. If provided, this can be used instead of module-specific thread settings. |
| `mapping_threads` | `4` | Number of threads to use for read mapping steps, such as BBMap. |
| `assembly_threads` | `4` | Number of threads to use for subtractive assembly if `threads` is not provided. |
| `bbmap_version` | `39.81` | Version of BBMap to install/use. |
| `megahit_version` | `1.2.9` | Version of MEGAHIT to install/use. |
| `spades_version` | `4.2.0` | Version of SPAdes/metaSPAdes to install/use. |
| `bbmap_extra_args` | `""` | Additional custom arguments to pass directly to BBMap. |
| `bbmap_minid` | `0.99` | Minimum sequence identity for BBMap mapping during subtraction. |
| `bbmap_ambig` | `random` | How BBMap handles ambiguous read mappings. |
| `bbmap_xmx` | `null` | Optional Java memory setting for BBMap. If unset, no custom BBMap memory value is used. |
| `megahit_preset` | `meta-large` | MEGAHIT preset to use for assembly. Options can be `meta-large` or `meta-sensitive`. |
| `megahit_threads` | `null` | Optional MEGAHIT-specific thread override. If provided, this overrides the general assembly thread setting for MEGAHIT. |
| `metaspades_memory_gb` | `0` | Memory limit in GB passed to metaSPAdes. Use `0` to leave metaSPAdes memory unset. |
| `run_second_pass_binning_refinement` | `true` | Whether to run a second-pass binning and refinement workflow after subtractive assembly. If `false`, Module 5 publishes the original Module 4 refined MAGs to its `final_mag_database` output. |
| `secondpass_metabat2` | `true` | Enables MetaBAT2 during second-pass binning. |
| `secondpass_quickbin` | `true` | Enables QuickBin during second-pass binning. |
| `secondpass_maxbin2` | `true` | Enables MaxBin2 during second-pass binning. |
| `module3_script` | `${samwise_dir}/module_3_binning.nf` | Path to the Module 3 binning Nextflow script used for second-pass binning. |
| `module4_script` | `${samwise_dir}/module_4_binrefinement.nf` | Path to the Module 4 bin refinement Nextflow script used for final/second-pass refinement. |
| `nextflow_exe` | `nextflow` | Nextflow executable used to launch nested/second-pass workflows. |
| `secondpass_working_dir` | `null` | Optional custom working directory for second-pass binning/refinement. If unset, defaults to `<outdir>/second_pass`. |
| `final_joint_working_dir` | `null` | Optional custom working directory for final joint refinement. If unset, defaults to `<outdir>/final_joint_refinement`. |
| `dependencies_dir` | `${samwise_dir}/dependencies` | Directory containing external dependency files used by downstream refinement steps. |
| `tigrfam_hmm` | `null` | Path to the TIGRFAM HMM database file. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `pfam_hmm` | `null` | Path to the Pfam HMM database file. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `magscot_script` | `null` | Path to the MAGSCOT script. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `magscot_extra_args` | `""` | Additional custom arguments to pass directly to MAGSCOT during final refinement. |
| `publish_reference_mode` | `copy` | How reference files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_unmapped_mode` | `symlink` | How unmapped read files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_assemblies_mode` | `symlink` | How final assembly files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_final_mags_mode` | `copy` | How final MAG files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results root, derived from the resolved `working_dir`. Usually does not need to be set directly. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `module3_outdir` | `<results_dir>/module_3_binning` | Expected Module 3 output directory. Usually derived automatically and does not need to be set directly. |
| `module4_outdir` | `<results_dir>/module_4_binrefinement` | Expected Module 4 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_5_subassembly` | Module 5 output directory. Usually derived automatically and does not need to be set directly. |
| `secondpass_dir` | `<outdir>/second_pass` | Derived second-pass output directory. Uses `secondpass_working_dir` if provided. Usually does not need to be set directly. |
| `final_joint_dir` | `<outdir>/final_joint_refinement` | Derived final joint refinement output directory. Uses `final_joint_working_dir` if provided. Usually does not need to be set directly. |

![SAMWISE step6](images/step_6.png)

`module_6_magannotate.nf` prepares normalized MAG FASTAs and runs the selected final quality and annotation tools. It prefers Module 5's `final_mag_database` and its manifest; if that is unavailable, it falls back to Module 4 `refined_bins`.

The module can run CheckM2 quality assessment, dRep dereplication, GTDB-Tk taxonomy, EggNOG-mapper functional annotation, and microTrait trait prediction. dRep requires `--run_checkm2 true`, because it uses the CheckM2 quality report as genome information. Only selected tools are installed and run. Database paths can be supplied explicitly or downloaded into Module 6's database directory when automatic download is enabled.

## Usage:
```bash
nextflow run module_6_magannotate.nf \
--working_dir ./output_samwise \
--run_checkm2 true \
--run_gtdbtk true \
--run_eggnog true \
--threads 32

# IMPORTANT: pplacer + GTDB-Tk can overallocate memory on HPC systems.
# Module 6 therefore defaults --gtdbtk_pplacer_cpus to 1. Increase it
# cautiously if the scheduler allocation has adequate memory.

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

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 6 outputs are written to `<working_dir>/module_6_magannotate`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_mag_dir` | `null` | Directory containing MAG FASTA files to annotate. Nextflow stages this directory as a formal workflow input. If unset, Module 6 prefers Module 5's `final_mag_database`, then falls back to Module 4 refined bins. |
| `input_mag_manifest` | `null` | Input MAG manifest describing MAG files to annotate. Nextflow stages this file as a formal workflow input. If unset and Module 5's final MAG database is selected, Module 6 automatically uses Module 5's `summary/final_mag_database_manifest.tsv`. |
| `mag_extension` | `fa` | Required prepared-MAG extension. Module 6 normalizes all inputs to uncompressed `.fa`; this value must remain `fa`. |
| `threads` | `null` | Total number of threads to use. If unset, tool-specific defaults may be used. |
| `tool_env_dir` | `null` | Optional custom path for the conda/mamba environment containing Module 6 tools. |
| `auto_install` | `true` | Whether to automatically install required tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `conda_pkgs_dir` | `null` | Optional workflow-local conda/mamba package cache directory. |
| `run_drep` | `true` | Whether to dereplicate MAGs with dRep. Requires `run_checkm2=true`. |
| `run_checkm2` | `true` | Whether to run CheckM2 for MAG quality assessment. |
| `run_gtdbtk` | `true` | Whether to run GTDB-Tk for taxonomic classification. |
| `run_eggnog` | `false` | Whether to run EggNOG-mapper for functional annotation. |
| `run_microtrait` | `false` | Whether to run microTrait genomic trait prediction. |
| `drep_version` | `null` | Version of dRep to install/use. |
| `drep_env_dir` | `null` | Optional custom environment directory for dRep. |
| `drep_threads` | `null` | dRep thread count; defaults to the global `threads` setting when unset. |
| `drep_extra_args` | `-sa 0.99 -comp 50 -con 10` | Additional arguments passed to dRep. |
| `checkm2_version` | `null` | Version of CheckM2 to install/use. If `null`, the environment/tool setup may use its default version. |
| `gtdbtk_version` | `2.7.2` | Version of GTDB-Tk to install/use. |
| `checkm2_db_path` | `null` | Path to an existing CheckM2 database file. If provided, this database is used directly. |
| `checkm2_db_dir` | `null` | Directory where the CheckM2 database should be stored or checked. If unset, defaults to `<outdir>/databases/checkm2`. |
| `checkm2_zenodo_record` | `14897628` | Zenodo record ID used for downloading the CheckM2 database. |
| `checkm2_auto_download_db` | `true` | Whether to automatically download the CheckM2 database if it is missing. |
| `checkm2_extension` | `fa` | Required extension for prepared MAGs passed to CheckM2; must remain `fa`. |
| `gtdbtk_data_path` | `null` | Path to an existing GTDB-Tk database directory. If provided, this database is used directly. |
| `gtdbtk_db_dir` | `null` | Directory where the GTDB-Tk database should be stored or checked. If unset, defaults to `<outdir>/databases/gtdbtk`. |
| `gtdbtk_auto_download_db` | `true` | Whether to automatically download the GTDB-Tk database if it is missing. |
| `gtdbtk_extension` | `fa` | Required extension for prepared MAGs passed to GTDB-Tk; must remain `fa`. |
| `gtdbtk_download_url` | `https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz` | URL used to download the GTDB-Tk database package. |
| `eggnog_env_dir` | `null` | Optional custom environment directory for EggNOG-mapper. |
| `eggnog_mapper_version` | `2.1.13` | Version of EggNOG-mapper to install/use. |
| `eggnog_data_path` | `null` | Path to an existing EggNOG-mapper data directory. If provided, this database/data path is used directly. |
| `eggnog_data_dir` | `null` | Directory where the EggNOG-mapper database should be stored or checked. If unset, defaults to `<outdir>/databases/eggnog`. |
| `eggnog_auto_download_db` | `true` | Whether to automatically download the EggNOG-mapper database if it is missing. |
| `eggnog_download_args` | `-y` | Arguments passed to the EggNOG-mapper database download command. |
| `eggnog_fixurl` | `true` | Whether to use/apply the EggNOG-mapper URL fix package during setup. |
| `eggnog_fixurl_package` | `eggnog-mapper-fixurl` | Package used to fix EggNOG-mapper database download URLs. |
| `eggnog_method` | `diamond` | Search method used by EggNOG-mapper. Common options include `diamond` and other EggNOG-supported methods. |
| `eggnog_itype` | `metagenome` | Input type passed to EggNOG-mapper. |
| `eggnog_genepred` | `prodigal` | Gene prediction method used by EggNOG-mapper. |
| `eggnog_trans_table` | `11` | Translation table used for gene prediction/annotation. |
| `eggnog_output_prefix` | `samwise_eggnog` | Prefix used for EggNOG-mapper output files. |
| `eggnog_extra_args` | `""` | Additional custom arguments to pass directly to EggNOG-mapper. |
| `eggnog_mmseqs_db` | `null` | Optional MMseqs database path for EggNOG-mapper workflows that use MMseqs. |
| `eggnog_fail_nonfatal` | `false` | If enabled, EggNOG-mapper failures are treated as non-fatal, allowing the workflow to continue. |
| `microtrait_env_dir` | `null` | Optional custom microTrait environment directory. |
| `microtrait_type` | `genomic` | microTrait analysis type. |
| `microtrait_output_prefix` | `samwise_microtrait` | Prefix for microTrait outputs. |
| `microtrait_auto_download_db` | `true` | Whether microTrait dependencies/databases may be downloaded during setup. |
| `microtrait_fail_nonfatal` | `false` | If enabled, an unrecoverable microTrait failure is recorded without failing the complete workflow. |
| `publish_tool_outputs_mode` | `copy` | How annotation and quality-control tool outputs are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results root, derived from the resolved `working_dir`. Usually does not need to be set directly. |
| `outdir` | `<results_dir>/module_6_magannotate` | Module 6 output directory. Usually derived automatically and does not need to be set directly. |
| `module5_final_mag_dir` | `<results_dir>/module_5_subassembly/final_mag_database` | Candidate MAG directory from Module 5 subtractive assembly. |
| `module5_final_mag_manifest` | `<results_dir>/module_5_subassembly/summary/final_mag_database_manifest.tsv` | Module 5 final-MAG manifest used automatically when Module 5's final MAG directory is selected. |
| `module4_refined_mag_dir` | `<results_dir>/module_4_binrefinement/refined_bins` | Candidate refined MAG directory from Module 4 bin refinement. |
| `checkm2_db_outdir` | `<outdir>/databases/checkm2` | Derived CheckM2 database output directory. Uses `checkm2_db_dir` if provided. Usually does not need to be set directly. |
| `gtdbtk_db_outdir` | `<outdir>/databases/gtdbtk` | Derived GTDB-Tk database output directory. Uses `gtdbtk_db_dir` if provided. Usually does not need to be set directly. |
| `eggnog_db_outdir` | `<outdir>/databases/eggnog` | Derived EggNOG-mapper database output directory. Uses `eggnog_data_path` or `eggnog_data_dir` if provided. Usually does not need to be set directly. |

![SAMWISE step7](images/step_7.png)

# Module 7: Genome-Scale Metabolic Model Generation

Module 7 constructs one genome-scale metabolic model (GEM) for each MAG represented in Module 6 eggNOG predictions. It uses gapseq to reconstruct and gap-fill models, can adapt selected models to MAG-specific empirical growth requirements, and validates each final SBML model with MEMOTE.

The workflow implementation is [`module_7_gems.nf`](module_7_gems.nf). It consumes Module 6 results by default and publishes a per-MAG output manifest containing stable result paths rather than temporary Nextflow work-directory paths.

## What this module does

1. **Creates or validates the pinned gapseq/MEMOTE environment.**
2. **Prepares per-MAG protein FASTAs** by splitting the Module 6 unified eggNOG predicted-protein FASTA according to MAG IDs in the Module 6 manifest.
3. **Runs gapseq `doall`** for every MAG using the selected medium and bacterial or archaeal template.
4. **Optionally runs gapseq `adapt`** for only the MAGs named in an adaptation manifest; unlisted MAGs retain their `doall` model.
5. **Runs optional MEMOTE validation** on each final model as an HTML snapshot report, JSON test result, or both.
6. **Writes summary files** describing every input, model, report, and log produced for each MAG.

## Usage

### Default run

This expects Module 6 outputs below the working directory, uses the bundled comprehensive medium and `Bacteria` template, and creates one MEMOTE HTML snapshot report per final model.

```bash
nextflow run module_7_gems.nf \
  --working_dir ./output_samwise
```

### Use an archaeal template and custom medium

```bash
nextflow run module_7_gems.nf \
  --working_dir ./output_samwise \
  --template_organism Archaea \
  --media_csv /absolute/path/to/media.csv
```

### Adapt selected MAGs and generate both MEMOTE outputs

```bash
nextflow run module_7_gems.nf \
  --working_dir ./output_samwise \
  --run_gapseq_adapt true \
  --adapt_manifest /absolute/path/to/adapt_manifest.tsv \
  --memote_mode both
```

### Run JSON MEMOTE tests only

```bash
nextflow run module_7_gems.nf \
  --working_dir ./output_samwise \
  --memote_mode run
```

### Build GEMs without MEMOTE validation

```bash
nextflow run module_7_gems.nf \
  --working_dir ./output_samwise \
  --run_memote false
```

### Use Module 6 files stored elsewhere

`protein_fasta_dir` is a historic parameter name. It accepts the path to a unified predicted-protein FASTA **file**, not a directory.

```bash
nextflow run module_7_gems.nf \
  --working_dir ./module_7_results \
  --protein_fasta_dir /absolute/path/to/samwise_eggnog.emapper.genepred.fasta \
  --input_manifest /absolute/path/to/eggnog_input_manifest.tsv
```

## Inputs and defaults

### Module 6 inputs

When `working_dir` is supplied, Module 7 reads these files by default:

| Input | Default path |
|---|---|
| Unified eggNOG predicted-protein FASTA | `<working_dir>/module_6_magannotate/eggnog/samwise_eggnog.emapper.genepred.fasta` |
| eggNOG input manifest | `<working_dir>/module_6_magannotate/summary/eggnog_input_manifest.tsv` |

If only `samwise_dir` is provided, `working_dir` inherits that same directory. Override one or both defaults with `--protein_fasta_dir` and `--input_manifest`.

The unified FASTA and manifest must be present and non-empty. The manifest must include a non-empty `mag_id` column. Module 7 associates proteins with MAGs from FASTA headers: the text before the first `|` (or whitespace, if there is no `|`) must equal a `mag_id` in the manifest. MAGs with no matching proteins are recorded as warnings; the run fails only if no proteins match any MAG.

### Media and template

- `--media_csv` selects the CSV medium passed to gapseq. It must be an existing, non-empty file.
- If `media_csv` is not given, Module 7 uses [`background/media/gapseq_all_nutrients.csv`](background/media/gapseq_all_nutrients.csv), the bundled comprehensive medium.
- The bundled minimal M9 glucose aerobic medium is [`background/media/gapseq_M9_glucose_aerobic.csv`](background/media/gapseq_M9_glucose_aerobic.csv). Select it explicitly with `--media_csv`.
- `--template_organism` must be exactly `Bacteria` (default) or `Archaea`.

### MAG-specific adaptation manifest

Adaptation is off by default. To enable it, both of the following are required:

```bash
--run_gapseq_adapt true \
--adapt_manifest /absolute/path/to/adapt_manifest.tsv
```

The adaptation manifest is tab-separated and must contain `mag_id` and `adapt_compounds` headers. It may have additional columns. Each data row requires non-empty values, each `mag_id` must occur in the Module 6 manifest exactly once, and every comma-separated compound value must use `cpd#####:(TRUE|FALSE)` syntax.

```tsv
mag_id	adapt_compounds
MAGScoT_cleanbin_000023	cpd00076:TRUE,cpd00027:TRUE
MAGScoT_cleanbin_000038	cpd00027:TRUE
```

[`background/test_adapt_manifest.tsv`](background/test_adapt_manifest.tsv) provides a format example. Only MAGs listed in the manifest enter the adaptation step. If gapseq reports that a model already grows and produces no `-adapt` files, Module 7 publishes copies of the original `doall` model as that MAG's adaptation output and records the `no_changes` adaptation status in its task metadata.

`adaptation_compounds` has been retired and is rejected if supplied. A global compound list is not supported because adaptation requirements are MAG-specific.

## Environment and tool versions

Module 7 requires these fixed versions:

| Tool | Required version |
|---|---|
| gapseq | `1.4.0` |
| MEMOTE | `0.17.0` |

The workflow confirms both executable versions before it constructs models. The package pins are internal workflow settings and are not user-configurable.

With no environment override, Module 7 manages an environment at `<outdir>/conda_envs/gapseq`, where `<outdir>` is `<results_dir>/module_7_gems`. If it is missing and `--auto_install true` (the default), the workflow creates it with `mamba` or `conda`. A stale workflow-managed environment can be removed and rebuilt automatically.

Provide either `--gapseq_env_dir` or `--tool_env_dir` to use a user-managed environment:

- `--gapseq_env_dir /path/to/environment` uses that exact environment path.
- `--tool_env_dir /path/to/tool-parent` uses `/path/to/tool-parent/gapseq`.

A user-managed environment must already exist and provide exactly the required gapseq and MEMOTE versions. Module 7 never creates, deletes, or modifies it. `gapseq_env_dir` takes precedence if both overrides are provided. `--conda_pkgs_dir` optionally chooses a package cache for the workflow-managed environment.

The workflow uses these pinned command forms:

```bash
gapseq doall <protein.faa> <media.csv> <Bacteria|Archaea>

gapseq adapt \
  -m <model.RDS> \
  -w <cpd#####:TRUE|FALSE,...> \
  -c <model-rxnWeights.RDS> \
  -g <model-rxnXgenes.RDS> \
  -b <model-all-Reactions.tbl> \
  -f <adapt-output-dir>

memote report snapshot --filename <report.html> <model.xml>
memote run --ignore-git --filename <result.json> <model.xml>
```

The `gapseq doall` interface is positional. The `threads` parameter asks the Nextflow scheduler for CPUs for each `doall` task, but Module 7 does not translate it into a gapseq command-line argument. `gapseq_extra_args` is retired and fails preflight validation if supplied.

## MEMOTE modes

MEMOTE is enabled by default with `--memote_mode snapshot`.

| `memote_mode` | Command | Output per final model |
|---|---|---|
| `snapshot` | `memote report snapshot` | HTML report |
| `run` | `memote run --ignore-git` | JSON result |
| `both` | Both commands | HTML report and JSON result |

Set `--run_memote false` to skip every MEMOTE task. `--memote_extra_args` appends options to whichever MEMOTE command is selected; supply only options that MEMOTE `0.17.0` supports.

## Parameters

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 7 outputs are published to `<working_dir>/module_7_gems`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `input_manifest` | `null` | Override path for the Module 6 eggNOG input manifest TSV. |
| `protein_fasta_dir` | `null` | Override path for the Module 6 unified predicted-protein FASTA; this is a file path despite the name. |
| `media_csv` | `null` | User medium CSV. The bundled comprehensive medium is used when unset. |
| `template_organism` | `Bacteria` | gapseq template: exactly `Bacteria` or `Archaea`. |
| `run_gapseq_adapt` | `false` | Enable adaptation for MAGs listed in `adapt_manifest`. |
| `adapt_manifest` | `null` | Adaptation TSV. Required when `run_gapseq_adapt` is `true`. |
| `adaptation_compounds` | retired | Do not supply. A global adaptation-compound list is rejected; use `adapt_manifest` for MAG-specific requirements. |
| `run_memote` | `true` | Enable MEMOTE validation for final models. |
| `memote_mode` | `snapshot` | `snapshot`, `run`, or `both`. |
| `threads` | `null` | CPU count requested from Nextflow for each gapseq `doall` task. If `threads` is unset or `null`, Nextflow defaults allocation to 1 thread per `doall` process. |
| `tool_env_dir` | `null` | Parent directory for a user-managed `gapseq` environment. |
| `gapseq_env_dir` | `null` | Exact path to a user-managed environment; takes precedence over `tool_env_dir`. |
| `auto_install` | `true` | Create a missing workflow-managed environment with `mamba` or `conda`. |
| `conda_pkgs_dir` | `null` | Optional conda/mamba package cache for the managed environment. |
| `gapseq_extra_args` | retired | Do not supply. gapseq `doall` uses a fixed positional interface in this workflow. |
| `memote_extra_args` | `""` | Additional options appended to the selected MEMOTE command. |
| `publish_gems_mode` | `copy` | Publish mode for prepared FASTAs and GEM directories: `copy`, `symlink`, or `move`. |
| `publish_reports_mode` | `copy` | Publish mode for MEMOTE reports: `copy`, `symlink`, or `move`. |
| `results_dir` | derived | Internal results root, derived from the resolved `working_dir`. Normally do not set directly. |
| `outdir` | `<results_dir>/module_7_gems` | Internal Module 7 output location. Normally do not set directly. |
| `module6_output_dir` | `<results_dir>/module_6_magannotate` | Default Module 6 output location used to derive input paths. |
| `media_minimal` | bundled M9/glucose CSV | Internal path for the bundled minimal medium; use `media_csv` to select it. |
| `media_comprehensive` | bundled all-nutrients CSV | Internal path for the bundled comprehensive medium used by default. |

## Outputs

All published outputs are beneath `<results_dir>/module_7_gems/`:

```text
module_7_gems/
├── inputs/
│   ├── protein_fastas/<mag_id>.faa
│   ├── gapseq_input_manifest.tsv
│   └── gapseq_inputs_stats.tsv
├── setup/
│   └── gapseq_setup_status.env
├── gems/
│   ├── doall/<mag_id>_gapseq_doall/
│   └── adapt/<mag_id>_gapseq_adapt/             # selected MAGs only
├── reports/
│   ├── snapshot/<mag_id>_memote_snapshot.html   # snapshot or both
│   └── run/<mag_id>_memote_run.json             # run or both
├── logs/
└── summary/
    ├── module_7_gems_manifest.tsv
    └── module_7_gems_summary.tsv
```

`logs/` contains the execution logs for `PREPARE_GAPSEQ_INPUTS`, `SETUP_GAPSEQ`, `RUN_GAPSEQ_DOALL`, `RUN_GAPSEQ_ADAPT`, `RUN_MEMOTE_*`, and `WRITE_OUTPUT_MANIFEST`.

`summary/module_7_gems_manifest.tsv` has one row per MAG. It records the prepared protein FASTA, `doall` model and support-artifact paths, optional adaptation paths, selected final model stage and paths, optional MEMOTE output paths, and associated logs.

`summary/module_7_gems_summary.tsv` reports the total MAG count, adapted MAG count, doall-only MAG count, and counts of snapshot- and full-run-validated MAGs.

## Validation

Before production use, run Module 7 on a small representative Module 6 output set and inspect the generated `summary/module_7_gems_manifest.tsv`, gapseq logs, and MEMOTE reports. This confirms that the installed gapseq and MEMOTE versions, selected medium, and available compute resources are appropriate for the local environment.

![SAMWISE AuxModules](images/AuxModules.png)

## Auxiliary Module 1: Assembly annotation

`AuxModule_1_assemblyAnnotate.nf` functionally annotates assemblies rather than MAGs. It discovers available assembly directories from Module 2, Module 2B, and Module 5; filters scaffolds shorter than `--min_scaffold_bp`; removes exact duplicate and reverse-complement-equivalent full-length scaffolds; then runs EggNOG-mapper on the combined retained sequences. At least one of those assembly-producing modules must have been run.

```bash
nextflow run AuxModule_1_assemblyAnnotate.nf \
--working_dir ./output_samwise \
--threads 36 \
--min_scaffold_bp 1000
```

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root containing Module 2, Module 2B, and/or Module 5 outputs. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `min_scaffold_bp` | `1000` | Minimum scaffold length retained for annotation. |
| `threads` | `null` | Global EggNOG thread override. |
| `auto_install` | `true` | Install required tools when unavailable. |
| `eggnog_data_path` | `null` | Existing EggNOG data directory to use. |
| `eggnog_data_dir` | `null` | EggNOG data directory to check or populate. |
| `eggnog_auto_download_db` | `true` | Download EggNOG data when no valid configured or Module 6 database is available. |
| `eggnog_method` | `diamond` | EggNOG-mapper search method. |
| `eggnog_fail_nonfatal` | `false` | Continue with a recorded status if EggNOG fails. |
| `outdir` | `<working_dir>/AuxModule_1_assemblyAnnotate` | Auxiliary Module 1 output directory. |

Auxiliary Module 1 reuses a valid Module 6 EggNOG database at `module_6_magannotate/databases/eggnog` before downloading another copy. Results include `filtered_assemblies/`, `inputs/`, `eggnog/`, and summary tables under `summary/`.

## Auxiliary Module 2: MVP viral analysis

`AuxModule_2_mvp.nf` prepares Module 2 individual assemblies, optional Module 2B coassemblies, and optional Module 5 subtractive assemblies for the MVP viral workflow. It derives compatible metadata and read inputs, interleaving paired reads with BBTools only when required, then invokes selected MVP modules in canonical order. Module 5 is optional: absent or empty subtractive manifests are skipped safely.

```bash
nextflow run AuxModule_2_mvp.nf \
--working_dir ./output_samwise \
--mvp_modules "0,1,2,3,4,5,100" \
--threads 36 \
--memory_gb 0
```

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it inherits `working_dir`. |
| `mvp_modules` | `0,1,2,3,4,5,100` | Comma-, semicolon-, or whitespace-separated MVP stages; valid values are `0,1,2,3,4,5,6,7,99,100`. |
| `include_individual_assemblies` | `true` | Include Module 2 assemblies. |
| `include_coassemblies` | `true` | Include Module 2B assemblies when present. |
| `include_subtractive_assemblies` | `true` | Include Module 5 assemblies when present. |
| `assembly_manifest` | `null` | Override Module 2 assembly manifest. |
| `trimmed_manifest` | `null` | Override Module 1 trimmed-read manifest. |
| `coassembly_manifest` | `null` | Override Module 2B assembly manifest. |
| `subtractive_assembly_manifest` | `null` | Override Module 5 subtractive-assembly manifest. |
| `threads` | `4` | MVP and preparation thread count. |
| `memory_gb` | `0` | Memory value passed to MVP; `0` delegates memory selection to MVP. |
| `install_databases` | `false` | Allow MVP to install its databases. |
| `genomad_db_path` | `null` | Existing geNomad database path. |
| `checkv_db_path` | `null` | Existing CheckV database path. |
| `outdir` | `<results_dir>/AuxModule_2_mvp` | MVP output directory. |

MVP stages are not reordered by the supplied list. If an earlier MVP stage is omitted, any outputs it requires must already exist under the Module 2 output directory from a prior MVP run.

BETA AI AGENT: 

If you would like to test out the AI Agent that can help you interrogate your genomes and their metabolisms, simply set up the OpenAI agent by running:

`bash setupAgentOpenai.sh`

This will generate a conda environment called `langchain-chat-openai` that you then need to activate with `conda activate langchain-chat-openai`. 

After this, you need to set up your environment file (see env.txt example in samwise-main/agent/) to include your API key as well as user settings. Once that is done, change the file name to .env instead of env.txt so that the agent can find it. 

Then, you can run the agent using `python chatOpenai.py` within the /agent/ folder (you need to cd /agent/ if you have not already).

_____________________________________________________________________________________
DISCLAIMER
This material was prepared as an account of work sponsored by an agency of the
United States Government.  Neither the United States Government nor the United
States Department of Energy, nor Battelle, nor any of their employees, nor any
jurisdiction or organization that has cooperated in the development of these
materials, makes any warranty, express or implied, or assumes any legal
liability or responsibility for the accuracy, completeness, or usefulness or
any information, apparatus, product, software, or process disclosed, or
represents that its use would not infringe privately owned rights.
 
Reference herein to any specific commercial product, process, or service by
trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or any agency thereof, or Battelle Memorial Institute. The
views and opinions of authors expressed herein do not necessarily state or
reflect those of the United States Government or any agency thereof.
 
                 PACIFIC NORTHWEST NATIONAL LABORATORY
                              operated by
                                BATTELLE
                                for the
                   UNITED STATES DEPARTMENT OF ENERGY
                    under Contract DE-AC05-76RL01830
_____________________________________________________________________________________
LICENSE
Copyright Battelle Memorial Institute 2026
 
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 
1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 
2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
 
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
