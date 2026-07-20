![SAMWISE title](images/SAMWISE_title.png)

# Welcome to SAMWISE!

SAMWISE is an automated, end-to-end metagenomic read processing program. Here is a quick conceptual rundown of what this software can enable you to do via Nextflow DSL2 workflows.

![SAMWISE workflow](images/SAMWISE_FULL-git.png)

---

## Requirements

- [Nextflow](https://www.nextflow.io/)
- `mamba` (or `conda`) for automatic package installation
- `Java`, required by Nextflow (auto installed if installing with `mamba` / `conda`

To start with SAMWISE, you will want to make sure that you have `mamba` (or `conda`) installed. We recommend mamba, and you can follow the instructions here: https://conda-forge.org/download/

Then, you need to install NextFlow - this can be done via `mamba` / `conda`: https://anaconda.org/channels/bioconda/packages/nextflow/overview
We recommend that you install NextFlow into its own, separate environment from your base environment. For example, with `mamba install -n nextflow -c bioconda nextflow` Then, when running SAMWISE, make sure that you activate your NextFlow environment with `mamba activate nextflow`!

Once NextFlow is isntalled, go ahead and clone this repo or download it / extract. You can click on `clone repo` in the top right on GitHub or just download the whole thing. Then, change directory into the directory of the cloned repo: `cd ./samwise-main`

`samwise-main` is what will hold all of the nextflow .nf files, and is what we recommend get set as the `--working_dir` flag. SAMWISE will auto-generate all module folders as needed.

Now, you are ready to proceed with SAMWISE!

---
![SAMWISE quickstart](images/quick_start.png)
---

`"In a hole in the ground there lived a hobbit... Not a nasty, dirt..."`

Alright alright - you want to run SAMWISE quickly and do not want to read through the full docs. Here is how I would run this as an sbatch script on a server.

NOTE: Your reads MUST be in one of the naming formats (_R1, _R2, _1, _2, _interleaved) and must
have extensions (.fq or .fastq - gzipped or not gzipped is fine). See module_0 info below!

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

nextflow run module_2b_coassembly.nf \
--working_dir ./samwise-main \
--coassembly_groups ./coassembly_manifest.txt \
--threads 36 \
--memory_gb 0

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

#Note: If running only 1 binner, you need to pass the MAG manifest from Module 3 directly with: --input_mag_manifest

# To save compute time, you can pre-download the databases before you run module 6.
# For this, see: CheckM2: https://zenodo.org/records/14897628, gtdbtk: https://ecogenomics.github.io/GTDBTk/installing/index.html,
# eggnog: https://github.com/eggnogdb/eggnog-mapper; command: download_eggnog_data.py --data_dir /path/to/eggnog-data
# If you already pre-downloaded the gtdb, checkm2, and eggnog databases, 
# you can directly pass the paths as arguments:
# --checkm2_db_path /path/to/uniref100.KO.1.dmnd
# --gtdbtk_data_path /path/to/gtdbtk/database_directory
# --eggnog_data_path /path/to/eggnog/database_directory

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
| `results_dir` | null | Internal results directory. Uses `--working_dir` if provided, otherwise `--output_dir`, otherwise `.`. Usually does not need to be set directly. |
| `module0_outdir` | `<results_dir>/module_0_readprocess` | Expected Module 0 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_1_readtrimming` | Module 1 output directory. Usually derived automatically and does not need to be set directly. |

```
*IMPORTANT*
Currently, rarefied assemblies are set to run as paralell processes to single assemblies to speed things up.
In theory, they should play nice. However, if you run into issues with clobbering memory, we will be working
on adding a flag so that the rarefied assemblies run only after single assemblies are complete.
```

![SAMWISE step2](images/step_2.png)

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
| `megahit_preset` | `meta-large` | MEGAHIT preset to use for assembly. Default is `meta-large`. Could also use `meta-sensitive` |
| `publish_assemblies_mode` | `symlink` | How final assembly files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results directory. Uses `--working_dir` if provided, otherwise `--output_dir`, otherwise `.`. Usually does not need to be set directly. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_2_readassembly` | Module 2 output directory. Usually derived automatically and does not need to be set directly. |

```
*IMPORTANT*
The default assembly outputs get written into the NextFlow work directories to save space. If you want it to
write out the output assemblies into a more accessible location, you can set publish_assemblies_mode to be
`copy`. Argument `move` here would also work but may cause issues with NextFlow not finding what it needs.
```

![SAMWISE step2b](images/step_2b.png)

`module_2b_coassembly.nf` performs **grouped co-assembly** from Module 1 trimmed reads using **MEGAHIT only**. This module is designed to run alongside the normal Module 2 assembly workflow. It produces Module-3-compatible manifests so that Module 3 can automatically bin co-assemblies using the exact concatenated reads that were used to generate each co-assembly.

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
| `working_dir` | `null` | Main working/results directory for the pipeline. If provided, Module 2B outputs are written to `<working_dir>/module_2_readassembly`. |
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
| `tool_env_dir` | `null` | Optional custom path for the conda environment containing Module 2B assembly tools. |
| `threads` | `null` | Global thread override. If provided, this can be used instead of module-specific thread settings. |
| `assembly_threads` | `4` | Number of threads to use for assembly if `--threads` is not provided. |
| `memory_gb` | `0` | Global memory limit in GB for assembly processes. Use `0` to leave memory unset. |
| `megahit_threads` | `null` | Optional MEGAHIT-specific thread override. If provided, this overrides the general assembly thread setting for MEGAHIT. |
| `megahit_preset` | `meta-large` | MEGAHIT preset to use for assembly. Default is `meta-large`. |
| `publish_assemblies_mode` | `symlink` | How final assembly files are published to the output directory. Options can be `symlink`, `copy`, or `move`.|
| `results_dir` | `null` | Internal results directory. Uses `--working_dir` if provided, otherwise `--output_dir`, otherwise `.`. Usually does not need to be set directly. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_2_readassembly` | Module 2B output directory. Usually derived automatically and does not need to be set directly. |

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
| `working_dir` | `null` | Main working/results directory for the pipeline. If provided, Module 3 outputs are written to `<working_dir>/module_3_binning`. |
| `input_assembly_manifest` | `null` | Input assembly manifest file. Used to provide assemblies for binning. |
| `input_trimmed_manifest` | `null` | Input trimmed-read manifest file. Used to provide trimmed reads for read mapping/coverage generation. |
| `input_coassembly_assembly_manifest` | `null` | Input coassembly assembly manifest file, typically produced by Module 2B/coassembly. |
| `input_coassembly_trimmed_manifest` | `null` | Input coassembly trimmed-read manifest file, typically used with Module 2B/coassembly outputs. |
| `include_module2` | `true` | Whether to include/use assemblies from Module 2. |
| `include_module2b` | `true` | Whether to include/use coassemblies from Module 2B. |
| `output_dir` | `null` | Alternative output directory used only if `--working_dir` is not provided. |
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
| `results_dir` | null | Internal results directory. Uses `--working_dir` if provided, otherwise `--output_dir`, otherwise `.`. Usually does not need to be set directly.|
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `module2_outdir` | `<results_dir>/module_2_readassembly` | Expected Module 2 output directory. Usually derived automatically and does not need to be set directly. |
| `module2b_outdir` | `<results_dir>/module_2b_coassembly` | Expected Module 2B/coassembly output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_3_binning` | Module 3 output directory. Usually derived automatically and does not need to be set directly. |

![SAMWISE step4](images/step_4.png)

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

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Main working/results directory for the pipeline. If provided, Module 4 outputs are written to `<working_dir>/module_4_binrefinement`. |
| `output_dir` | `null` | Alternative output directory used only if `working_dir` is not provided. |
| `input_binning_manifest` | `null` | Input binning manifest file, typically produced by Module 3. This should describe the bins to be refined. |
| `dependencies_dir` | `${projectDir}/dependencies` | Directory containing external dependency files used by Module 4. |
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
| `results_dir` | Derived | Internal results directory. Uses `working_dir` if provided, otherwise `output_dir`, otherwise `.`. Usually does not need to be set directly. |
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
--megahit \
--metaspades \
--secondpass_metabat2 true \
--secondpass_quickbin true \
--secondpass_maxbin2 true \
--run_second_pass_binning_refinement true

# --run_second_pass_binning_refinement specifies whether or not you want it to re-bin after subassembly - some users may want to disable this if they want to make sure subassemblies are worth performing after looking at the assembly stats, but most should leave on. Default is true.

# use `--` for any additional flags as well\
```
| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Main working/results directory for the pipeline. If provided, Module 5 outputs are written to `<working_dir>/module_5_subtractiveassembly`. |
| `output_dir` | `null` | Alternative output directory used only if `working_dir` is not provided. |
| `input_trimmed_manifest` | `null` | Input trimmed-read manifest file, typically produced by Module 1. |
| `input_original_binning_manifest` | `null` | Input original binning manifest file, typically produced by Module 3. |
| `input_refined_manifest` | `null` | Input refined-bin manifest file, typically produced by Module 4. |
| `megahit` | `false` | Enables subtractive assembly with MEGAHIT. |
| `metaspades` | `false` | Enables subtractive assembly with metaSPAdes. |
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
| `metaspades_memory_gb` | `0` | Memory limit in GB for metaSPAdes. Use `0` to leave memory unset. |
| `run_second_pass_binning_refinement` | `true` | Whether to run a second-pass binning and refinement workflow after subtractive assembly. |
| `secondpass_metabat2` | `true` | Enables MetaBAT2 during second-pass binning. |
| `secondpass_quickbin` | `true` | Enables QuickBin during second-pass binning. |
| `secondpass_maxbin2` | `true` | Enables MaxBin2 during second-pass binning. |
| `module3_script` | `${projectDir}/module_3_binning.nf` | Path to the Module 3 binning Nextflow script used for second-pass binning. |
| `module4_script` | `${projectDir}/module_4_binrefinement.nf` | Path to the Module 4 bin refinement Nextflow script used for final/second-pass refinement. |
| `nextflow_exe` | `nextflow` | Nextflow executable used to launch nested/second-pass workflows. |
| `secondpass_working_dir` | `null` | Optional custom working directory for second-pass binning/refinement. If unset, defaults to `<outdir>/second_pass`. |
| `final_joint_working_dir` | `null` | Optional custom working directory for final joint refinement. If unset, defaults to `<outdir>/final_joint_refinement`. |
| `dependencies_dir` | `${projectDir}/dependencies` | Directory containing external dependency files used by downstream refinement steps. |
| `tigrfam_hmm` | `null` | Path to the TIGRFAM HMM database file. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `pfam_hmm` | `null` | Path to the Pfam HMM database file. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `magscot_script` | `null` | Path to the MAGSCOT script. If not provided, the workflow may look for it in `dependencies_dir`, depending on module logic. |
| `magscot_extra_args` | `""` | Additional custom arguments to pass directly to MAGSCOT during final refinement. |
| `publish_reference_mode` | `copy` | How reference files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_unmapped_mode` | `symlink` | How unmapped read files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_assemblies_mode` | `symlink` | How final assembly files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_final_mags_mode` | `copy` | How final MAG files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results directory. Uses `working_dir` if provided, otherwise `output_dir`, otherwise `.`. Usually does not need to be set directly. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. Usually derived automatically and does not need to be set directly. |
| `module3_outdir` | `<results_dir>/module_3_binning` | Expected Module 3 output directory. Usually derived automatically and does not need to be set directly. |
| `module4_outdir` | `<results_dir>/module_4_binrefinement` | Expected Module 4 output directory. Usually derived automatically and does not need to be set directly. |
| `outdir` | `<results_dir>/module_5_subtractiveassembly` | Module 5 output directory. Usually derived automatically and does not need to be set directly. |
| `secondpass_dir` | `<outdir>/second_pass` | Derived second-pass output directory. Uses `secondpass_working_dir` if provided. Usually does not need to be set directly. |
| `final_joint_dir` | `<outdir>/final_joint_refinement` | Derived final joint refinement output directory. Uses `final_joint_working_dir` if provided. Usually does not need to be set directly. |

![SAMWISE step6](images/step_6.png)

This module performs the following steps:

1. **Checks for dependencies and installs them if necessary** \
   -Module 2 will download required databases for each tool if needed, but arguments can be passed to directly point to dbs. \
   -Module 2 will also dereplicate genomes prior to running final characterization.
3. **Runs CheckM2, GTDB-tk, and Eggnog (or DRAM2)** \
   -We note that for right now, DRAM2 has been replaced with eggnog v2 since DRAM2 is undergoing significant development and is currently not fully installable. Once development is finished, we will update our module to incldue both DRAM2 and eggnog (and can updated to eggnog v3 once available as well).

## Usage:
```bash
nextflow run module_6_magannotate.nf \
--working_dir ./output_samwise \
--run_checkm2 true \
--run_gtdbtk true \
--run_eggnog true \
--threads 32

# IMPORTANT: For some reason, pplacer + gtdbtk do not play very well
# in HPC-like systems in the way that they try and allocate memory.
# As such, for this step pplacer thread defaults are set to 2. You
# can try and add more threads using arg. --pplacer_cpus but if
# this step fails, this is probably why.

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
| `working_dir` | `null` | Main working/results directory for the pipeline. If provided, Module 6 outputs are written to `<working_dir>/module_6_magannotation`. |
| `output_dir` | `null` | Alternative output directory used only if `working_dir` is not provided. |
| `input_mag_dir` | `null` | Directory containing MAG/bin FASTA files to annotate. |
| `input_mag_manifest` | `null` | Input MAG manifest file describing MAG/bin files to annotate. |
| `mag_extension` | `fa` | File extension used to detect MAG files in `input_mag_dir`. |
| `threads` | `null` | Total number of threads to use. If unset, tool-specific defaults may be used. |
| `tool_env_dir` | `null` | Optional custom path for the conda/mamba environment containing Module 6 tools. |
| `auto_install` | `true` | Whether to automatically install required tools using `mamba` or `conda` if they are not found. If set to `false`, required tools must already be available. |
| `conda_pkgs_dir` | `null` | Optional workflow-local conda/mamba package cache directory. |
| `run_checkm2` | `true` | Whether to run CheckM2 for MAG quality assessment. |
| `run_gtdbtk` | `true` | Whether to run GTDB-Tk for taxonomic classification. |
| `run_eggnog` | `false` | Whether to run EggNOG-mapper for functional annotation. |
| `checkm2_version` | `null` | Version of CheckM2 to install/use. If `null`, the environment/tool setup may use its default version. |
| `gtdbtk_version` | `2.7.2` | Version of GTDB-Tk to install/use. |
| `checkm2_db_path` | `null` | Path to an existing CheckM2 database file. If provided, this database is used directly. |
| `checkm2_db_dir` | `null` | Directory where the CheckM2 database should be stored or checked. If unset, defaults to `<outdir>/databases/checkm2`. |
| `checkm2_zenodo_record` | `14897628` | Zenodo record ID used for downloading the CheckM2 database. |
| `checkm2_auto_download_db` | `true` | Whether to automatically download the CheckM2 database if it is missing. |
| `checkm2_extension` | `fa` | File extension used for MAG files passed to CheckM2. |
| `gtdbtk_data_path` | `null` | Path to an existing GTDB-Tk database directory. If provided, this database is used directly. |
| `gtdbtk_db_dir` | `null` | Directory where the GTDB-Tk database should be stored or checked. If unset, defaults to `<outdir>/databases/gtdbtk`. |
| `gtdbtk_auto_download_db` | `true` | Whether to automatically download the GTDB-Tk database if it is missing. |
| `gtdbtk_extension` | `fa` | File extension used for MAG files passed to GTDB-Tk. |
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
| `publish_mags_mode` | `copy` | How MAG files are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `publish_tool_outputs_mode` | `copy` | How annotation and quality-control tool outputs are published to the output directory. Options can be `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results directory. Uses `working_dir` if provided, otherwise `output_dir`, otherwise `.`. Usually does not need to be set directly. |
| `outdir` | `<results_dir>/module_6_magannotation` | Module 6 output directory. Usually derived automatically and does not need to be set directly. |
| `module5_final_mag_dir` | `<results_dir>/module_5_subtractiveassembly/final_mag_database` | Candidate MAG directory from Module 5 subtractive assembly. |
| `module4_refined_mag_dir` | `<results_dir>/module_4_binrefinement/refined_bins` | Candidate refined MAG directory from Module 4 bin refinement. |
| `checkm2_db_outdir` | `<outdir>/databases/checkm2` | Derived CheckM2 database output directory. Uses `checkm2_db_dir` if provided. Usually does not need to be set directly. |
| `gtdbtk_db_outdir` | `<outdir>/databases/gtdbtk` | Derived GTDB-Tk database output directory. Uses `gtdbtk_db_dir` if provided. Usually does not need to be set directly. |
| `eggnog_db_outdir` | `<outdir>/databases/eggnog` | Derived EggNOG-mapper database output directory. Uses `eggnog_data_path` or `eggnog_data_dir` if provided. Usually does not need to be set directly. |
