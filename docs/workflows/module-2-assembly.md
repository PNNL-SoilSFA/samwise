# Module 2: Assembly

![Module 2 workflow](../images/step_2.png)

Module 2 assembles the reads produced by Module 1. It supports independent
single-sample and rarefied assemblies with MEGAHIT, metaSPAdes, or both.
Both Module 2 and the separate [Module 2B grouped co-assembly
workflow](module-2b-coassembly.md) write manifests that Module 3 can consume.

The repository `README.md` documents the main commands and parameters only. The
input matching rules, validation behavior, generated group files, and manifest
field values described below are implementation behavior documented here; they
are not all described in the README.

## Inputs and prerequisites

Complete [Module 1: Read trimming](module-1-read-trimming.md) first. Module 2
and Module 2B consume its published trimmed-read manifest by default:

```text
<results_dir>/module_1_readtrimming/summary/trimmed_manifest.tsv
```

Use `--input_manifest` to provide another trimmed-read manifest. The complete
Module 1 handoff, output layout, and validation guidance are documented in the
[Module 1 next-step note](module-1-read-trimming.md#next-step-module-2-assembly).
The manifest must be tab-separated with one row per trimmed sample:

```text
sample_id  safe_sample_id  layout  read1  read2  interleaved  merged  fastp_html  fastp_json
```

The fields used by assembly are:

| Column | Required content |
|---|---|
| `sample_id` | Original sample identifier. Module 2 derives the assembly identifier by removing non-alphanumeric characters. |
| `safe_sample_id` | Filesystem-safe sample identifier used for output filenames and fallback identifier derivation. |
| `layout` | Exactly `paired` or `interleaved`. |
| `read1`, `read2` | Paths to the two reads when `layout=paired`. |
| `interleaved` | Path to the interleaved FASTQ when `layout=interleaved`. |

`merged`, `fastp_html`, and `fastp_json` are retained by the Module 1
contract but are not used to run assembly. Relative read paths are resolved
relative to the Nextflow launch directory; absolute paths are used as given.

## Usage

Select at least one assembler and at least one assembly mode. A standard run
with both assemblers is:

```bash
nextflow run module_2_readassembly.nf \
  --working_dir ./output_samwise \
  --threads 6 \
  --memory_gb 0 \
  --megahit \
  --metaspades
```

Enable rarefaction explicitly. This example creates two subsets per sample
for each selected assembler in addition to the single assemblies:

```bash
nextflow run module_2_readassembly.nf \
  --working_dir ./output_samwise \
  --threads 6 \
  --memory_gb 0 \
  --megahit \
  --metaspades \
  --rarefied_assembly true \
  --rarefaction_splits 2
```

`--single_assembly` defaults to `true`; set it to `false` to run only
rarefied assemblies. `--rarefied_assembly` defaults to `false`. The workflow
validates `--rarefaction_splits` on every run, so its value must be at least
`2` even when rarefied assembly is disabled. When rarefaction is enabled, split
labels are generated alphabetically:
`a` through `z`, then `aa`, `ab`, and so on. The label is appended to the
assembly sample identifier and output basename.

Assembly strategy codes in `summary/assembly_manifest.tsv` are:

| Code | Assembly |
|---|---|
| `A` | Single MEGAHIT |
| `B` | Single metaSPAdes |
| `C` | Rarefied MEGAHIT |
| `D` | Rarefied metaSPAdes |

If MEGAHIT is the selected assembler and more than one thread is not
supported by the host, use `--megahit_threads 1`; this overrides the global
thread setting for MEGAHIT only.

## Parameters

| Parameter | Default | Description |
|---|---:|---|
| `working_dir` | `null` | Results root. Outputs go to `<working_dir>/module_2_readassembly`. If omitted, uses `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source tree containing the workflow, `bin/`, and dependencies. If omitted, uses the launched workflow directory. |
| `input_manifest` | `null` | Input trimmed-read manifest. Defaults to Module 1's `summary/trimmed_manifest.tsv`. |
| `megahit` | `false` | Enable MEGAHIT. At least one of `megahit` or `metaspades` is required. |
| `metaspades` | `false` | Enable metaSPAdes. |
| `single_assembly` | `true` | Run one assembly using all input reads. |
| `rarefied_assembly` | `false` | Run assemblies on each rarefied read subset. |
| `rarefaction_splits` | `2` | Number of subsets. Must be at least `2` on every run, including when rarefaction is disabled. |
| `megahit_version` | `1.2.9` | MEGAHIT version for automatic installation. |
| `spades_version` | `4.2.0` | SPAdes/metaSPAdes version for automatic installation. |
| `auto_install` | `true` | Install missing tools with `mamba` or `conda`; set `false` to require tools already in the environment. |
| `tool_env_dir` | `null` | Custom conda environment path. Otherwise uses `<outdir>/conda_envs/module2_tools`. |
| `threads` | `null` | Global task thread override. |
| `assembly_threads` | `4` | Task threads when `threads` is unset. |
| `megahit_threads` | `null` | MEGAHIT-specific thread override. |
| `memory_gb` | `0` | Assembly memory limit. A positive value is passed as a memory limit; `0` lets the assembler auto-detect available memory with safety headroom. |
| `clean_partial_assembler_outputs` | `true` | Remove stale partial assembler outputs before each task. |
| `megahit_preset` | `meta-large` | MEGAHIT preset. `meta-sensitive` is another supported choice. |
| `publish_assemblies_mode` | `symlink` | Publication mode for final FASTA files: `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results root derived from `working_dir`; normally do not set directly. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. |
| `outdir` | `<results_dir>/module_2_readassembly` | Module 2 output directory. |

## Outputs

The output directory contains:

```text
module_2_readassembly/
|-- assemblies/                 # *.renamed.fa, using publish_assemblies_mode
|-- header_maps/                # old-to-new contig headers
|-- logs/                       # per-assembly logs
|-- summary/
|   |-- assembly_manifest.tsv
|   |-- assembly_stats_summary.tsv
|   `-- per_assembly_stats/
|-- setup/module2_tools_status.env
`-- conda_envs/module2_tools/   # when the workflow creates the environment
```

`summary/assembly_manifest.tsv` has this complete schema:

```text
sample_id  safe_sample_id  assembly_sample_id  assembler  assembly_mode  rarefaction_label  assembly_strategy  renamed_fasta
```

`assembly_mode` is `single` or `rarefied`; `rarefaction_label` is empty for a
single assembly. `renamed_fasta` is the published FASTA path. Per-assembly
statistics and `assembly_stats_summary.tsv` add these columns:

```text
sample_id  safe_sample_id  assembly_sample_id  assembler  assembly_mode  rarefaction_label  assembly_strategy  assembly_status  assembly_warning  contigs  total_bp  max_contig_bp  n50_bp  renamed_fasta
```

Contigs are renamed for downstream use and header maps record the original and
new headers. The default `symlink` mode keeps final FASTA files in Nextflow
work directories; use `copy` for accessible independent copies. `move` is
supported but can interfere with Nextflow's ability to find staged outputs.

## SLURM execution

Without a config file, assembly tasks use the local executor. To submit each
single or rarefied assembly as its own SLURM job:

```bash
nextflow run module_2_readassembly.nf \
  -c ./bin/module_2_slurm.config \
  --working_dir ./output_samwise \
  --slurm_account ChargeAccountID \
  --max_parallel_assemblies 8 \
  --threads 36 \
  --memory_gb 0 \
  --megahit \
  --metaspades \
  --rarefied_assembly true \
  --rarefaction_splits 2
```

The profile defaults are `slurm_account=null`, `slurm_partition=null`,
`slurm_qos=null`, `slurm_extra_options=''`, `max_parallel_assemblies=8`,
`assembly_time='48h'`, `assembly_max_retries=2`, `slurm_exclusive=true`,
`slurm_request_all_memory=true`, and `slurm_before_script=''`. It submits
assembly tasks with one exclusive node, `--mem=0` by default, and caps both
the process forks and global executor queue. Setup and summary tasks remain
local. `slurm_extra_options` and `slurm_before_script` provide site-specific
customization.

Keep `--memory_gb 0` with this profile. The profile requests the whole node
and lets the assembler detect node memory; a nonzero value would create a
conflicting second SLURM memory request. Set `max_parallel_assemblies` to the
number of simultaneous nodes allowed by the allocation. Reports are written
under `module_2_readassembly/pipeline_info/`.
