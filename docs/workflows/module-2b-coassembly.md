# Module 2B: Grouped coassembly

![Module 2B workflow](../images/step_2b.png)

Module 2B co-assembles selected Module 1 samples into groups using MEGAHIT
only. Each group becomes one coassembly and uses strategy code `G`. It is
intended to run alongside Module 2 and produces both assembly and trimmed-read
manifests for Module 3. metaSPAdes is not supported for coassembly because of
its memory requirements.

## Coassembly group manifest

Pass a required tab-separated file with two columns:

```text
read_or_sample_id  group_id
sample1             group_1
sample2             group_1
sample3             group_2
```

The repository example is `examples/coassembly_manifest_example.txt`. A
header is optional. If present, the first row is skipped only when its first
field is one of `read`, `read_id`, `read_or_sample_id`, `fastq`, `fastq_id`,
`fastq_file`, `fastq_file_id`, `sample`, or `sample_id`, and its second field
is one of `group`, `group_id`, or `coassembly_group`.

For each non-empty data row:

- Both columns must be non-empty and separated by a tab.
- Column 1 may match a Module 1 `sample_id`, `safe_sample_id`, read path,
  read basename, or basename with `.fastq.gz`, `.fq.gz`, `.fastq`, or `.fq`
  removed. Relative paths are resolved from the Nextflow launch directory.
- A sample may occur only once, or repeated assignments must use the same
  group. Conflicting assignments are an error.
- Unmatched rows are logged as warnings. If no rows match, the workflow stops.
- Only samples listed in the group file are included; other Module 1 samples
  are ignored.

The generated per-group file is `<safe_group_id>.coassembly_reads.tsv` with
this complete schema:

```text
group_id  safe_group_id  sample_id  safe_sample_id  layout  read1  read2  interleaved
```

Group IDs are sanitized by replacing every character other than alphanumeric,
period, underscore, and hyphen with `_`, then removing leading and trailing
underscores. An empty result becomes `unnamed_group`. If two groups produce the
same sanitized ID, the later group in sorted group-ID order gets `_2`, `_3`,
and so on. At least one valid group must be created. The generated group TSV
preserves both the original `group_id` and the sanitized `safe_group_id`.

Before assembly, Module 2B validates the selected Module 1 rows. `layout`
must be `paired` or `interleaved`. Paired rows require both files, both files
must exist, and their FASTQ record counts must match. Interleaved rows require
an existing file with an even number of FASTQ records. FASTQ records must be
complete. Paired reads are written R1 then R2; interleaved reads are copied as
provided. The resulting group FASTQ is a compressed interleaved file.

## Usage

```bash
nextflow run module_2b_coassembly.nf \
  --working_dir ./output_samwise \
  --coassembly_groups ./coassembly_manifest.txt \
  --threads 6 \
  --memory_gb 0
```

## Parameters

| Parameter | Default | Description |
|---|---:|---|
| `working_dir` | `null` | Results root. Outputs go to `<working_dir>/module_2b_coassembly`. If omitted, uses `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source tree. If omitted, uses the launched workflow directory. |
| `input_manifest` | `null` | Module 1 trimmed-read manifest; defaults to `summary/trimmed_manifest.tsv`. |
| `coassembly_groups` | `null` | Required two-column tab-separated group manifest. |
| `megahit_version` | `1.2.9` | MEGAHIT version for automatic installation. |
| `auto_install` | `true` | Install MEGAHIT with `mamba` or `conda` when unavailable. |
| `tool_env_dir` | `null` | Custom conda environment path. Otherwise uses `<outdir>/conda_envs/module2b_tools`. |
| `threads` | `null` | Global task thread override. |
| `assembly_threads` | `4` | Task threads when `threads` is unset. |
| `megahit_threads` | `null` | MEGAHIT-specific thread override. |
| `memory_gb` | `0` | Assembly memory limit; `0` lets MEGAHIT use auto-detected task memory. |
| `megahit_preset` | `meta-large` | MEGAHIT preset. |
| `publish_assemblies_mode` | `symlink` | Publication mode for coassembly FASTAs: `symlink`, `copy`, or `move`. |
| `publish_coassembly_reads_mode` | `symlink` | Publication mode for generated interleaved group FASTQs. |
| `results_dir` | `null` | Internal results root derived from `working_dir`. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. |
| `outdir` | `<results_dir>/module_2b_coassembly` | Module 2B output directory. |

## Outputs and manifests

The output directory contains:

```text
module_2b_coassembly/
|-- assemblies/                         # *.renamed.fa
|-- header_maps/
|-- logs/
|-- coassembly/
|   |-- groups/                         # prepared group TSVs
|   `-- interleaved_inputs/             # generated group FASTQs
|-- summary/
|   |-- assembly_manifest.tsv
|   |-- assembly_stats_summary.tsv
|   |-- coassembly_group_summary.tsv
|   |-- coassembly_trimmed_manifest.tsv
|   `-- per_assembly_stats/
|-- setup/module2b_tools_status.env
`-- conda_envs/module2b_tools/          # when created by the workflow
```

The Module 2B `assembly_manifest.tsv` has the same schema as Module 2 and one
row per group. Its values are `assembler=megahit`, `assembly_mode=coassembly`,
an empty `rarefaction_label`, and `assembly_strategy=G`. In this generated
manifest, `sample_id` and `safe_sample_id` are both the sanitized group ID;
`assembly_sample_id` is that value with non-alphanumeric characters removed
(falling back to `coassembly` only if that is empty). This reflects the
implementation's assembly job channel, which derives both IDs from the
generated `<safe_group_id>.coassembly_reads.tsv` filename rather than carrying
the original group ID into the assembly process.

`coassembly_group_summary.tsv` has:

```text
group_id  safe_group_id  sample_count  paired_sample_count  interleaved_sample_count  group_reads_tsv
```

`coassembly_trimmed_manifest.tsv` is the Module-3-compatible read manifest:

```text
sample_id  safe_sample_id  layout  read1  read2  interleaved  merged  fastp_html  fastp_json
```

Each coassembly row uses the sanitized group ID as both `sample_id` and
`safe_sample_id`, and the published generated interleaved FASTQ in
`interleaved`. The other read/report fields are empty. The original group ID is
retained in the prepared group TSV and `coassembly_group_summary.tsv`, but is
not emitted in these assembly-stage manifest records. This FASTQ is the exact
concatenated read input used to create that group's assembly, so Module 3 can
map/bin the coassembly with the matching reads.

The coassembly FASTA, header map, logs, per-assembly statistics, and manifest
use the same publication and statistics conventions as Module 2. The default
`symlink` mode keeps large files in work directories; use `copy` when durable,
directly browsable files are required. `move` is supported but may interfere
with Nextflow staging.

## SLURM execution

Use the separate profile to run one SLURM job per group:

```bash
nextflow run module_2b_coassembly.nf \
  -c ./bin/module_2b_slurm.config \
  --working_dir ./output_samwise \
  --coassembly_groups ./coassembly_manifest.txt \
  --slurm_account ChargeAccountID \
  --max_parallel_coassemblies 4 \
  --threads 36 \
  --memory_gb 0
```

The profile defaults are `slurm_account=null`, `slurm_partition=null`,
`slurm_qos=null`, `slurm_extra_options=''`, `max_parallel_coassemblies=4`,
`coassembly_time='48h'`, `coassembly_max_retries=2`,
`slurm_exclusive=true`, `slurm_request_all_memory=true`, and
`slurm_before_script=''`. Setup, group preparation, and summary tasks remain
local; each `ASSEMBLE_COASSEMBLY` task uses one exclusive node and the global
queue is capped at the configured group count.

Keep `--memory_gb 0` with this profile so the exclusive-node `--mem=0`
request and MEGAHIT memory auto-detection remain consistent. Coassembly
parallelism is bounded by the number of groups: with one group there is no
parallel speed-up, and the largest group determines a lower bound on runtime.
Reports are written under `module_2b_coassembly/pipeline_info/`.

## Operational notes

- Module 2B supports MEGAHIT only; metaSPAdes coassembly remains unavailable.
- Coassembly group rows that do not match the trimmed manifest are warnings,
  not errors, unless every row is unmatched. Check the preparation log when
  validating a group file.
- The SLURM profiles assume shared filesystem access and site-compatible
  `--exclusive` and `--mem=0` options. Cluster-specific settings may require
  `slurm_partition`, `slurm_qos`, `slurm_extra_options`, or disabling
  `slurm_request_all_memory`.
