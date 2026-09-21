# Module 1: Read trimming

![Module 1 workflow](../images/step_1.png)

Module 1 trims reads that have been discovered and validated by [Module 0: Read
processing](module-0-read-processing.md). It runs `fastp` on paired-end and
interleaved reads, runs FastQC on the trimmed reads, and writes per-sample and
combined trimming statistics for downstream modules.

## Prerequisite

Run Module 0 first. By default, Module 1 expects the Module 0 manifest at:

```text
<working_dir>/module_0_readprocess/naming/read_manifest.tsv
```

If `--working_dir` is omitted, the resolved SAMWISE source directory is used as
the results root. Use `--input_manifest` to provide a different manifest:

```bash
nextflow run module_1_readtrimming.nf \
  --working_dir ./output_samwise \
  --input_manifest ./somewhere/read_manifest.tsv
```

The manifest must be tab-separated and include the fields used by Module 1:
`sample_id`, `layout`, and the corresponding `read1`/`read2` fields for
`paired` rows or `interleaved` for `interleaved` rows. The referenced input
files must exist. Module 0 is responsible for read naming, pairing/layout
detection, and optional FASTQ structure validation.

## Usage

```bash
nextflow run module_1_readtrimming.nf \
  --working_dir ./output_samwise \
  --threads 6
```

The default manifest discovery makes additional input flags unnecessary when
running after Module 0 in the same results directory.

## Parameters

| Parameter | Default | Description |
|---|---:|---|
| `working_dir` | `null` | Results root. Outputs are written to `<working_dir>/module_1_readtrimming`; when omitted, it resolves to `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing workflows, `bin/`, and `dependencies/`; when omitted, it resolves to the directory containing the launched workflow. |
| `input_manifest` | `null` | Override for the input read manifest. Otherwise Module 1 uses `<results_dir>/module_0_readprocess/naming/read_manifest.tsv`. |
| `threads` | `null` | Global CPU override for both `fastp` and FastQC. When set, it takes precedence over `fastp_threads` and `fastqc_threads`. |
| `fastp_threads` | `4` | CPUs per `fastp` trimming task when `threads` is not set. |
| `fastqc_threads` | `2` | CPUs per FastQC task when `threads` is not set. |
| `fastp_version` | `0.23.4` | `fastp` version requested during managed-tool installation. |
| `fastqc_version` | `0.12.1` | FastQC version requested during managed-tool installation. |
| `auto_install` | `true` | If required tools are unavailable, install them with `mamba` or `conda`. If `false`, both tools must already be available. |
| `tool_env_dir` | `null` | Optional conda environment path for Module 1 tools. Otherwise the environment is created at `<outdir>/conda_envs/module1_tools`. |
| `skip_fastqc` | `false` | Skip FastQC on trimmed reads when `true`. Tool setup still checks/installs both `fastp` and FastQC. |
| `publish_trimmed_mode` | `symlink` | Publication mode for trimmed FASTQ files. Supported values are `symlink`, `copy`, and `move`. |
| `compression` | `4` | Compression level for trimmed FASTQ output. |
| `detect_adapter_for_pe` | `true` | Enable `fastp` adapter detection for paired-end reads. |
| `enable_correction` | `false` | Enable `fastp` base correction for paired-end data. |
| `cut_front` | `true` | Enable quality trimming from the front of reads. |
| `cut_tail` | `true` | Enable quality trimming from the tail of reads. |
| `cut_window_size` | `4` | Sliding-window size for quality trimming. |
| `cut_mean_quality` | `30` | Minimum mean quality in the trimming window. |
| `qualified_quality_phred` | `30` | Phred threshold for a base to be considered qualified. |
| `unqualified_percent` | `40` | Maximum percentage of unqualified bases allowed in a read. |
| `n_base_limit` | `5` | Maximum number of `N` bases allowed before filtering. |
| `length_required` | `75` | Minimum read length after trimming/filtering; shorter reads are discarded. |
| `trim_poly_g` | `false` | Enable poly-G tail trimming. |
| `trim_poly_x` | `false` | Enable poly-X tail trimming. |
| `results_dir` | `working_dir` | Internal results root derived from the resolved `working_dir`. |
| `module0_outdir` | `<results_dir>/module_0_readprocess` | Derived Module 0 output directory used for default manifest discovery. |
| `outdir` | `<results_dir>/module_1_readtrimming` | Module 1 output directory. |

## Tool setup and environment

The setup process writes `setup/module1_tools_status.env` and selects tools in
this order:

1. Reuse `tool_env_dir`, or the default managed environment, if it exists and
   both executables pass `--version` checks.
2. Otherwise use working `fastp` and `fastqc` executables already on `PATH`.
3. Otherwise, when `auto_install=true`, create the managed environment with
   `mamba` if available, or `conda` if available, using the `conda-forge` and
   `bioconda` channels.

An existing incomplete or broken managed environment is removed before
reinstallation. With `auto_install=false`, the workflow fails if either tool
is missing or unusable. Trimming and FastQC tasks read the selected tool paths
from the status file and add the managed environment to `PATH` when needed.

## Output layout

For `<outdir> = <results_dir>/module_1_readtrimming`, published outputs are:

```text
module_1_readtrimming/
├── conda_envs/
│   └── module1_tools/                 # default managed environment, if created
├── setup/
│   └── module1_tools_status.env
├── trimmed_reads/
│   ├── <safe_id>_R1_trimmed.fastq.gz
│   ├── <safe_id>_R2_trimmed.fastq.gz
│   └── <safe_id>_interleaved_trimmed.fastq.gz
├── fastp_reports/
│   ├── <safe_id>_fastp.html
│   ├── <safe_id>_fastp.json
│   └── <safe_id>_fastp.log
├── fastqc_reports/
│   ├── <read>_fastqc.html
│   └── <read>_fastqc.zip
└── summary/
    ├── <safe_id>_trimming_stats.tsv
    ├── trimming_stats_summary.tsv
    └── trimmed_manifest.tsv
```

`trimmed_reads` uses `publish_trimmed_mode`; setup, FastQC, fastp reports, and
summary files are published with copy mode. FastQC files are absent when
`skip_fastqc=true`.

`summary/trimmed_manifest.tsv` has these columns:

```text
sample_id  safe_sample_id  layout  read1  read2  interleaved  merged  fastp_html  fastp_json
```

The read and report paths in this manifest are absolute published paths so the
manifest can be consumed from later Nextflow task directories. Paired rows
populate `read1` and `read2`; interleaved rows populate `interleaved`.

Each per-sample stats file and `trimming_stats_summary.tsv` report the fastp
before/after counts, reads removed by filters, adapter-trimmed reads, final
FASTQ record counts, and the associated fastp JSON path. Paired final counts
are measured separately for read 1 and read 2; interleaved counts are measured
from the single interleaved output.

## Validation notes

- Module 1 checks that the resolved manifest exists before processing it and
  fails when no trimmed-manifest records or trimming-stat files are produced.
- The manifest layout value must be exactly `paired` or `interleaved`; each
  layout requires its corresponding path fields.
- Module 0 should be used to validate FASTQ structure. Module 1 additionally
  counts output FASTQ records and fails if a FASTQ line count is not divisible
  by four.
- Sample identifiers are converted to safe output names by replacing any
  character outside `[A-Za-z0-9._-]` with `_`. Different identifiers can
  therefore map to the same safe name; use unique, filesystem-safe sample IDs.
- Check `setup/module1_tools_status.env` and the FastQC/fastp reports when
  diagnosing tool or trimming failures.

## Next step: Module 2 assembly

After a successful run, Module 2 and Module 2B consume
`summary/trimmed_manifest.tsv` from this output directory. The manifest is the
handoff between read trimming and assembly; it contains the trimmed read paths,
layout, sample IDs, and associated report paths.

Continue with [Module 2: Assembly](module-2-assembly.md), or use [Module 2B:
Grouped coassembly](module-2b-coassembly.md) when assigning samples to
coassembly groups. Both workflows accept `--input_manifest` when the default
trimmed manifest is not being used.
