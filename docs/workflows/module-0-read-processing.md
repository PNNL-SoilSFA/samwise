# Module 0: Read processing

![Module 0 workflow](../images/step_0.png)

Module 0 performs initial read preprocessing and validation. It checks read
filenames, detects paired-end or interleaved layouts, optionally validates
FASTQ structure, and runs FastQC on every read that passes the applicable
checks.

## Usage

```bash
nextflow run module_0_readprocess.nf \
  --input_dir ./reads_dir \
  --threads 6 \
  --working_dir ./output_samwise
```

`--input_dir` is required. For separate source and results directories, run
the workflow by its absolute path:

```bash
nextflow run /path/to/samwise/module_0_readprocess.nf \
  --samwise_dir /path/to/samwise \
  --working_dir /path/to/samwise-results \
  --input_dir /path/to/reads_dir
```

`samwise_dir` is the SAMWISE installation/source directory containing the
workflows, `bin/`, `background/`, `dependencies/`, and related files.
`working_dir` is the results and environment directory. The fallback rules
are:

- Supplying both parameters keeps source code under `samwise_dir` and writes
  results, generated Conda environments, and downloaded tool databases below
  `working_dir`.
- Supplying only `working_dir` uses the directory containing the launched
  workflow as `samwise_dir`.
- Supplying only `samwise_dir` uses that source directory as `working_dir`.

## Arguments

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 0 outputs are written to `<working_dir>/module_0_readprocess`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it defaults to the directory containing the launched workflow. |
| `input_dir` | `null` | Required directory where the metagenomic reads are stored. |
| `threads` | `null` | Total number of threads to use. If provided, this overrides `fastqc_threads`. |
| `skip_validate` | `false` | Skip FASTQ structure validation. Naming and pairing checks still run; every named read is passed to FastQC. |
| `ignore_invalid_fastq` | `false` | Continue after FASTQ structure validation failures. Each invalid read is excluded from the Module 0 manifest, its validation report includes `READ IGNORED - BAD FASTQ FILE`, and `validation_reports/ignored_bad_fastq_files.tsv` records the exclusion. The default (`false`) stops the workflow on the first invalid FASTQ. Has no effect with `--skip_validate true`. |
| `file_pattern` | `*.{fastq.gz,fq.gz,fastq,fq}` | Glob/text pattern used to detect input files inside `input_dir`. The default detects `fastq.gz`, `fq.gz`, `fastq`, and `fq` files. We recommend leaving this unchanged. |
| `fastqc_threads` | `2` | Number of threads to use specifically for FastQC if the global `threads` argument is not passed. |
| `fastqc_version` | `0.12.1` | FastQC version to install if a different version is desired. |
| `auto_install` | `true` | Controls whether SAMWISE installs required packages, such as FastQC. If set to `false`, FastQC must already be available in your environment. |
| `tool_env_dir` | `null` | Optional custom directory for the Module 0 tool environment. If unset, the module uses `<module_0_outdir>/conda_envs/module0_tools`. |
| `outdir` | `<working_dir>/module_0_readprocess` | Derived Module 0 output directory. |

## Input Reads

Reads must use one of the supported paired-end or interleaved naming formats,
with either compressed or uncompressed FASTQ extensions:

- Paired-end: `SampleA_R1.fastq.gz` and `SampleA_R2.fastq.gz`
- Paired-end: `SampleB_1.fq` and `SampleB_2.fq`
- Interleaved: `SampleC_interleaved.fastq`
- Illumina-style paired-end names are also supported, such as
  `sample_S1_L001_R1_001.fastq.gz` and `sample_S1_L001_R2_001.fastq.gz`.

Valid extensions are `.fq`, `.fastq`, `.fq.gz`, and `.fastq.gz`. R1/R2 files
must use the same naming style, and each R1 must have a matching R2. A sample
cannot have both paired-end and interleaved files. Files with FASTQ-like
extensions that do not match a supported naming convention, non-FASTQ files,
unpaired reads, duplicate read files, and mixed R1/R2 styles fail naming
validation.

For example, this is a valid input tree:

```text
reads_dir/
├── SampleA_R1.fastq.gz
├── SampleA_R2.fastq.gz
├── SampleB_1.fq
├── SampleB_2.fq
└── SampleC_interleaved.fastq
```

The workflow uses the complete filename prefix before the recognized read
suffix as the sample ID: `SampleA`, `SampleB`, and `SampleC` in the example
above. This is the implemented prefix-based behavior; punctuation, including
periods, is not independently removed or used as a truncation point. For
`sample_S1_L001_R1_001.fastq.gz`, the sample ID is `sample_S1_L001`.
Choose one unambiguous sample identifier before `_R1`, `_R2`, `_1`, `_2`, or
`_interleaved`.

Do not create names whose distinct reads normalize to the same prefix and read
role. For example, `Sample_R1.fastq` and `Sample_R1_001.fastq` both use sample
ID `Sample` and read role R1, so they are reported as duplicate R1 files. A
paired sample combined with an interleaved file using the same prefix also
fails. Duplicate normalized sample IDs or duplicate read roles cause Module 0
to fail naming validation. Use distinct prefixes and the supported suffixes
instead.

## FASTQ Validation

Unless `--skip_validate true` is set, each named read is checked for:

- A supported extension and readable compressed or uncompressed content.
- FASTQ records in groups of four lines.
- Header lines beginning with `@`.
- Non-empty sequence, plus, and quality lines.
- Matching sequence and quality lengths.
- A non-empty file.

With the default `--ignore_invalid_fastq false`, an invalid FASTQ fails the
workflow. With `--ignore_invalid_fastq true`, invalid reads are omitted from
the Module 0 read manifest and are not sent to FastQC. The workflow still
fails if no valid samples remain. With `--skip_validate true`, structure
validation is skipped, but filename, extension, naming, and pairing checks
still run and the named reads are passed directly to FastQC.

## FastQC

Module 0 checks the configured `tool_env_dir` first, then a working `fastqc`
executable on `PATH`. If neither is available and `--auto_install true` is
set, SAMWISE creates the environment at `tool_env_dir`, or at
`<module_0_outdir>/conda_envs/module0_tools` when it is unset. It uses `mamba`
when available and otherwise `conda`, and installs the requested
`fastqc_version` plus Perl. Existing but incomplete or failing module-local
environments are recreated. With `--auto_install false`, setup fails unless a
working FastQC executable is already available.

FastQC uses `threads` when it is provided; otherwise it uses
`fastqc_threads` (default `2`). FastQC produces an HTML report and ZIP report
for every retained read.

## Outputs

Module 0 publishes outputs below
`<working_dir>/module_0_readprocess`:

```text
module_0_readprocess/
├── naming/
│   ├── read_manifest.tsv
│   └── read_naming_report.txt
├── validation_reports/
│   ├── *_validation.txt
│   ├── *_validation_skipped.txt
│   └── ignored_bad_fastq_files.tsv
├── fastqc_reports/
│   ├── *_fastqc.html
│   └── *_fastqc.zip
└── setup/
    └── module0_tools_status.env
```

The `naming/read_manifest.tsv` file is tab-separated with the stable header
`sample_id`, `layout`, `read1`, `read2`, `interleaved`. Each row contains either
the absolute R1/R2 paths for a `paired` sample or the absolute interleaved path
for an `interleaved` sample; unused path columns are empty. This is the
filtered manifest consumed by FastQC and downstream modules. The validation
reports document per-read structure checks or skipped checks. The ignored-read
table has the header `sample_id`, `read_file`, `reason`, `message` and records
excluded invalid reads. Nextflow work directories may also contain
intermediate files; the published directories above are the Module 0 outputs
intended for downstream modules.

## Validation Notes

Module 0 must complete naming validation before reads are made available to
downstream steps. FastQC runs only on reads retained in the filtered manifest.
Review `read_naming_report.txt`, the validation reports, and
`ignored_bad_fastq_files.tsv` when a read is missing from the manifest or when
`--ignore_invalid_fastq true` was used.
