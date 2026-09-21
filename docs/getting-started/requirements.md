# Requirements

SAMWISE requires:

- [Nextflow](https://www.nextflow.io/)
- `mamba` or `conda` for workflow-managed environments
- Java, which Nextflow requires

We recommend installing Nextflow into a separate environment:

```bash
mamba create -n nextflow -c bioconda nextflow
mamba activate nextflow
```

Clone the repository after installing the prerequisites:

```bash
cd samwise
```

## Source and results directories

Use `--samwise_dir` for the SAMWISE source directory and `--working_dir` for
workflow results, generated environments, and downloaded tool databases:

```bash
nextflow run /path/to/samwise/module_0_readprocess.nf \
  --samwise_dir /path/to/samwise \
  --working_dir /path/to/samwise-results
```

Keeping source and results separate is recommended for repeatable runs and
cleaner repository checkouts.
