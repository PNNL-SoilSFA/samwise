# Requirements

SAMWISE requires:

- [Nextflow](https://www.nextflow.io/)
- `mamba` (or `conda`) for automatic package installation
- Java, which Nextflow requires and is installed automatically when installing
  with `mamba` or `conda`

Install `mamba` (or `conda`) first. Mamba is recommended; installation
instructions are available at <https://conda-forge.org/download/>.

Install Nextflow from the Bioconda package. We recommend using a separate
environment rather than the base environment:

```bash
mamba install -n nextflow -c bioconda nextflow
mamba activate nextflow
```

Nextflow can also be installed through the
[Bioconda Nextflow package](https://anaconda.org/channels/bioconda/packages/nextflow/overview).

Clone SAMWISE after installing the prerequisites, or download and extract the
repository in the location of your choice:

```bash
git clone https://github.com/PNNL-SoilSFA/samwise.git
```

## Source and results directories

Use `--samwise_dir` for the SAMWISE source directory and `--working_dir` for
workflow results, generated environments, and downloaded tool databases:

```bash
nextflow run /path/to/samwise/module_0_readprocess.nf \
  --samwise_dir /path/to/samwise \
  --working_dir /path/to/samwise-results
```

The downloaded `samwise` directory contains the Nextflow workflows, helper
scripts, and bundled dependencies. The directory options have these fallback
rules:

- Supplying both parameters keeps source code under `--samwise_dir` and
  results, created environments, and downloaded databases under
  `--working_dir`.
- Supplying only `--working_dir` uses the directory containing the launched
  workflow as `samwise_dir`; results and code are under `--working_dir`.
- Supplying only `--samwise_dir` uses that source directory as the default
  `working_dir`; results and code are under `--samwise_dir`.

For separate source and results directories, invoke the workflow by its
absolute path as shown above.
