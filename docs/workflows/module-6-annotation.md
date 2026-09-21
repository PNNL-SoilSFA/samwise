# Module 6: MAG Annotation

![Module 6 workflow](../images/step_6.png)

`module_6_magannotate.nf` prepares normalized MAG FASTAs and runs the selected
quality, dereplication, taxonomy, functional annotation, and trait tools.

## Usage

```bash
nextflow run module_6_magannotate.nf \
  --working_dir ./output_samwise \
  --run_drep true \
  --run_checkm2 true \
  --run_gtdbtk true \
  --run_eggnog true \
  --run_microtrait true \
  --threads 32
```

Only selected tools are installed and run. `--auto_install true` uses mamba or
conda when an environment is absent. Tool environments and downloaded data are
kept below the Module 6 output directory unless an explicit path is supplied.

`dRep` requires CheckM2 because the workflow converts the CheckM2 quality report
to the `genome,completeness,contamination` CSV passed to dRep as `--genomeInfo`.

## MAG Inputs

Input selection is applied in this order:

1. `--input_mag_dir`, when supplied.
2. Module 5 `final_mag_database`, when it exists and contains MAG FASTAs.
3. Module 4 `refined_bins`, as the fallback.

When Module 5 is selected and no manifest is supplied, the workflow uses
`summary/final_mag_database_manifest.tsv`. An explicit `--input_mag_manifest`
always takes precedence over automatic manifest selection. A manifest may use
`final_mag_fasta`, `refined_bin_fasta`, `bin_fasta`, or `prepared_fasta`; IDs may
come from `final_mag_id`, `refined_bin_id`, `bin_id`, or `mag_id`.

When using a Module 3 `summary/binning_manifest.tsv`, also pass the matching
`--input_mag_dir`. That directory must contain the MAG FASTA basenames listed in
the manifest's `bin_fasta` column. The manifest does not locate or stage its
source directory; Module 6 stages the supplied directory and matches manifest
rows by basename. Do not combine a Module 3 manifest with a different MAG
directory, or its rows will be skipped as missing.

The preparation step accepts `.fa` and `.fa.gz` files, copies or decompresses
them into `refined_genomes/`, sanitizes MAG IDs to safe filename characters, and
adds numeric suffixes for duplicate names. All prepared files are uncompressed
`.fa`. Directory discovery therefore recognizes only `.fa` and `.fa.gz`; a
manifest-referenced file must likewise be available in the supplied directory
under its basename. `mag_extension`, `drep_extension`, `checkm2_extension`,
`gtdbtk_extension`, and `microtrait_extension` are implementation-constrained
to `fa` (a leading dot is accepted and stripped). These values cannot be used
to select another FASTA suffix. Missing manifest entries are recorded and
skipped; there must still be at least one usable MAG.

## Tools and Databases

The workflow runs tools in this order when enabled:

1. CheckM2 predicts MAG quality.
2. dRep dereplicates the prepared MAGs using the CheckM2-derived genome info.
3. GTDB-Tk classifies the prepared or dereplicated MAGs.
4. EggNOG-mapper creates a unified, MAG-prefixed contig FASTA and functional annotations.
5. microTrait runs per-MAG genomic trait prediction and merges its results.

### CheckM2

An explicit `--checkm2_db_path` is used directly. Otherwise the workflow looks
for `uniref100.KO.1.dmnd` below `--checkm2_db_dir` (default
`<outdir>/databases/checkm2`). If it is absent and auto-download is enabled,
the workflow first tries `checkm2 database --download --path <dir>` and then
falls back to the configured Zenodo record (`14897628`).

### GTDB-Tk

An explicit `--gtdbtk_data_path` is used directly and must be a GTDB-Tk database
root containing `metadata/metadata.txt`, `metadata.txt`, or `VERSION`.
Otherwise data are stored in `--gtdbtk_db_dir` (default
`<outdir>/databases/gtdbtk`). Automatic download uses the configured GTDB-Tk
package URL, by default:

```text
https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
```

The archive is validated, extracted with its top-level directory stripped, and
the resulting database root is validated before classification.

GTDB-Tk writes to a task-local `gtdbtk_out` directory. Published files are
placed directly in `<outdir>/gtdbtk/` (not in an additional `gtdbtk_out/`
subdirectory). The GTDB-Tk database is separate, under
`<outdir>/databases/gtdbtk` or the explicitly supplied `gtdbtk_db_dir`; the
download archive is removed after successful extraction. At run start, cleanup
removes only the managed result directories `<outdir>/drep_out`,
`checkm2`, `gtdbtk`, `eggnog`, and `microtrait`. It does not remove databases,
`setup/`, `logs/`, or `summary/`.

**HPC warning:** pplacer and GTDB-Tk can overallocate memory. The workflow
defaults `--gtdbtk_pplacer_cpus` to `1`; increase it cautiously only when the
scheduler allocation has adequate memory.

### EggNOG-mapper

`--eggnog_data_path` takes precedence over `--eggnog_data_dir`; otherwise data
are stored in `<outdir>/databases/eggnog`. A usable diamond-mode database must
contain `eggnog.db` and at least one `.dmnd` file. If needed,
`download_eggnog_data.py --data_dir <dir> <eggnog_download_args>` downloads the
database. The default download argument is `-y`. The optional
`eggnog-mapper-fixurl` package is installed and run by default to repair
database download URLs. An explicit `--eggnog_mmseqs_db` is used when supplied;
otherwise the workflow checks `<data_dir>/mmseqs/mmseqs.db` and
`<data_dir>/mmseqs.db`.

### microTrait

microTrait setup clones `microtrait_github_repo` (default
`https://github.com/ukaraoz/microtrait`) and installs gRodon from
`grodon_github_repo` (default `jlw-ecoevo/gRodon`). Optional Git refs pin either
clone. With `microtrait_auto_download_db true`, HMM databases are deployed by
`microtrait::prep.hmmmodels()`; otherwise they must already exist in the R
library.

## Commands

The commands recorded in the tool logs are equivalent to:

```bash
checkm2 predict --threads <threads> --input <mag_dir> -x fa \
  --output-directory checkm2_out --database_path <checkm2_db>

dRep dereplicate drep_out --processors <drep_threads> -g <mag_files> \
  --genomeInfo drep_genomeInfo_used.csv -sa 0.99 -comp 50 -con 10

gtdbtk classify_wf --genome_dir <mag_dir> --out_dir gtdbtk_out \
  --extension fa --cpus <threads> --pplacer_cpus 1

emapper.py -m diamond --cpu <threads> -i module6_eggnog_input.fasta \
  --itype metagenome --genepred prodigal --trans_table 11 \
  --data_dir <eggnog_data_dir> --output samwise_eggnog \
  --output_dir eggnog_out --excel

Rscript local_microTrait_runner.R microtrait_inputs samwise_microtrait genomic
Rscript microTrait_merger.R microtrait_inputs microtrait_out
```

Additional dRep, GTDB-Tk, and EggNOG arguments are appended from their
corresponding `*_extra_args` parameters and must be valid shell-style argument
text. The workflow supplies the fixed input, output, extension, and thread
options shown above; extra arguments should not duplicate or contradict them.

## Parameters

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results/environment root. Defaults to `samwise_dir`; outputs go to `<working_dir>/module_6_magannotate`. |
| `samwise_dir` | `null` | SAMWISE source tree containing workflows, `bin/`, and `dependencies/`; defaults to the launched workflow directory. |
| `input_mag_dir` | `null` | Explicit directory of `.fa` or `.fa.gz` MAG FASTAs; otherwise Module 5 then Module 4 is selected. Required as the matching source directory when `input_mag_manifest` is a Module 3 manifest. |
| `input_mag_manifest` | `null` | Explicit MAG manifest; otherwise the Module 5 final-MAG manifest is used when applicable. With a Module 3 manifest, its `bin_fasta` basenames must exist in `input_mag_dir`. |
| `mag_extension` | `fa` | Prepared MAG extension; implementation requires `fa` (with or without a leading dot). Input files may be `.fa` or `.fa.gz`; prepared files are always `.fa`. |
| `threads` | `null` | Global thread count; tool defaults are CheckM2 8, GTDB-Tk 32, EggNOG 30, and microTrait 20 when unset. |
| `tool_env_dir` | `null` | Base directory for automatically managed tool environments. |
| `auto_install` | `true` | Install missing tools with mamba or conda. |
| `conda_pkgs_dir` | `null` | Optional conda/mamba package cache root. |
| `run_drep` | `true` | Run dRep; requires `run_checkm2 true`. |
| `run_checkm2` | `true` | Run CheckM2 quality assessment. |
| `run_gtdbtk` | `true` | Run GTDB-Tk taxonomy. |
| `run_eggnog` | `false` | Run EggNOG-mapper functional annotation. |
| `run_microtrait` | `false` | Run microTrait genomic trait prediction. |
| `drep_version` | `null` | dRep package version. |
| `drep_env_dir` | `null` | Custom dRep environment directory. |
| `drep_extension` | `fa` | dRep input extension; must be `fa`. |
| `drep_threads` | `null` | dRep threads; otherwise uses `threads`, then 16. |
| `drep_extra_args` | `-sa 0.99 -comp 50 -con 10` | Shell-style arguments appended to dRep; do not duplicate or contradict fixed workflow options. |
| `checkm2_version` | `null` | CheckM2 package version. |
| `checkm2_db_path` | `null` | Existing CheckM2 `uniref100.KO.1.dmnd` file. |
| `checkm2_db_dir` | `null` | CheckM2 database directory; otherwise `<outdir>/databases/checkm2`. |
| `checkm2_zenodo_record` | `14897628` | Zenodo record for CheckM2 fallback download. |
| `checkm2_auto_download_db` | `true` | Download CheckM2 data when missing. |
| `checkm2_extension` | `fa` | CheckM2 input extension; must be `fa`. |
| `checkm2_tmp_dir` | `null` | Shared filesystem location for CheckM2 temporary data. |
| `gtdbtk_version` | `2.7.2` | GTDB-Tk package version. |
| `gtdbtk_data_path` | `null` | Existing GTDB-Tk database root. |
| `gtdbtk_db_dir` | `null` | GTDB-Tk database directory; otherwise `<outdir>/databases/gtdbtk`. |
| `gtdbtk_auto_download_db` | `true` | Download GTDB-Tk data when missing. |
| `gtdbtk_extension` | `fa` | GTDB-Tk input extension; must be `fa`. |
| `gtdbtk_download_url` | GTDB release 232 URL | GTDB-Tk database package URL. |
| `gtdbtk_pplacer_cpus` | `1` | pplacer CPUs; keep low on HPC systems. |
| `gtdbtk_extra_args` | `""` | Shell-style arguments appended to GTDB-Tk; do not duplicate or contradict fixed workflow options. |
| `eggnog_env_dir` | `null` | Custom EggNOG-mapper environment directory. |
| `eggnog_mapper_version` | `2.1.13` | EggNOG-mapper package version. |
| `eggnog_data_path` | `null` | Existing EggNOG data directory; takes precedence over `eggnog_data_dir`. |
| `eggnog_data_dir` | `null` | EggNOG data directory; otherwise `<outdir>/databases/eggnog`. |
| `eggnog_auto_download_db` | `true` | Download EggNOG data when incomplete. |
| `eggnog_download_args` | `-y` | Arguments for `download_eggnog_data.py`. |
| `eggnog_fixurl` | `true` | Install and run the EggNOG URL fixer. |
| `eggnog_fixurl_package` | `eggnog-mapper-fixurl` | URL-fixer package. |
| `eggnog_method` | `diamond` | EggNOG search method. |
| `eggnog_itype` | `metagenome` | EggNOG input type. |
| `eggnog_genepred` | `prodigal` | EggNOG gene prediction method. |
| `eggnog_trans_table` | `11` | Translation table. |
| `eggnog_output_prefix` | `samwise_eggnog` | EggNOG output prefix. |
| `eggnog_extra_args` | `""` | Shell-style arguments appended to EggNOG-mapper; do not duplicate or contradict fixed workflow options. |
| `eggnog_mmseqs_db` | `null` | Optional MMseqs database path. |
| `eggnog_fail_nonfatal` | `false` | Record EggNOG failure and continue when true. |
| `microtrait_env_dir` | `null` | Custom microTrait environment directory. |
| `microtrait_type` | `genomic` | microTrait analysis type. |
| `microtrait_extension` | `fa` | microTrait input extension; implementation requires `fa` (with or without a leading dot). |
| `microtrait_output_prefix` | `samwise_microtrait` | microTrait output prefix. |
| `microtrait_github_repo` | `https://github.com/ukaraoz/microtrait` | microTrait source repository. |
| `microtrait_git_ref` | `null` | Optional microTrait Git ref. |
| `grodon_github_repo` | `jlw-ecoevo/gRodon` | gRodon source repository. |
| `grodon_git_ref` | `null` | Optional gRodon Git ref. |
| `microtrait_auto_download_db` | `true` | Deploy microTrait HMM databases automatically. |
| `microtrait_runner_script` | `<samwise_dir>/bin/local_microTrait_runner.R` | microTrait runner script. |
| `microtrait_merger_script` | `<samwise_dir>/bin/microTrait_merger.R` | microTrait merger script. |
| `microtrait_fail_nonfatal` | `false` | Record unrecoverable microTrait failure and continue when true. |
| `publish_tool_outputs_mode` | `copy` | Tool output publication mode: `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Derived results root, normally the resolved `working_dir`. |
| `outdir` | `<results_dir>/module_6_magannotate` | Derived Module 6 output directory. |
| `module5_final_mag_dir` | `<results_dir>/module_5_subassembly/final_mag_database` | Module 5 MAG candidate directory. |
| `module5_final_mag_manifest` | `<results_dir>/module_5_subassembly/summary/final_mag_database_manifest.tsv` | Automatic Module 5 manifest. |
| `module4_refined_mag_dir` | `<results_dir>/module_4_binrefinement/refined_bins` | Module 4 fallback MAG directory. |
| `checkm2_db_outdir` | `<outdir>/databases/checkm2` | Derived CheckM2 database location. |
| `gtdbtk_db_outdir` | `<outdir>/databases/gtdbtk` | Derived GTDB-Tk database location. |
| `eggnog_db_outdir` | `<outdir>/databases/eggnog` | Derived EggNOG location, overridden by `eggnog_data_path` or `eggnog_data_dir`. |

## Output Contract

Only the five managed tool-result directories are removed at the start of a
run: `drep_out/`, `checkm2/`, `gtdbtk/`, `eggnog/`, and `microtrait/`.
Database directories and `setup/`, `logs/`, and `summary/` are preserved.
Results are published below `<outdir>`:

| Path | Contents |
|---|---|
| `summary/` | Input manifest and statistics, `module6_run_summary.tsv`, per-tool status files, `modified_quality_report.csv`, and the combined quality/taxonomy summary when dRep, CheckM2, and GTDB-Tk all run. |
| `drep_out/` | dRep results and copied `dereplicated_genomes/`. |
| `checkm2/` | CheckM2 result tree. |
| `gtdbtk/` | Published contents of the GTDB-Tk task's `gtdbtk_out/` directory; no extra `gtdbtk_out/` level is retained. |
| `eggnog/` | EggNOG result tree and published Excel annotation when produced. |
| `microtrait/` | microTrait merged trait CSVs plus retained per-genome RDS/CSV provenance. |
| `logs/` | Per-process logs, including database/setup and tool logs. |
| `setup/` | Per-tool setup status files describing environments and resolved databases. |

The input manifest records `mag_id`, source and prepared FASTA paths, contig
counts, and total base pairs. Each tool status TSV records completion or failure,
exit status, output path, and a tool-specific message. `module6_run_summary.tsv`
links the input statistics and all emitted status files.

The combined publication summary joins dRep genome info with the exact GTDB-Tk
`gtdbtk.bac120.summary.tsv` and/or `gtdbtk.ar53.summary.tsv` files. It retains
the GTDB-Tk taxonomy columns and writes one row per dereplicated genome.

## Failure Behavior and Validation

The workflow fails before processing when no annotation tool is selected, when
dRep is enabled without CheckM2, when no usable MAG is found, or when any
required extension is not `fa`. Setup failures and failures from dRep, CheckM2,
or GTDB-Tk are fatal. Database paths are checked for the expected files or
markers before the corresponding tool runs.

EggNOG failures are fatal by default and become `failed_nonfatal` status rows
with exit status zero when `--eggnog_fail_nonfatal true`. microTrait failures
are likewise fatal by default. With `--microtrait_fail_nonfatal true`, missing
or incomplete per-MAG RDS output and merger failures are recorded as
`failed_nonfatal` and the workflow continues. A microTrait runner that exits
nonzero after producing every expected per-MAG RDS file is allowed to continue
to the merger even when nonfatal mode is false.

Before production runs, validate the selected input directory or manifest,
database paths, available disk space, and scheduler memory. Inspect the
published `summary/` status files and `module6_run_summary.tsv`; a completed
Nextflow task alone is not a substitute for checking tool status rows and the
expected output files.
