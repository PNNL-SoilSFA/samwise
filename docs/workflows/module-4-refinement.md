# Module 4: Bin refinement

![Module 4 workflow](../images/step_4.png)

`module_4_binrefinement.nf` collects Module 3 bins, predicts genes with
Prodigal, searches marker proteins with HMMER, runs MAGScoT, and reconstructs
refined MAG FASTAs. The intended input is compatible bin output from more than
one binner. Module 3 normally supplies `quickbin`, `metabat2`, and `maxbin2`.
Using only one binner does not provide the multiple-binner comparison required
for meaningful refinement.

## Usage

```bash
nextflow run module_4_binrefinement.nf \
  --working_dir ./output_samwise \
  --magscot_threshold 0

# Use `--` for additional MAGScoT flags passed through by the workflow.
```

`--threads` is also accepted as an HMMER-only thread override, for example:

```bash
nextflow run module_4_binrefinement.nf \
  --working_dir ./output_samwise \
  --threads 36
```

The original MAGScoT threshold is `0.5`. SAMWISE defaults to `0` so that all
possible MAGs are retained for downstream GTDB and CheckM evaluation; use
`--magscot_threshold 0.5` when the stricter original threshold is desired.
The value must be between `0` and `1`.

## Input Manifest

By default, the workflow reads:

`<results_dir>/module_3_binning/summary/binning_manifest.tsv`

Use `--input_binning_manifest` to provide another manifest. It must be a
tab-separated file with these columns:

| Column | Description |
|---|---|
| `bin_id` | Identifier for the source bin. If empty, a generated identifier is used. |
| `binner` | Binner name, normally `quickbin`, `metabat2`, or `maxbin2`. |
| `bin_fasta` | Path to the bin FASTA file; `.gz` FASTA is supported. |

Rows without `bin_fasta` are ignored. Relative `bin_fasta` paths are resolved
relative to the Nextflow launch directory, not the task directory. Missing
source FASTAs are skipped and counted in `mag_collection_stats.tsv`.

## Processing

1. Gather the manifest's bin FASTAs into `gathered_bins/` and create a combined
   `mag_contigs.fa` plus combined and per-binner contig-to-bin tables.
2. Run Prodigal in metagenome mode (`-p meta`) on the combined FASTA, producing
   protein, nucleotide, and GFF files.
3. Run HMMER `hmmsearch` against the GTDB rel207 TIGRFAM and Pfam databases
   with `--cut_nc`, using `threads` or `hmm_threads` CPUs. These parameters
   control HMMER only; they do not set the CPU count for Prodigal, MAGScoT, or
   the setup process. The resulting marker table is passed to MAGScoT.
4. Run `bin/MAGScoT.py` with the combined contig-to-bin table, marker table,
   GTDB rel207 default marker profiles, the configured threshold, and any
   `magscot_extra_args`.
5. Reconstruct one FASTA per refined MAGScoT assignment and write the refined
   bin manifest and statistics.

## Dependencies

The default dependency root is `${samwise_dir}/dependencies`. The workflow
uses:

- `python`, `pandas`, `r-base`, `r-optparse`, `r-dplyr`, `r-readr`, `r-funr`,
  and `r-digest`.
- HMMER (`hmmsearch`), Prodigal, and GNU Parallel.
- `gtdbtk_rel207_tigrfam.hmm` and `gtdbtk_rel207_Pfam-A.hmm` for marker
  searches.
- `gtdb_rel207_default_markers.tsv` in the selected MAGScoT profiles
  directory.
- `bin/MAGScoT.py`.

If the configured tool environment is absent or fails its checks, the workflow
creates it with `mamba` or `conda` when `auto_install` is `true`. With
`auto_install` set to `false`, the required tools and packages must already be
available. HMM and MAGScoT paths can be overridden explicitly with the path
parameters below.

## Best-Effort Behavior

MAGScoT is intentionally best-effort. Empty contig-to-bin input or an empty
HMM table causes MAGScoT to be skipped with a status record. A nonzero MAGScoT
exit is recorded as `failed_nonfatal` and does not fail the workflow. If no
refined assignments are available, the reconstruction step writes an empty
refined manifest and stats instead of refined FASTAs.

Missing required input manifests, tools, HMM databases, or the MAGScoT script
are not silently accepted. The default profile file is also required: if
contig and HMM inputs are non-empty and
`gtdb_rel207_default_markers.tsv` is missing from the selected profiles
directory, the MAGScoT process exits nonzero and the workflow fails. This
profile validation occurs before MAGScoT's nonzero execution status is treated
as `failed_nonfatal`. Empty contig or HMM input is checked first and skips
MAGScoT without validating the script or profile. Prodigal and HMMER failures
are workflow failures. Always review the Module 4 summary and refined manifest
before running Module 5 or Module 6; those modules require usable refined MAGs.

## Parameters

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 4 outputs are written to `<working_dir>/module_4_binrefinement`. If omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing the workflows, `bin/`, and `dependencies/`. If omitted, it defaults to the directory containing the launched workflow. |
| `input_binning_manifest` | `null` | Input binning manifest, normally Module 3's `summary/binning_manifest.tsv`. |
| `dependencies_dir` | `${samwise_dir}/dependencies` | Directory containing Module 4 dependency files. |
| `tigrfam_hmm` | `null` | Path to the TIGRFAM HMM database. If unset, the workflow searches the dependency directory for the GTDB rel207 file. |
| `pfam_hmm` | `null` | Path to the Pfam HMM database. If unset, the workflow searches the dependency directory for the GTDB rel207 file. |
| `magscot_script` | `null` | Path to the MAGScoT script. If unset, the workflow checks `${samwise_dir}/bin/MAGScoT.py`, then `${samwise_dir}/bin/magscot.py`; it does not search the dependency directories for the script. |
| `magscot_profiles_dir` | `null` | Directory containing MAGScoT marker profiles. If unset, the workflow selects the first existing path in this order: `${dependencies_dir}/`, `${samwise_dir}/dependencies/`, `${dependencies_dir}/MAGScoT_profiles`, then `${samwise_dir}/dependencies/MAGScoT_profiles`. The selected directory must contain `gtdb_rel207_default_markers.tsv`; the lookup does not search within an already-selected directory. |
| `auto_install` | `true` | Automatically install required tools with `mamba` or `conda` when needed. |
| `tool_env_dir` | `null` | Optional custom path for the Module 4 conda environment. |
| `threads` | `null` | HMMER-only thread override; takes precedence over `hmm_threads` and does not change Prodigal, MAGScoT, or setup threads. |
| `hmm_threads` | `8` | HMMER thread count when `threads` is not provided. |
| `r_base_version` | `null` | `r-base` version to install/use. |
| `hmmer_version` | `null` | HMMER version to install/use. |
| `prodigal_version` | `null` | Prodigal version to install/use. |
| `parallel_version` | `null` | GNU Parallel version to install/use. |
| `magscot_extra_args` | `""` | Additional arguments passed directly to MAGScoT. |
| `magscot_threshold` | `0` | MAGScoT minimum completeness threshold; must be between `0` and `1`. |
| `publish_gathered_bins_mode` | `copy` | Publication mode for gathered bins: `symlink`, `copy`, or `move`. |
| `publish_refined_bins_mode` | `copy` | Publication mode for refined bins: `symlink`, `copy`, or `move`. |
| `results_dir` | derived | Internal results root derived from the resolved `working_dir`. |
| `module3_outdir` | `<results_dir>/module_3_binning` | Expected Module 3 output directory. |
| `outdir` | `<results_dir>/module_4_binrefinement` | Module 4 output directory. |

## Outputs and Paths

The primary downstream contract is:

- `summary/magscot_refined_bins_manifest.tsv`, with columns
  `refined_bin_id`, `refined_bin_fasta`, `contig_count`, and `total_bp`.
- `refined_bins/`, containing the reconstructed `*.fa` files referenced by
  that manifest.

Other published outputs include:

- `gathered_bins/` and `summary/mag_collection_manifest.tsv` /
  `mag_collection_stats.tsv`.
- `concat/mag_contigs.fa`, with the combined contig-to-bin table also retained
  in the MAGScoT publication inputs.
- `prodigal/` Prodigal outputs and `hmm/` HMMER outputs.
- `magscot/` MAGScoT inputs and outputs.
- `setup/`, `logs/`, and `summary/` status and run-summary files, including
  `module4_run_summary.tsv`, `prodigal_status.tsv`, `hmm_status.tsv`,
  `magscot_status.tsv`, and `magscot_refined_bins_stats.tsv`.

Published bin directories under `outdir` are cleaned before the run. The
publication modes apply to `gathered_bins/` and `refined_bins/`; status,
summary, logs, and intermediate outputs are published by the workflow in copy
mode.

## Validation

Before continuing downstream, confirm that:

1. `summary/module4_run_summary.tsv` and the individual status files show the
   expected completed or intentionally skipped states.
2. `summary/magscot_refined_bins_manifest.tsv` exists and contains usable
   non-empty `refined_bin_fasta` entries.
3. The referenced FASTAs exist under `refined_bins/` and have contigs.
4. `mag_collection_stats.tsv` does not show an unexpected number of missing
   source bins, and a completed HMMER run reports TIGRFAM/Pfam hits.

An empty refined manifest is a valid best-effort result, but it is not a usable
input for Module 5 or Module 6 and should be investigated before proceeding.
