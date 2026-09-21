# Module 7: Genome-Scale Metabolic Model Generation

![Module 7 workflow](../images/step_7.png)

`module_7_gems.nf` constructs one genome-scale metabolic model (GEM) for each
MAG represented in the Module 6 eggNOG predictions. It uses gapseq to
reconstruct and gap-fill models, optionally adapts selected models to
MAG-specific empirical growth requirements, and validates final SBML models
with MEMOTE.

## Workflow stages

The workflow runs these stages:

1. `SETUP_GAPSEQ` creates or validates the pinned gapseq/MEMOTE environment.
2. `PREPARE_GAPSEQ_INPUTS` splits the unified Module 6 predicted-protein FASTA
   into one `.faa` file per MAG.
3. `RUN_GAPSEQ_DOALL` runs gapseq `doall` for every MAG using the selected
   medium and organism template.
4. `RUN_GAPSEQ_ADAPT` optionally adapts only MAGs listed in an adaptation
   manifest. Unlisted MAGs retain their `doall` model.
5. `RUN_MEMOTE_SNAPSHOT` and/or `RUN_MEMOTE_RUN` optionally validate each final
   model.
6. `WRITE_OUTPUT_MANIFEST` writes stable per-MAG paths and run summary files.

## Usage

### Default run

This uses Module 6 files below the working directory, the bundled comprehensive
medium, the `Bacteria` template, and one MEMOTE HTML snapshot per final model.

```bash
nextflow run module_7_gems.nf \
  --working_dir ./output_samwise
```

### Archaeal template and custom medium

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

`protein_fasta_dir` is a historic parameter name. It accepts the unified
predicted-protein FASTA **file**, not a directory.

```bash
nextflow run module_7_gems.nf \
  --working_dir ./module_7_results \
  --protein_fasta_dir /absolute/path/to/samwise_eggnog.emapper.genepred.fasta \
  --input_manifest /absolute/path/to/eggnog_input_manifest.tsv
```

## Inputs

### Module 6 inputs

With `working_dir`, the defaults are:

| Input | Default path |
|---|---|
| Unified eggNOG predicted-protein FASTA | `<working_dir>/module_6_magannotate/eggnog/samwise_eggnog.emapper.genepred.fasta` |
| eggNOG input manifest | `<working_dir>/module_6_magannotate/summary/eggnog_input_manifest.tsv` |

If only `samwise_dir` is provided, `working_dir` inherits that directory.
Override either default with `--protein_fasta_dir` or `--input_manifest`.

Both files must be present, regular, and non-empty. The Module 6 manifest must
contain non-empty `mag_id` values. For each FASTA header, the text before the
first `|`, or before whitespace when there is no `|`, must equal a manifest
`mag_id`. MAGs with no matching proteins are recorded as warnings; the run
fails only when no proteins match any MAG.

### Media and template

- `--media_csv` selects the existing, non-empty CSV medium passed to gapseq.
- If it is unset, the bundled comprehensive medium
  `background/media/gapseq_all_nutrients.csv` is used.
- The bundled minimal M9 glucose aerobic medium is
  `background/media/gapseq_M9_glucose_aerobic.csv`; select it explicitly with
  `--media_csv`.
- `--template_organism` must be exactly `Bacteria` (default) or `Archaea`.

### Adaptation manifest

Adaptation is off by default. Enabling it requires both:

```bash
--run_gapseq_adapt true \
--adapt_manifest /absolute/path/to/adapt_manifest.tsv
```

The manifest is tab-separated and must contain `mag_id` and
`adapt_compounds` headers. Additional columns are allowed. Each data row must
have non-empty values, each `mag_id` must occur in the Module 6 manifest, and
each comma-separated compound must match
`cpd#####:(TRUE|FALSE)`.

```tsv
mag_id	adapt_compounds
MAGScoT_cleanbin_000023	cpd00076:TRUE,cpd00027:TRUE
MAGScoT_cleanbin_000038	cpd00027:TRUE
```

`background/test_adapt_manifest.tsv` is a format example. Only listed MAGs
that have mapped proteins enter adaptation; a listed MAG with no mapped
proteins is omitted from the prepared input manifest and remains `doall`-only.
If gapseq reports that a model already grows and produces no `-adapt` files,
Module 7 copies the original `doall` model into the adaptation output and
records `no_changes` in adaptation task metadata. The adaptation task requires
exactly one adapted RDS and one adapted XML, or neither; a partial or multiple
result is an error. Supplying `adapt_manifest` while adaptation is disabled
produces a warning and skips adaptation.

The former `adaptation_compounds` parameter is rejected. A global compound
list is not supported because requirements are MAG-specific.

## Environment and commands

Module 7 requires gapseq `1.4.0` and MEMOTE `0.17.0`. The install
specifications are configurable through `gapseq_package` and `memote_package`,
but the workflow still checks the installed executable versions and will only
run with gapseq `1.4.0` and MEMOTE `0.17.0`. Changing a package specification
does not opt into support for another tool version.

Without an environment override, the workflow manages
`<outdir>/conda_envs/gapseq`, where `<outdir>` is
`<results_dir>/module_7_gems`. With the default `--auto_install true`, a
missing environment is created with `mamba` or `conda`. A stale
workflow-managed environment may be removed and rebuilt automatically.

Use a pre-existing user-managed environment with either:

- `--gapseq_env_dir /path/to/environment` for the exact environment path.
- `--tool_env_dir /path/to/tool-parent` for `/path/to/tool-parent/gapseq`.

User-managed environments must already contain exactly the required versions;
Module 7 never creates, deletes, or modifies them. `gapseq_env_dir` takes
precedence when both overrides are supplied. `--conda_pkgs_dir` selects the
`CONDA_PKGS_DIRS` location. Setup creates that directory even when a
user-managed environment is selected, but does not install packages into or
modify that environment; with no override, the default is the output-local
`conda_pkgs/gapseq` directory.

The fixed command interfaces are:

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

`threads` requests CPUs from Nextflow for each `doall` task; it is not passed
as a gapseq argument. When unset, each `doall` process requests one CPU.
`gapseq_extra_args` is retired and fails preflight validation.

## MEMOTE modes

MEMOTE is enabled by default with `--memote_mode snapshot`.

| `memote_mode` | Command | Output per final model |
|---|---|---|
| `snapshot` | `memote report snapshot` | HTML report |
| `run` | `memote run --ignore-git` | JSON result |
| `both` | Both commands | HTML report and JSON result |

Set `--run_memote false` to skip all MEMOTE tasks. `--memote_extra_args` is
appended to the selected MEMOTE command and must contain only options
supported by MEMOTE `0.17.0`.

## Parameters

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root; outputs go to `<working_dir>/module_7_gems`. Inherits `samwise_dir` when unset. |
| `samwise_dir` | `null` | SAMWISE source tree containing workflows, `bin/`, and `dependencies/`; defaults to the launched workflow directory. |
| `input_manifest` | `null` | Module 6 eggNOG input manifest TSV override. |
| `protein_fasta_dir` | `null` | Module 6 unified predicted-protein FASTA override; despite its name, this is a file. |
| `media_csv` | `null` | User medium CSV; comprehensive bundled medium is used when unset. |
| `template_organism` | `Bacteria` | Exactly `Bacteria` or `Archaea`. |
| `run_gapseq_adapt` | `false` | Enable adaptation for MAGs listed in `adapt_manifest`. |
| `adapt_manifest` | `null` | Adaptation TSV; required when adaptation is enabled. |
| `adaptation_compounds` | retired | Do not supply; global adaptation compounds are rejected. |
| `run_memote` | `true` | Enable MEMOTE validation. |
| `memote_mode` | `snapshot` | `snapshot`, `run`, or `both`. |
| `threads` | `null` | CPUs requested per gapseq `doall`; unset means one. |
| `tool_env_dir` | `null` | Parent of a user-managed `gapseq` environment. |
| `gapseq_env_dir` | `null` | Exact user-managed environment path; takes precedence over `tool_env_dir`. |
| `auto_install` | `true` | Create a missing workflow-managed environment with `mamba` or `conda`. |
| `conda_pkgs_dir` | `null` | Optional package cache for the managed environment. |
| `gapseq_package` | `gapseq=1.4.0` | Conda/mamba install specification; executable validation still requires gapseq `1.4.0`. |
| `memote_package` | `memote=0.17.0` | Conda/mamba install specification; executable validation still requires MEMOTE `0.17.0`. |
| `gapseq_extra_args` | retired | Do not supply; `doall` uses a fixed positional interface. |
| `memote_extra_args` | `""` | Extra options appended to the selected MEMOTE command. |
| `publish_gems_mode` | `copy` | Publishing mode for prepared FASTAs and GEM directories: `copy`, `symlink`, or `move`. |
| `publish_reports_mode` | `copy` | Publishing mode for MEMOTE reports: `copy`, `symlink`, or `move`. |
| `results_dir` | derived | Internal results root derived from resolved `working_dir`; normally do not set. |
| `outdir` | `<results_dir>/module_7_gems` | Internal Module 7 output location; normally do not set. |
| `module6_output_dir` | `<results_dir>/module_6_magannotate` | Default Module 6 output location used to derive inputs. |
| `media_minimal` | bundled M9/glucose CSV | Internal path to the bundled minimal medium. |
| `media_comprehensive` | bundled all-nutrients CSV | Internal path to the default medium. |

## Outputs

All published outputs are beneath `<results_dir>/module_7_gems/`:

```text
module_7_gems/
├── conda_envs/gapseq/                         # workflow-managed environment, if used
├── conda_pkgs/gapseq/                         # CONDA_PKGS_DIRS; created even for user-managed envs
├── inputs/
│   ├── protein_fastas/<mag_id>.faa
│   ├── gapseq_input_manifest.tsv
│   └── gapseq_inputs_stats.tsv
├── setup/
│   └── gapseq_setup_status.env
├── gems/
│   ├── doall/<mag_id>_gapseq_doall/
│   └── adapt/<mag_id>_gapseq_adapt/           # selected MAGs only
├── reports/
│   ├── snapshot/<mag_id>_memote_snapshot.html  # snapshot or both
│   └── run/<mag_id>_memote_run.json            # run or both
├── logs/
│   └── ...
└── summary/
    ├── module_7_gems_manifest.tsv
    └── module_7_gems_summary.tsv
```

The managed environment is created only when no user-managed environment is
supplied. The `CONDA_PKGS_DIRS` directory is created for every run, including
user-managed environments, although no package installation is performed in
that case. `logs/` contains execution logs for `PREPARE_GAPSEQ_INPUTS`, `SETUP_GAPSEQ`,
`RUN_GAPSEQ_DOALL`, `RUN_GAPSEQ_ADAPT`, `RUN_MEMOTE_*`, and
`WRITE_OUTPUT_MANIFEST`.

`summary/module_7_gems_manifest.tsv` has one row per MAG and these columns:

```text
mag_id, protein_fasta,
doall_model_rds, doall_model_xml, doall_rxn_weights, doall_rxn_genes,
doall_reactions_tbl, doall_log,
adapted, adapt_model_rds, adapt_model_xml, adapt_log,
final_model_stage, final_model_rds, final_model_xml,
memote_snapshot_html, memote_snapshot_log,
memote_run_json, memote_run_log
```

Paths are stable published paths, not temporary Nextflow work-directory paths.
`adapted` is `true` only for MAGs that entered the adaptation process, and the
final model is the adaptation result for those MAGs or the `doall` result for
all other MAGs. The adaptation output can still contain the original model
when its status is `no_changes`. MEMOTE columns are empty when the relevant
mode was not run.

`summary/module_7_gems_summary.tsv` reports `total_mags`, `adapted_mags`,
`doall_only_mags`, `memote_snapshot_validated_mags`, and
`memote_run_validated_mags`.

## Validation guidance

Before production use, run Module 7 on a small representative Module 6 output
set. Confirm that:

- `setup/gapseq_setup_status.env` reports gapseq `1.4.0` and MEMOTE `0.17.0`.
- `inputs/gapseq_input_manifest.tsv` and `gapseq_inputs_stats.tsv` account for
  the expected MAGs and protein mappings.
- `summary/module_7_gems_manifest.tsv` has one complete row per expected MAG,
  with the intended final stage and stable paths.
- gapseq logs contain the selected medium, template, and successful required
  model/support artifacts; each `doall` model, XML, reaction weights, gene
  mapping, and reaction table is present and non-empty.
- Adapted MAGs have the expected final model and adaptation metadata. A
  `no_changes` result is valid only when the original model was copied; a
  partial or ambiguous adaptation result must fail.
- MEMOTE reports exist for the selected mode, refer to the final model, and
  have non-empty output files.
- `module_7_gems_summary.tsv` counts match the manifest and the selected
  adaptation and MEMOTE modes.

Record the selected medium, template, environment/version status, adaptation
manifest, and MEMOTE mode with analysis results. The workflow should fail
for missing, non-regular, or empty required input files; invalid templates or
MEMOTE modes; invalid adaptation manifests; retired non-empty parameters; and
unavailable required tool versions. During execution it also fails when no
proteins map to any MAG, a required gapseq artifact is missing or empty, an
adaptation result is partial or ambiguous, a final XML model is not uniquely
identified, or a MEMOTE command fails. The final manifest rejects duplicate
metadata and missing `doall` results.
