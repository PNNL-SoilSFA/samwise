# Auxiliary modules

![Auxiliary modules](../images/AuxModules.png)

These workflows are run separately from Modules 0-7. `samwise_dir` identifies
the SAMWISE source tree; `working_dir` identifies the results root. If either
is omitted, `samwise_dir` defaults to the launched workflow's directory and
`working_dir` defaults to `samwise_dir`.

## Auxiliary Module 1: Assembly annotation

`AuxModule_1_assemblyAnnotate.nf` functionally annotates assemblies, not MAGs.
It has no Module 1 read-trimming prerequisite. It discovers assembly
directories produced by Module 2, Module 2B, and Module 5, retains scaffolds
at least `min_scaffold_bp` bases long, removes exact and
reverse-complement-equivalent full-length duplicates, combines the retained
scaffolds, and runs EggNOG-mapper.

```bash
nextflow run AuxModule_1_assemblyAnnotate.nf \
  --working_dir ./output_samwise \
  --threads 36 \
  --min_scaffold_bp 1000
```

At least one of these directories must exist:

- `<working_dir>/module_2_readassembly/assemblies`
- `<working_dir>/module_2b_coassembly/assemblies`
- `<working_dir>/module_5_subassembly/assemblies`

The workflow fails if none exists or if `min_scaffold_bp` is less than 1.
Missing optional source directories are skipped. Managed published directories
are cleaned before a new run.

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. If omitted, inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source tree containing the workflow, `bin/`, and `dependencies/`. |
| `min_scaffold_bp` | `1000` | Minimum retained scaffold length. |
| `threads` | `null` | Global EggNOG thread override. |
| `auto_install` | `true` | Install required tools when unavailable. |
| `tool_env_dir` | `null` | Existing general tool environment directory. |
| `conda_pkgs_dir` | `null` | Conda package directory. |
| `publish_filtered_assemblies_mode` | `copy` | Nextflow publish mode for filtered assemblies. |
| `publish_tool_outputs_mode` | `copy` | Nextflow publish mode for EggNOG outputs. |
| `eggnog_env_dir` | `null` | Existing EggNOG environment directory. |
| `eggnog_mapper_version` | `2.1.13` | EggNOG-mapper version used for installation. |
| `eggnog_data_path` | `null` | Explicit EggNOG data directory. |
| `eggnog_data_dir` | `null` | EggNOG data directory to check or populate. |
| `eggnog_auto_download_db` | `true` | Download EggNOG data when no valid database is available. |
| `eggnog_download_args` | `-y` | Arguments passed to the EggNOG database download. |
| `eggnog_fixurl` | `true` | Enable the EggNOG URL-fix package. |
| `eggnog_fixurl_package` | `eggnog-mapper-fixurl` | URL-fix package name. |
| `eggnog_method` | `diamond` | EggNOG search method: `diamond` or `mmseqs`. |
| `eggnog_itype` | `metagenome` | EggNOG input type. |
| `eggnog_genepred` | `prodigal` | EggNOG gene prediction method. |
| `eggnog_trans_table` | `11` | Translation table. |
| `eggnog_output_prefix` | `samwise_assembly_eggnog` | EggNOG output prefix. |
| `eggnog_extra_args` | `""` | Additional arguments passed to EggNOG-mapper. |
| `eggnog_mmseqs_db` | `null` | Explicit MMseqs database path when using `mmseqs`. |
| `eggnog_fail_nonfatal` | `false` | Record EggNOG failure and continue instead of failing the workflow. |
| `outdir` | `<working_dir>/AuxModule_1_assemblyAnnotate` | Auxiliary Module 1 output directory. |

### Database and outputs

A valid database at
`<working_dir>/module_6_magannotate/databases/eggnog` is reused before another
database is downloaded. An explicit `eggnog_data_path` or `eggnog_data_dir`
takes precedence; otherwise the auxiliary database is under
`<outdir>/databases/eggnog`. `diamond` requires a valid DIAMOND database;
`mmseqs` requires a valid MMseqs database. Unsupported methods and missing
databases fail validation or execution.

The output contains:

- `filtered_assemblies/`: retained per-source assembly FASTA files.
- `inputs/eggnog_assembly_scaffolds.fasta`: combined EggNOG input.
- `eggnog/`: EggNOG-mapper output files.
- `summary/assembly_filtering_manifest.tsv`, `scaffold_filtering_manifest.tsv`, and `scaffold_filtering_stats.tsv`.
- `summary/eggnog_scaffold_manifest.tsv` and `eggnog_status.tsv`.
- `summary/assembly_annotation_run_summary.tsv` and cleanup status.
- `logs/`: preparation and EggNOG logs; `setup/` contains setup status.

## Auxiliary Module 2: MVP viral analysis

`AuxModule_2_mvp.nf` prepares compatible inputs from Module 2 individual
assemblies, optional Module 2B coassemblies, and optional Module 5 subtractive
assemblies, then invokes selected MVP stages.

```bash
nextflow run AuxModule_2_mvp.nf \
  --working_dir ./output_samwise \
  --mvp_modules "0,1,2,3,4,5,100" \
  --threads 36 \
  --memory_gb 0
```

### Inputs and manifests

By default, the workflow reads the following paths under `working_dir`:

- Module 2 assemblies: `<working_dir>/module_2_readassembly/summary/assembly_manifest.tsv`.
- Module 1 trimmed reads: `<working_dir>/module_1_readtrimming/summary/trimmed_manifest.tsv`.
- Module 2B assemblies: `<working_dir>/module_2b_coassembly/summary/assembly_manifest.tsv`.
- Module 2B reads: `<working_dir>/module_2b_coassembly/summary/coassembly_trimmed_manifest.tsv`.
- Module 5 assemblies: `<working_dir>/module_5_subassembly/summary/subtractive_assembly_summary_manifest.tsv`.
- Module 5 assembly files: `<working_dir>/module_5_subassembly/assemblies`.
- Module 5 unmapped reads: `<working_dir>/module_5_subassembly/unmapped_reads`.

Override these paths with `assembly_manifest`, `trimmed_manifest`,
`coassembly_manifest`, `coassembly_trimmed_manifest`,
`subtractive_assembly_manifest`, `subtractive_assembly_dir`, and
`subtractive_unmapped_reads_dir`. Individual, coassembly, and subtractive
inputs are enabled by their corresponding inclusion flags. At least one class
must be enabled. With the default `include_individual_assemblies=true`, the
Module 2 assembly manifest and Module 1 trimmed-read manifest must exist and
contain matching, non-empty files. Missing or empty Module 2B and Module 5
manifests are skipped safely even when their inclusion flags remain true;
records that are present must still resolve to non-empty assembly and read
files.

Paired reads are interleaved record-by-record with BBTools only when required;
existing interleaved reads are used directly. The prepared read manifest and
metadata are generated before MVP runs.

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root; defaults to `samwise_dir`, not the source tree unless the two are intentionally the same. |
| `samwise_dir` | `null` | SAMWISE source tree containing the workflow and dependencies; defaults to the launched workflow directory. |
| `mvp_modules` | `0,1,2,3,4,5,100` | Comma-, semicolon-, or whitespace-separated stages. |
| `include_individual_assemblies` | `true` | Include Module 2 assemblies. |
| `include_coassemblies` | `true` | Include Module 2B assemblies when present. |
| `include_subtractive_assemblies` | `true` | Include Module 5 assemblies when present. |
| `assembly_manifest` | `null` | Module 2 manifest override; otherwise `<working_dir>/module_2_readassembly/summary/assembly_manifest.tsv`. |
| `trimmed_manifest` | `null` | Module 1 manifest override; otherwise `<working_dir>/module_1_readtrimming/summary/trimmed_manifest.tsv`. Required for individual assemblies. |
| `coassembly_manifest` | `null` | Module 2B manifest override; otherwise `<working_dir>/module_2b_coassembly/summary/assembly_manifest.tsv`. |
| `coassembly_trimmed_manifest` | `null` | Module 2B read manifest override; otherwise `<working_dir>/module_2b_coassembly/summary/coassembly_trimmed_manifest.tsv`. |
| `subtractive_assembly_manifest` | `null` | Module 5 manifest override; otherwise `<working_dir>/module_5_subassembly/summary/subtractive_assembly_summary_manifest.tsv`. |
| `subtractive_assembly_dir` | `null` | Module 5 assembly directory override; otherwise `<working_dir>/module_5_subassembly/assemblies`. |
| `subtractive_unmapped_reads_dir` | `null` | Module 5 unmapped-read directory override; otherwise `<working_dir>/module_5_subassembly/unmapped_reads`. |
| `threads` | `4` | MVP and preparation thread count. |
| `memory_gb` | `0` | Nextflow memory request for MVP; `0` requests no explicit memory limit. |
| `auto_install` | `true` | Install MVP and BBTools when unavailable. |
| `mvp_version` | `1.1.5` | MVP version. |
| `tool_env_dir` | `null` | MVP environment directory. |
| `bbmap_version` | `40.02` | BBTools/BBMap version. |
| `bbmap_env_dir` | `null` | BBTools environment directory. |
| `install_databases` | `false` | Allow MVP to install databases. |
| `mvp_database_dir` | `null` | MVP database root; default is `<outdir>/00_DATABASES`. |
| `genomad_db_path` | `null` | Explicit geNomad database path. |
| `checkv_db_path` | `null` | Explicit CheckV database path. |
| `skip_check_errors` | `false` | Pass MVP's skip-check-errors option. |
| `min_seq_size` | `0` | MVP minimum sequence size; `0` leaves MVP default. |
| `genomad_relaxed` / `genomad_conservative` | `false` / `false` | Select one geNomad filter mode; both cannot be enabled. |
| `skip_modify_headers` | `false` | Skip MVP header modification. |
| `viral_min_genes` / `host_viral_genes_ratio` | `1` / `1` | MVP Module 02 thresholds. |
| `min_ani` / `min_tcov` / `min_qcov` | `95` / `85` / `0` | MVP Module 03 thresholds. |
| `read_type` | `short` | Read type: `short` or `long`. |
| `unfiltered_protein_file` | `false` | Use the unfiltered protein file option. |
| `interleaved` / `delete_mapping_intermediates` | `true` / `true` | MVP Module 04 options. The workflow always passes `--interleaved` because prepared reads are interleaved. |
| `covered_fraction` | `null` | MVP Module 05 covered-fraction override. |
| `normalization` / `filtration` | `RPKM` / `conservative` | Module 05 values: `RPKM` or `FPKM`; `relaxed` or `conservative`. |
| `functional_fasta_files` | `representative` | MVP Module 06 FASTA selection. |
| `phrogs_evalue` / `phrogs_score` | `0.01` / `60` | PHROGs thresholds. |
| `pfam_evalue` / `pfam_score` | `0.01` / `50` | PFAM thresholds. |
| `functional_ads` / `functional_rdrp` / `functional_dram` | `true` / `false` / `true` | Module 06 functional analyses. |
| `ads_evalue` / `ads_score` | `0.01` / `60` | ADS thresholds. |
| `rdrp_evalue` / `rdrp_score` | `0.001` / `50` | RdRP thresholds. |
| `delete_functional_intermediates` | `true` | Delete Module 06 intermediates. |
| `binning_sample_group` / `read_mapping_sample_group` | `null` / `null` | Module 07 sample-group overrides. |
| `keep_bam` / `delete_binning_intermediates` | `false` / `false` | Module 07 intermediate-file options. |
| `miuvig_identifier` / `miuvig_step` / `miuvig_template` | `null` / `null` / `null` | Module 99 settings. `miuvig_step` is `setup_metadata` or `prep_submission`; the latter requires a template. |
| `force` | `false` | Pass MVP force options. |
| `outdir` | `<working_dir>/AuxModule_2_mvp` | MVP work and output root. MVP stage directories are created directly below this directory. |

### Validation, ordering, and outputs

Valid stage values are `0,1,2,3,4,5,6,7,99,100`. Empty or invalid selections
fail. Selected stages always run in this canonical order, regardless of the
order supplied: `0` (setup), `1` (geNomad/CheckV), `2` (filter), `3`
(clustering), `4` (read mapping), `5` (vOTU table), `6` (functional
annotation), `7` (binning), `99` (MIUViG), then `100` (summary). Omitting an
earlier stage means its required outputs must already exist under
`<outdir>` from a previous run.

The workflow validates `read_type`, `normalization`, and `filtration`, rejects
conflicting geNomad modes, requires MIUViG settings when stage 99 is selected,
and requires a template for `prep_submission`. MVP and BBTools environments
are installed or reused and validated; disabling auto-install requires valid
existing environments. Database installation is disabled by default. MVP uses
`<outdir>/00_DATABASES` as the database root unless `mvp_database_dir` is
supplied. `genomad_db_path` and `checkv_db_path` override the corresponding
databases independently. When `install_databases=false`, the required
databases must already exist; when it is true, MVP may install missing
databases under the selected database root.

Outputs are written directly to `<outdir>` so completed MVP outputs remain
available if a later stage fails. Important managed outputs include:

- `interleaved_reads/`: generated interleaved FASTQ files when needed.
- `metadata/mvp_prepared_reads_manifest.tsv` and `metadata/mvp_metadata.tsv`.
- `summary/mvp_read_preparation_stats.tsv`, input summaries, and completion status.
- `logs/`: setup, preparation, command, and run logs.
- `00_DATABASES/`: MVP databases when installed there.
- MVP's stage directories and files, using the names expected by MVP.
