# Module 5: Subtractive Assembly

![Module 5 workflow](../images/step_5.png)

`module_5_subassembly.nf` maps Module 1 trimmed reads against the Module 4
refined MAG reference, assembles the unmapped reads, optionally bins the new
assemblies, and performs a final joint Module 4 refinement. It then publishes a
single final MAG database for downstream modules.

## Workflow Stages

The workflow runs these stages in order:

1. Check or install the Module 5 tools.
2. Build a combined FASTA reference from the Module 4 refined MAG manifest.
3. Map each Module 1 sample to that reference with BBMap and retain unmapped
   reads.
4. Run the selected subtractive assemblers for each sample. MEGAHIT uses
   assembly strategy `E`; metaSPAdes uses strategy `F`.
5. Write subtractive mapping, assembly, and Module 3-compatible manifests.
6. If enabled, run a second Module 3 binning pass on the subtractive
   assemblies.
7. Combine the original and subtractive binning manifests and run final joint
   Module 4 refinement when new bins are available.
8. Build `final_mag_database` from the final joint refined MAGs, or from the
   original Module 4 refined MAGs when the second pass is disabled or produces
   no new MAGs.

## Usage

Basic local execution:

```bash
nextflow run module_5_subassembly.nf \
  --working_dir ./output_samwise \
  --threads 20 \
  --megahit true \
  --metaspades true \
  --secondpass_metabat2 true \
  --secondpass_quickbin true \
  --secondpass_maxbin2 true \
  --run_second_pass_binning_refinement true
```

The README example also shows `--coassembly_groups ./coassembly_manifest.txt`.
That manifest is a required input for Module 2B grouped co-assembly, not for
the current Module 5 implementation: Module 5 does not declare, read, or use
`coassembly_groups`. Supplying it to Module 5 therefore has no effect. Use the
Module 2B workflow when grouped co-assemblies are required.

## Prerequisites and Manifests

By default, Module 5 resolves these inputs under `results_dir`:

| Input | Default path | Required contents |
|---|---|---|
| Trimmed reads | `<results_dir>/module_1_readtrimming/summary/trimmed_manifest.tsv` | Module 1 TSV with `sample_id`, `safe_sample_id`, `layout`, and valid `read1`/`read2` or `interleaved` paths. |
| Original bins | `<results_dir>/module_3_binning/summary/binning_manifest.tsv` | Module 3 binning manifest used to retain the original bins during final joint refinement. |
| Refined MAGs | `<results_dir>/module_4_binrefinement/summary/magscot_refined_bins_manifest.tsv` | Module 4 manifest with non-empty, existing `refined_bin_fasta` values. |

Override these paths with `input_trimmed_manifest`,
`input_original_binning_manifest`, and `input_refined_manifest`. Paired reads
must provide both `read1` and `read2`; interleaved reads must provide
`interleaved`. Relative FASTA paths in legacy manifests are resolved relative
to the workflow launch directory.

The external refinement dependencies are resolved from
`<samwise_dir>/dependencies` unless `dependencies_dir` is overridden. The
TIGRFAM HMM, Pfam HMM, and MAGSCOT script can be supplied explicitly with
`tigrfam_hmm`, `pfam_hmm`, and `magscot_script`.

The Module 5 tool environment contains Python, OpenJDK 17, and BBMap, plus
MEGAHIT when `megahit` is enabled and SPAdes/metaSPAdes when `metaspades` is
enabled. With `auto_install=true`, the workflow creates or repairs this
environment with mamba or conda under `<outdir>/conda_envs/module5_tools`
unless `tool_env_dir` is set.

## Assemblers and Second Pass

MEGAHIT is enabled by default and metaSPAdes is disabled by default. At least
one assembler must be enabled; both may be enabled together. The workflow
creates one assembly job per selected assembler and sample. Assembly failures
and missing/unusable contigs are recorded in summaries where possible rather
than silently treated as successful assemblies.

The second pass is enabled by default and all three binners are enabled by
default: MetaBAT2, QuickBin, and MaxBin2. When second pass is enabled, at
least one second-pass binner must be selected. Module 3 is run with the
subtractive trimmed-read and assembly manifests, and its output is combined
with the original Module 3 binning manifest.

If `run_second_pass_binning_refinement` is `false`, Module 5 skips second-pass
binning and final joint refinement and publishes the original Module 4 refined
MAGs to `final_mag_database`. If second-pass binning has no subtractive
assemblies or produces no new bins, the original bins are also used as the
final database. A second-pass nonzero Module 3 exit is treated as non-fatal
when it produces no bins; it is fatal when it produces bins but fails.

## Parameters

| Argument | Default | Description |
|---|---|---|
| `working_dir` | `null` | Results root. Module 5 writes to `<working_dir>/module_5_subassembly`; if omitted, it inherits `samwise_dir`. |
| `samwise_dir` | `null` | SAMWISE source directory containing workflows, `bin/`, and `dependencies/`; defaults to the launched workflow directory. |
| `coassembly_groups` | `null` | README invocation parameter for the two-column Module 2B read/sample-to-group manifest; not consumed by the current Module 5 implementation. |
| `input_trimmed_manifest` | `null` | Module 1 trimmed-read manifest override. |
| `input_original_binning_manifest` | `null` | Module 3 original binning manifest override. |
| `input_refined_manifest` | `null` | Module 4 refined-bin manifest override. |
| `megahit` | `true` | Enable subtractive MEGAHIT assembly. |
| `metaspades` | `false` | Enable subtractive metaSPAdes assembly. |
| `auto_install` | `true` | Install missing tools with mamba or conda; if false, tools must already exist. |
| `tool_env_dir` | `null` | Custom Module 5 conda environment path. |
| `threads` | `null` | Global thread override for mapping, assembly, and nested workflows. |
| `mapping_threads` | `4` | BBMap threads when `threads` is unset. |
| `assembly_threads` | `4` | Assembly threads when `threads` is unset. |
| `bbmap_version` | `39.81` | BBMap version. |
| `megahit_version` | `1.2.9` | MEGAHIT version. |
| `spades_version` | `4.2.0` | SPAdes/metaSPAdes version. |
| `bbmap_extra_args` | `""` | Extra BBMap arguments. |
| `bbmap_minid` | `0.99` | BBMap minimum identity for subtraction. |
| `bbmap_ambig` | `random` | BBMap ambiguous-mapping behavior. |
| `bbmap_xmx` | `null` | Optional BBMap Java heap size, such as `32g`; do not include `-Xmx`. |
| `megahit_preset` | `meta-large` | MEGAHIT preset: `meta-large` or `meta-sensitive`. |
| `megahit_threads` | `null` | MEGAHIT-specific thread override. |
| `metaspades_memory_gb` | `0` | metaSPAdes memory in GB; `0` leaves its memory argument unset. |
| `run_second_pass_binning_refinement` | `true` | Run second-pass binning and final joint refinement. |
| `secondpass_metabat2` | `true` | Enable MetaBAT2 in the second pass. |
| `secondpass_quickbin` | `true` | Enable QuickBin in the second pass. |
| `secondpass_maxbin2` | `true` | Enable MaxBin2 in the second pass. |
| `module3_script` | `${samwise_dir}/module_3_binning.nf` | Module 3 script for second-pass binning. |
| `module4_script` | `${samwise_dir}/module_4_binrefinement.nf` | Module 4 script for final joint refinement. |
| `nextflow_exe` | `nextflow` | Executable used for nested workflows. |
| `secondpass_working_dir` | `null` | Custom second-pass working directory. The default is `<outdir>/second_pass_binning` (not the README's `<outdir>/second_pass`). The nested Module 3 results are written below `<secondpass_dir>/module_3_binning`. |
| `final_joint_working_dir` | `null` | Custom final joint refinement working directory. |
| `dependencies_dir` | `${samwise_dir}/dependencies` | External dependency directory. |
| `tigrfam_hmm` | `null` | TIGRFAM HMM path. |
| `pfam_hmm` | `null` | Pfam HMM path. |
| `magscot_script` | `null` | MAGSCOT script path. |
| `magscot_extra_args` | `""` | Extra MAGSCOT arguments. |
| `publish_reference_mode` | `copy` | Reference publication mode: `symlink`, `copy`, or `move`. |
| `publish_unmapped_mode` | `symlink` | README-only parameter; not consumed by Module 5. Unmapped reads are always published with `copy`. |
| `publish_assemblies_mode` | `symlink` | README-only parameter; not consumed by Module 5. Assemblies are always published with `copy`. |
| `publish_final_mags_mode` | `copy` | Final MAG publication mode: `symlink`, `copy`, or `move`. |
| `results_dir` | `null` | Internal results root derived from the resolved `working_dir`. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Derived Module 1 output directory. |
| `module3_outdir` | `<results_dir>/module_3_binning` | Derived Module 3 output directory. |
| `module4_outdir` | `<results_dir>/module_4_binrefinement` | Derived Module 4 output directory. |
| `outdir` | `<results_dir>/module_5_subassembly` | Derived Module 5 output directory. |
| `secondpass_dir` | `<outdir>/second_pass_binning` | Derived second-pass directory; overridden by `secondpass_working_dir`. The nested Module 3 output path is `<secondpass_dir>/module_3_binning`, including `summary/binning_manifest.tsv`. |
| `final_joint_dir` | `<outdir>/final_joint_refinement` | Derived final joint directory; overridden by `final_joint_working_dir`. |

## SLURM Scheduler Options

Use the opt-in profile with `-c bin/module_5_slurm.config`:

```bash
nextflow run module_5_subassembly.nf \
  -c bin/module_5_slurm.config \
  --working_dir ./output_samwise \
  --slurm_account ChargeAccountID \
  --max_parallel_subassembly 8 \
  --threads 36 \
  --megahit true --metaspades true
```

The profile supports these scheduler parameters:

| Argument | Default | Description |
|---|---|---|
| `slurm_account` | `null` | SLURM account. |
| `slurm_partition` | `null` | SLURM partition. |
| `slurm_qos` | `null` | SLURM QoS. |
| `slurm_extra_options` | `''` | Additional `sbatch` options. |
| `max_parallel_subassembly` | `8` | Maximum concurrent mapping/assembly jobs and executor queue size. |
| `mapping_time` | `24h` | Wall time for each mapping job. |
| `subassembly_time` | `48h` | Wall time for each subtractive assembly job. |
| `subassembly_max_retries` | `2` | Retries for configured scheduler/node failure exit statuses. |
| `slurm_exclusive` | `true` | Request one exclusive node per heavy job. |
| `slurm_request_all_memory` | `true` | Add `--mem=0` to exclusive jobs. |
| `slurm_before_script` | `''` | Optional shell commands run before each SLURM task. |

Only `MAP_READS_TO_REFINED_MAGS` and `ASSEMBLE_SUBTRACTIVE` use SLURM in this
profile. Setup, reference preparation, summaries, manifest combination, and
final MAG publication are local. Second-pass binning and final joint
refinement are also local and use `--threads` on the launch host. The profile
does not set a Nextflow memory directive; keep `metaspades_memory_gb` and
`bbmap_xmx` within the node memory available to each task.

## Outputs and Manifests

Module 5 publishes outputs below `<outdir>`:

| Path | Contents |
|---|---|
| `reference/refined_mags_reference.fa` | Combined refined MAG reference. |
| `assemblies/` | Durable subtractive assembly FASTAs. |
| `unmapped_reads/` | Durable per-sample unmapped FASTQs. |
| `header_maps/` | Original-to-renamed contig header maps. |
| `summary/` | Mapping, assembly, second-pass, final-refinement, cleanup, combined-manifest, and final-MAG summaries. |
| `summary/subtractive_trimmed_manifest.tsv` | Module 3-compatible trimmed-read manifest for the second pass. |
| `summary/subtractive_module3_assembly_manifest.tsv` | Module 3-compatible subtractive assembly manifest. |
| `summary/second_pass_binning_status.tsv` | Second-pass status written by Module 5. The nested Module 3 binning manifest, when produced, is at `<secondpass_dir>/module_3_binning/summary/binning_manifest.tsv`; it is not published directly as `<outdir>/second_pass_binning.tsv`. |
| `final_mag_database/` | Plain `.fa` final MAG FASTAs. |
| `summary/final_mag_database_manifest.tsv` | Final MAG IDs, source refined MAGs, output paths, and FASTA statistics. |
| `summary/final_mag_database_stats.tsv` | Final database mode, copied/missing MAG counts, contigs, and base pairs. |
| `logs/` and `setup/` | Process logs and tool setup status. |

The final database manifest is the downstream publication contract. Its
`final_database_mode` is `final_joint_refinement` when new bins were refined,
or `original_refined_mags_only` when the original Module 4 MAGs were used.

## Validation and Fallback

`final_mag_database/` is populated from the manifest selected by the final
joint-refinement status. A `completed` status selects
`<final_joint_dir>/module_4_binrefinement/summary/magscot_refined_bins_manifest.tsv`;
`skipped_no_new_mags`, `skipped_no_bins`, or
`skipped_second_pass_disabled` selects the original Module 4 manifest at
`<module4_outdir>/summary/magscot_refined_bins_manifest.tsv`. A failed final
joint refinement is fatal and does not fall back to an older MAG set.

The durable publication is performed after the build process by Nextflow's
`publishDir`. If the published directory is unexpectedly empty, treat the
publication as incomplete: check `summary/final_mag_database_manifest.tsv`,
`summary/final_mag_database_stats.tsv`, and `build_final_mag_database.log`,
then inspect the selected source manifest and its referenced FASTAs. Do not
assume that an empty publication should be replaced with the final-joint tree.
Module 6 independently falls back to Module 4 `refined_bins` only when Module
5's published directory is absent or contains no FASTAs; it does not use the
final-joint tree as a special Module 5 fallback.
