# Module 3: Binning

![Module 3 workflow](../images/step_3.png)

## Purpose and downstream use

`module_3_binning.nf` assigns assembled scaffolds to metagenome-assembled genome
(MAG) bins. It filters assemblies, maps trimmed reads to the filtered
assemblies with BBMap, and runs any selected combination of QuickBin, MetaBAT2,
and MaxBin2.

The published bins and `summary/binning_manifest.tsv` are the Module 3 contract
for downstream workflows. Module 4 uses these bins for refinement. Select more
than one binner when using the downstream bin-refinement workflow; running one
binner is supported, but does not provide multiple binner results for
comparison/refinement.

## Usage

```bash
nextflow run module_3_binning.nf \
  --working_dir ./output_samwise \
  --threads 6 \
  --metabat2 \
  --quickbin \
  --maxbin2 \
  --min_scaffold_length 2500
```

At least one of `--metabat2`, `--quickbin`, or `--maxbin2` is required. Use
`--` for additional Nextflow flags.

## Inputs and manifest discovery

By default, Module 3 discovers inputs relative to the resolved results root:

| Input | Default path | Contents |
|---|---|---|
| Regular assemblies | `<results_dir>/module_2_readassembly/summary/assembly_manifest.tsv` | Module 2 assembly records, including the `renamed_fasta` path. |
| Regular trimmed reads | `<results_dir>/module_1_readtrimming/summary/trimmed_manifest.tsv` | Module 1 read records used for mapping. |
| Coassembly assemblies | `<results_dir>/module_2b_coassembly/summary/assembly_manifest.tsv` | Module 2B coassembly records, including `renamed_fasta`. |
| Coassembly trimmed reads | `<results_dir>/module_2b_coassembly/summary/coassembly_trimmed_manifest.tsv` | Module 2B reads used for coassembly and mapping. |

`results_dir` is derived from `working_dir`, and the output directory is
`<results_dir>/module_3_binning`. `samwise_dir` identifies the SAMWISE source
tree and defaults to the directory containing the launched workflow; it is
independent of `working_dir`.

Use these options to override discovery or exclude an input family:

| Option | Default | Description |
|---|---|---|
| `input_assembly_manifest` | `null` | Explicit regular Module 2 assembly manifest. |
| `input_trimmed_manifest` | `null` | Explicit regular Module 1 trimmed-read manifest. |
| `input_coassembly_assembly_manifest` | `null` | Explicit Module 2B coassembly assembly manifest. |
| `input_coassembly_trimmed_manifest` | `null` | Explicit Module 2B coassembly trimmed-read manifest. |
| `include_module2` | `true` | Discover/use regular Module 2 assemblies when available. |
| `include_module2b` | `true` | Discover/use Module 2B coassemblies when available. |

Explicit manifest paths are checked and must exist. Automatically discovered
assembly manifests require their corresponding trimmed-read manifest. Relative
paths are resolved from the Nextflow launch directory. At least one assembly
manifest and one read manifest must be selected.

Assembly manifests are tab-separated and must provide at least the fields used
by Module 3, including `sample_id`, `safe_sample_id`, `assembly_sample_id`,
`assembler`, `assembly_mode`, `rarefaction_label`, `assembly_strategy`, and
`renamed_fasta`. Read manifests must provide `sample_id`, `layout`, `read1`,
`read2`, and `interleaved` fields. Supported layouts are `paired` and
`interleaved`. Assemblies and reads are matched by `sample_id`; multiple
assemblies for one sample are retained.

## Tools and versions

With `auto_install` enabled, Module 3 creates or validates separate conda
environments under `<outdir>/conda_envs` (or under `tool_env_dir`): a mapping /
QuickBin environment, a MetaBAT2/JGI environment, and a MaxBin2 environment as
needed. The installer is `mamba` when available, otherwise `conda`. Existing
environments are checked and replaced if required tools are missing.

| Option | Default | Description |
|---|---|---|
| `auto_install` | `true` | Install missing tools with `mamba` or `conda`; when `false`, required environments/tools must already exist. |
| `tool_env_dir` | `null` | Base directory for the tool environments instead of `<outdir>/conda_envs`. |
| `seqkit_version` | `2.8.2` | SeqKit version for assembly length filtering. |
| `bbmap_version` | `39.81` | BBMap/BBTools version for read mapping and QuickBin. |
| `samtools_version` | `1.23.1` | Samtools version for BAM conversion, sorting, and indexing. |
| `metabat2_version` | `2.18` | MetaBAT2 version. |
| `maxbin2_version` | `2.2.7` | MaxBin2 version. |

The setup status is published at `setup/module3_tools_status.env`.

## Threads and filtering

`threads` is a global override. When it is unset, mapping uses
`mapping_threads` and each binner uses `binning_threads`. Filtering and coverage
generation do not define separate user-facing thread controls.

| Option | Default | Description |
|---|---|---|
| `threads` | `null` | Override both mapping and binner task CPUs. |
| `mapping_threads` | `4` | CPUs for BBMap and the Samtools mapping conversion/sort steps. |
| `binning_threads` | `4` | CPUs for QuickBin, MetaBAT2, and MaxBin2. |
| `min_scaffold_length` | `2500` | Minimum contig/scaffold length retained by SeqKit before mapping and binning. |

Assemblies with no contigs at least `min_scaffold_length` are skipped for
mapping and binning. Per-assembly filtering statistics are still written.

## BBMap mapping and coverage

Each retained assembly is mapped separately. Paired reads use `read1` and
`read2`; interleaved reads use `interleaved`. BBMap writes SAM, which is
converted and sorted into an indexed BAM. MetaBAT2 and MaxBin2 also use JGI
depth/coverage output generated from that BAM. QuickBin consumes the sorted BAM
directly.

| Option | Default | Description |
|---|---|---|
| `bbmap_minid` | `0.95` | Minimum BBMap sequence identity. |
| `bbmap_maxindel` | `10` | Maximum allowed BBMap indel length. |
| `bbmap_ambig` | `random` | BBMap handling of ambiguous mappings. |
| `bbmap_mateqtag` | `true` | Pass `mateqtag=t` to BBMap. |
| `bbmap_extra_args` | `""` | Additional BBMap arguments appended to the mapping command. |
| `bbmap_xmx` | `4g` | Java heap value passed as `-Xmx` to BBMap. |

BBMap also uses `build=1`, `fastareadlen=500`, `overwrite=true`, `nodisk=t`,
and the task CPU count. BAM files and indexes are published under
`mapping/`; mapping logs are under `logs/`.

## Binner controls

| Option | Default | Description |
|---|---|---|
| `metabat2` | `false` | Run MetaBAT2. |
| `quickbin` | `false` | Run QuickBin. |
| `maxbin2` | `false` | Run MaxBin2. |
| `metabat2_extra_args` | `""` | Additional MetaBAT2 arguments. MetaBAT2 also receives `-m min_scaffold_length`. |
| `maxbin2_extra_args` | `""` | Additional MaxBin2 arguments. MaxBin2 also receives `-min_contig_length min_scaffold_length`. |
| `quickbin_mincluster` | `50k` | QuickBin minimum cluster size. |
| `quickbin_minseed` | `2500` | QuickBin minimum seed contig length. |
| `quickbin_stringency` | `normal` | QuickBin stringency argument. |
| `quickbin_gzip` | `false` | QuickBin gzip option. |
| `quickbin_chaff` | `false` | Enable QuickBin chaff handling. |
| `quickbin_clade` | `false` | Enable QuickBin clade mode. |
| `quickbin_sketch` | `false` | Enable QuickBin sketch mode. |
| `quickbin_server` | `false` | Enable QuickBin server mode. |
| `quickbin_xmx` | `null` | Optional QuickBin Java heap value; `-Xmx` is added when supplied. |
| `quickbin_extra_args` | `""` | Additional QuickBin arguments. |
| `quickbin_use_positional_bam` | `false` | Pass the BAM positionally instead of as `reads=<bam>`. |

QuickBin uses `min_scaffold_length` as its `mincontig` value. There is no
separate `quickbin_mincontig` parameter in the implementation.

Binner process failures are recorded as non-fatal statuses in per-binner
statistics where supported. Review the status and logs rather than assuming
that a successful Nextflow completion means every selected binner produced
bins.

## Outputs and publication modes

The output root is `<working_dir>/module_3_binning` unless the derived paths are
overridden. Publication modes accept `symlink`, `copy`, or `move`.

| Option | Default | Description |
|---|---|---|
| `publish_filtered_assemblies_mode` | `symlink` | Publication mode for `filtered_assemblies/`. |
| `publish_bam_mode` | `symlink` | Publication mode for `mapping/` BAMs and indexes. |
| `publish_bins_mode` | `symlink` | Publication mode for binner output bins. |
| `results_dir` | `<working_dir>` | Internal results root derived from `working_dir`. |
| `module1_outdir` | `<results_dir>/module_1_readtrimming` | Expected Module 1 output directory. |
| `module2_outdir` | `<results_dir>/module_2_readassembly` | Expected Module 2 output directory. |
| `module2b_outdir` | `<results_dir>/module_2b_coassembly` | Expected Module 2B output directory. |
| `outdir` | `<results_dir>/module_3_binning` | Module 3 output directory. |

Published paths include:

| Path | Contents |
|---|---|
| `filtered_assemblies/` | Length-filtered FASTA files. |
| `mapping/` | Sorted, indexed BBMap BAM files. |
| `coverage/metabat2_depth/` | JGI depth files for MetaBAT2. |
| `coverage/maxbin2_abundance/` | MaxBin2 abundance files derived from depth output. |
| `coverage/quickbin_cov/` | QuickBin coverage files. |
| `bins/metabat2/<binning_id>/` | MetaBAT2 bins. |
| `bins/maxbin2/<binning_id>/` | MaxBin2 bins. |
| `bins/quickbin/<binning_id>/` | QuickBin bins. |
| `summary/` | Combined manifests and summary tables. |
| `summary/per_binner_stats/` | Per-binner statistics tables. |
| `logs/` | Filtering, mapping, coverage, and binner logs. |
| `setup/` | Tool setup status. |

The main summary files are:

- `summary/binning_manifest.tsv`: one row per discovered bin, with binner,
  bin ID, bin size, and published FASTA path.
- `summary/binning_stats_summary.tsv`: per-assembly and per-binner counts,
  total bin base pairs, largest bin, exit status, and status message.
- `summary/filtered_assembly_stats_summary.tsv`: original and retained contig
  counts and base pairs for each assembly.

## Validation checklist

Before starting, verify that:

1. At least one binner flag is enabled.
2. At least one assembly manifest and its matching trimmed-read manifest are
   available after discovery or explicit overrides.
3. Assembly rows contain a non-empty `renamed_fasta` path.
4. Read rows use `paired` or `interleaved` layout and provide the required read
   path fields.
5. `mamba` or `conda` is available when `auto_install=true`, or the configured
   tool environments already contain the required executables when
   `auto_install=false`.
6. `summary/binning_manifest.tsv`, `summary/binning_stats_summary.tsv`, the
   per-binner statistics, and logs are reviewed before passing bins to Module 4.

The workflow fails early for missing explicit manifests, missing paired
manifests, unsupported layouts, missing required tools, or missing final
binning statistics. An assembly with no retained contigs is skipped, and a
binner that exits non-zero can be represented as a non-fatal status; inspect the
published summaries and logs for those cases.
