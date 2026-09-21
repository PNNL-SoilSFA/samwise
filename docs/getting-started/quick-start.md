# Quick start

![SAMWISE quickstart](../images/quick_start.png)

The following is a complete SLURM-oriented quick start for the main workflow.
Run it from the SAMWISE source directory, replace the account ID, and adjust paths and
resources for your allocation. Remove `-c`, `--slurm_account`, and the
parallelization option when SLURM is not available; these options are optional.
The included `bin/.config` files support SLURM-based parallel execution without
modification.

Reads must use one of the supported `_R1`/`_R2`, `_1`/`_2`, or `_interleaved`
naming forms and have `.fq` or `.fastq` extensions, optionally gzipped.

```bash
# Pre-process reads
nextflow run module_0_readprocess.nf \
  -c ./bin/module_0_slurm.config \
  --working_dir ./samwise-main \
  --slurm_account ChargeAccountID \
  --max_parallel_fastqc 10 \
  --input_dir ./reads_dir \
  --threads 36

# Trim reads
nextflow run module_1_readtrimming.nf \
  -c ./bin/module_1_slurm.config \
  --working_dir ./samwise-main \
  --slurm_account ChargeAccountID \
  --max_parallel_trimming 10 \
  --threads 36

# Assemble reads in parallel across SLURM nodes
nextflow run module_2_readassembly.nf \
  -c ./bin/module_2_slurm.config \
  --working_dir ./samwise-main \
  --slurm_account ChargeAccountID \
  --max_parallel_assemblies 8 \
  --threads 36 \
  --memory_gb 0 \
  --megahit \
  --metaspades \
  --rarefied_assembly TRUE \
  --rarefaction_splits 2

# Optional: coassemble each group in parallel across SLURM nodes. This requires
# ./coassembly_manifest.txt; omit this block for a standard non-coassembly run.
nextflow run module_2b_coassembly.nf \
  -c ./bin/module_2b_slurm.config \
  --working_dir ./samwise-main \
  --coassembly_groups ./coassembly_manifest.txt \
  --slurm_account ChargeAccountID \
  --threads 36 \
  --memory_gb 0 \
  --max_parallel_coassemblies 4

# Bin assemblies
nextflow run module_3_binning.nf \
  --working_dir ./samwise-main \
  --threads 36 \
  --metabat2 \
  --quickbin \
  --maxbin2 \
  --min_scaffold_length 2500

# Refine MAGs; this requires more than one binner
nextflow run module_4_binrefinement.nf \
  --working_dir ./samwise-main \
  --threads 36

# Run a subtractive assembly. Module 5 does not consume coassembly_groups.
nextflow run module_5_subassembly.nf \
  -c ./bin/module_5_slurm.config \
  --working_dir ./samwise-main \
  --slurm_account ChargeAccountID \
  --max_parallel_subassembly 8 \
  --threads 36 \
  --megahit \
  --metaspades \
  --secondpass_metabat2 true \
  --secondpass_quickbin true \
  --secondpass_maxbin2 true

# Annotate MAGs
nextflow run module_6_magannotate.nf \
  --working_dir ./samwise-main \
  --run_checkm2 true \
  --run_gtdbtk true \
  --run_eggnog true \
  --run_microtrait true \
  --threads 36

# Optional downstream modules
nextflow run module_7_gems.nf \
  --working_dir ./samwise-main

nextflow run AuxModule_1_assemblyAnnotate.nf \
  --working_dir ./samwise-main \
  --threads 36 \
  --min_scaffold_bp 1000

nextflow run AuxModule_2_mvp.nf \
  --working_dir ./samwise-main \
  --mvp_modules "0,1,2,3,6" \
  --threads 36 \
  --memory_gb 0
```

For Module 2 and Module 2B, `--memory_gb 0` lets the SLURM profile assign a
whole node and detect available memory. Set the corresponding max-parallel
option to the number of jobs allowed by your allocation. See the module notes
for the `coassembly_manifest.txt` format.

To save compute time, pre-download the Module 6 databases:

- CheckM2: <https://zenodo.org/records/14897628>
- GTDB-Tk: <https://ecogenomics.github.io/GTDBTk/installing/index.html>
- eggNOG: <https://github.com/eggnogdb/eggnog-mapper>; use
  `download_eggnog_data.py --data_dir /path/to/eggnog-data`

Pass pre-downloaded database paths with `--checkm2_db_path`,
`--gtdbtk_data_path`, and `--eggnog_data_path`, respectively:

```bash
--checkm2_db_path /path/to/uniref100.KO.1.dmnd
--gtdbtk_data_path /path/to/gtdbtk/database_directory
--eggnog_data_path /path/to/eggnog/database_directory
```

If only one binner is run, Module 4 cannot perform its normal refinement
workflow. If MAGScoT is not used, pass the Module 3 MAG manifest directly to
Module 6 with `--input_mag_manifest`.

If Module 5's published `final_mag_database` directory is unexpectedly empty,
treat publication as incomplete. Check
`./samwise-main/module_5_subassembly/summary/final_mag_database_manifest.tsv`,
`final_mag_database_stats.tsv`, and the build log before inspecting the source
manifest and FASTAs. Module 6 falls back to Module 4 only when Module 5 is
absent or contains no FASTAs; an empty GTDB-Tk result directory should likewise
be investigated using the Module 6 summary and logs.

If a module halts because of a timeout or another interruption, rerun that
module with `-resume`.

Before a production run, validate a small representative dataset and inspect
the generated manifests, logs, and reports.
