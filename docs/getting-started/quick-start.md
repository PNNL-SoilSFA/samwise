# Quick start

![SAMWISE quick start](../images/quick_start.png)

The following commands show the usual order for a complete SAMWISE run. Adjust
the input, resource, and scheduler options for your environment.

```bash
nextflow run module_0_readprocess.nf \
  --input_dir ./reads_dir \
  --working_dir ./output_samwise \
  --threads 6

nextflow run module_1_readtrimming.nf \
  --working_dir ./output_samwise \
  --threads 6

nextflow run module_2_readassembly.nf \
  --working_dir ./output_samwise \
  --threads 6 \
  --megahit

nextflow run module_3_binning.nf \
  --working_dir ./output_samwise \
  --threads 6 \
  --metabat2

nextflow run module_4_binrefinement.nf \
  --working_dir ./output_samwise \
  --threads 6

nextflow run module_6_magannotate.nf \
  --working_dir ./output_samwise \
  --threads 6 \
  --run_checkm2 true \
  --run_gtdbtk true \
  --run_eggnog true

nextflow run module_7_gems.nf \
  --working_dir ./output_samwise
```

Use `-resume` to continue a module after an interrupted run. Before production
use, run a small representative dataset and inspect the generated manifests,
logs, and reports.

The full command examples, scheduler profiles, and parameter tables remain in
the [root README](https://github.com/PNNL-SoilSFA/samwise/blob/main/README.md)
while this documentation migration is underway.
