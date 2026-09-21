# Module 6: MAG annotation

![Module 6 workflow](../images/step_6.png)

`module_6_magannotate.nf` annotates final MAGs with CheckM2, GTDB-Tk,
EggNOG-mapper, and optional microTrait analysis.

```bash
nextflow run module_6_magannotate.nf \
  --working_dir ./output_samwise \
  --run_checkm2 true \
  --run_gtdbtk true \
  --run_eggnog true \
  --run_microtrait true \
  --threads 36
```

Large external databases can be downloaded before execution and supplied with
the database path options documented in the root README.
