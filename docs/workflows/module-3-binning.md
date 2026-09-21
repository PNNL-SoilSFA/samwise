# Module 3: Binning

![Module 3 workflow](../images/step_3.png)

`module_3_binning.nf` assigns assembled scaffolds to bins using the selected
binning tools.

```bash
nextflow run module_3_binning.nf \
  --working_dir ./output_samwise \
  --threads 36 \
  --metabat2 \
  --quickbin \
  --maxbin2 \
  --min_scaffold_length 2500
```

Use more than one binner when the downstream bin-refinement module is required.
