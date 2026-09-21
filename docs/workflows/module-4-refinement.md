# Module 4: Bin refinement

![Module 4 workflow](../images/step_4.png)

`module_4_binrefinement.nf` refines bins produced by multiple binners.

```bash
nextflow run module_4_binrefinement.nf \
  --working_dir ./output_samwise \
  --threads 36
```

This module requires compatible outputs from more than one binning method.
