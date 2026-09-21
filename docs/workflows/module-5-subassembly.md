# Module 5: Subtractive assembly

![Module 5 workflow](../images/step_5.png)

`module_5_subassembly.nf` performs grouped subtractive assembly and can run a
second binning and refinement pass.

```bash
nextflow run module_5_subassembly.nf \
  --working_dir ./output_samwise \
  --coassembly_groups ./coassembly_manifest.txt \
  --threads 36 \
  --megahit \
  --metaspades \
  --secondpass_metabat2 true \
  --secondpass_quickbin true \
  --secondpass_maxbin2 true
```
