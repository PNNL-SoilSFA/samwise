# Module 2: Assembly

![Module 2 workflow](../images/step_2.png)

`module_2_readassembly.nf` assembles trimmed reads. MEGAHIT and metaSPAdes can
be selected according to the available compute resources.

```bash
nextflow run module_2_readassembly.nf \
  --working_dir ./output_samwise \
  --threads 36 \
  --memory_gb 0 \
  --megahit \
  --metaspades
```

For grouped coassembly, use `module_2b_coassembly.nf` with a coassembly
manifest:

```bash
nextflow run module_2b_coassembly.nf \
  --working_dir ./output_samwise \
  --coassembly_groups ./coassembly_manifest.txt \
  --threads 36
```
