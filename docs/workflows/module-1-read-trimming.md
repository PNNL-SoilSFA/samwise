# Module 1: Read trimming

![Module 1 workflow](../images/step_1.png)

`module_1_readtrimming.nf` trims validated reads with `fastp`, runs FastQC on
the trimmed reads, and writes trimming statistics.

```bash
nextflow run module_1_readtrimming.nf \
  --working_dir ./output_samwise \
  --threads 6
```

The module normally discovers the Module 0 manifest from the working directory.
Use `--input_manifest` when an alternate manifest is required. The full
parameter reference is maintained in the [root README](https://github.com/PNNL-SoilSFA/samwise/blob/main/README.md).
