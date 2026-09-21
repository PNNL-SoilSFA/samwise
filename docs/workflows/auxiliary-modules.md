# Auxiliary modules

![Auxiliary modules](../images/AuxModules.png)

## Assembly annotation

`AuxModule_1_assemblyAnnotate.nf` annotates assemblies from Module 2, Module
2B, or Module 5 with EggNOG-mapper.

```bash
nextflow run AuxModule_1_assemblyAnnotate.nf \
  --working_dir ./output_samwise \
  --threads 36 \
  --min_scaffold_bp 1000
```

## MVP viral analysis

`AuxModule_2_mvp.nf` prepares assemblies and invokes selected MVP stages.

```bash
nextflow run AuxModule_2_mvp.nf \
  --working_dir ./output_samwise \
  --mvp_modules "0,1,2,3,6" \
  --threads 36 \
  --memory_gb 0
```

The auxiliary modules discover compatible upstream outputs under the working
directory. See the root README for all stage and database options.
