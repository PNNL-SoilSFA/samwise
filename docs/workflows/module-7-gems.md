# Module 7: GEM generation

![Module 7 workflow](../images/step_7.png)

`module_7_gems.nf` prepares genome-scale metabolic models from Module 6
outputs. It can run gapseq adaptation and MEMOTE validation.

```bash
nextflow run module_7_gems.nf \
  --working_dir ./output_samwise
```

Before production use, run a small representative Module 6 output set and
inspect:

- `summary/module_7_gems_manifest.tsv`
- `summary/module_7_gems_summary.tsv`
- gapseq logs
- MEMOTE reports

The selected medium, database paths, and local gapseq/MEMOTE versions should be
recorded with analysis results.
