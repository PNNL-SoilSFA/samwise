# Module 0: Read processing

![Module 0 workflow](../images/step_0.png)

`module_0_readprocess.nf` validates input reads, detects paired-end or
interleaved layouts, optionally validates FASTQ structure, and runs FastQC.

```bash
nextflow run module_0_readprocess.nf \
  --input_dir ./reads_dir \
  --working_dir ./output_samwise \
  --threads 6
```

Important options include:

- `--input_dir`: directory containing input reads
- `--threads`: global thread allocation
- `--skip_validate`: skip FASTQ validation when inputs are already trusted
- `--auto_install`: allow installation of required tools
- `--working_dir`: results root
- `--samwise_dir`: source directory when source and results are separate

The module writes a read manifest and FastQC outputs below
`<working_dir>/module_0_readprocess`.
