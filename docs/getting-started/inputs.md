# Input preparation

Module 0 validates read files and detects paired-end or interleaved layouts.
Supported extensions are `.fq`, `.fastq`, `.fq.gz`, and `.fastq.gz`.

Use one of these naming patterns:

```text
reads_dir/
├── SampleA_R1.fastq.gz
├── SampleA_R2.fastq.gz
├── SampleB_1.fq
├── SampleB_2.fq
└── SampleC_interleaved.fastq
```

Sample identifiers are derived from the filename. Avoid periods and other
delimiters inside the sample identifier because SAMWISE normalizes names for
downstream manifests. For example, `Sample.1_interleaved.fastq` and
`Sample.2_interleaved.fastq` can be interpreted as the same sample.

For coassembly, provide a manifest that maps input samples to groups. An
example is available at
[`examples/coassembly_manifest_example.txt`](https://github.com/PNNL-SoilSFA/samwise/blob/main/examples/coassembly_manifest_example.txt).
