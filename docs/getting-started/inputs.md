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

SAMWISE derives each sample ID from the complete filename prefix before the
recognized read suffix. For the example above, the sample names are `SampleA`,
`SampleB`, and `SampleC`. Illumina-style paired reads are supported, including
`sample_S1_L001_R1_001.fastq.gz` and `sample_S1_L001_R2_001.fastq.gz`; their
sample ID is `sample_S1_L001`.

This is prefix-based normalization: punctuation in the prefix is retained and
the implementation does not truncate a sample ID at the first period. For
example, these are distinct sample IDs:

```text
Sample.1_interleaved.fastq
Sample.2_interleaved.fastq
```

They normalize to `Sample.1` and `Sample.2`. Make sure each sample has a
unique prefix before `_R1`, `_R2`, `_1`, `_2`, or `_interleaved`, and do not
provide duplicate files for the same normalized sample and read role.

The root README contains an older statement that periods are removed from
sample IDs. The current Module 0 implementation retains them, as described
here.

For coassembly, provide a manifest that maps input samples to groups. An
example is available at
[`examples/coassembly_manifest_example.txt`](https://github.com/PNNL-SoilSFA/samwise/blob/main/examples/coassembly_manifest_example.txt).
