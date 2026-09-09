# module_7_goals.md

## Brief Background

SAMWISE is an automated, end-to-end metagenomic read processing program. Please see samwise_README.md for details.

There are currently 6 NextFlow DSL2 workflows in the SAMWISE pipeline.


## Goal

Our goal is to develop an additional module_7_gems that extends the SAMWISE pipeline to automated generation of Genome-scale Metabolic Models (GEMs) using the gapseq algorithm (https://github.com/jotech/gapseq) and includes standardized evaluation with MEMOTE (https://github.com/opencobra/memote).

We want to build the additional module using the same conventions as established SAMWISE modules (Nextflow workflows).
The repository for SAMWISE is: https://tanuki.pnnl.gov/josue.rodriguez/samwise.git
Please look at the structure of the Nextflow workflows to understand what we are trying to mimic here. Example: ../module_6_magannotate.nf Each module references the manifest.tsv file from the output of the prior module to execute the next analysis step in the workflow.

## Module requirements

Module 7 will entail the following steps:
1. Create a new mamba or conda environment
2. Install gapseq and memote
3. Activate the environment
4. Run the gapseq 'doall' command to generate and gapfill genome-scale metabolic models for all MAGs (bins) generated in prior steps of the SAMWISE workflow
5. If the user has empirical growth data, they can choose to refine the model using the gapseq 'adapt' command
6. Run memote 'report snapshot' and memote 'run' to generate reports for each genome-scale metabolic model.

Gapseq doall requires as input:
1. Filepath to fasta file of predicted protein sequences (this output is generated from module 6, and can be found within eggnog/ subdirectory of the output)
2. Filepath to media csv file for gapfilling [should be user input, default to complex media file if missing]
3. Which template to use (Bacteria|Archaea) [should be user input, default to Bacteria if missing]

Gapseq adapt requires as input:
1. Filepath to gapseq-generated genome-scale metabolic model (.RDS extension, not the .xml file)
2. A comma-separated string of compounds that the model should be adjusted to using or not using (format == cpd#####:TRUE|FALSE) [should be user input]
3. Filepath to rxnWeights.RDS file from generated model
4. Filepath to rxnXgenes.RDS file from generated model
5. Filepath to all-Reactions.tbl file from generated model

Memote requires as input:
1. Filepath to generated model
2. Filename for generated html report

Example scripts for running gapseq doall:
module load python/miniconda25.5.1
source /share/apps/python/miniconda25.5.1/etc/profile.d/conda.sh
conda create -c conda-forge -c bioconda -n gapseq
conda activate gapseq
gapseq doall ./faa_protein_files/Rhodococcus_MSC16_mgL.faa M9-all5-medium.csv Bacteria


Example scripts for running gapseq adapt:
module load python/miniconda25.5.1
source /share/apps/python/miniconda25.5.1/etc/profile.d/conda.sh
conda activate gapseq
gapseq adapt -m ./gapseq_doall_output/Ensifer_adhaerens_MSC03_mgL.RDS -w cpd00027:TRUE,cpd00122:TRUE,cpd00023:TRUE,cpd00053:TRUE,cp
d00794:TRUE -c ./gapseq_doall_output/Ensifer_adhaerens_MSC03_mgL-rxnWeights.RDS -g ./gapseq_doall_output/Ensifer_adhaerens_MSC03
_mgL-rxnXgenes.RDS -b ./gapseq_doall_output/Ensifer_adhaerens_MSC03_mgL-all-Reactions.tbl


Example scripts for running memote report snapshot and memote run
memote report snapshot --filename "out_report.html" ./gapseq_doall_output/Rhodococcus_MSC16_mgL-adapt.xml 
memote run ./gapseq_doall_output/Rhodococcus_MSC16_mgL-adapt.xml
