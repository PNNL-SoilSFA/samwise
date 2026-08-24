#!/usr/bin/env Rscript
### Merge microTrait per-genome RDS files in a folder into summary CSVs
# RED 2024; robert.danczak@pnnl.gov

# This script reads the per-genome microTrait result objects
# (<genome>.microtrait.rds) produced by local_microTrait_runner.R and merges
# them into three tidy CSVs: genome features, growth traits, and functional
# (granularity-3) traits.

# Usage:
#   Rscript microTrait_merger.R <rds_dir> [output_dir]
#     <rds_dir>    directory containing the microTrait *.microtrait.rds files
#     [output_dir] where to write the merged CSVs (defaults to <rds_dir>)

# requirements
# (These specific packages cover map_df [purrr], gather/select [tidyr/dplyr]
#  and the %>% pipe; loading them individually rather than the full tidyverse
#  keeps the conda environment lean.)
library(dplyr)
library(tidyr)
library(purrr)

# parse arguments
args = commandArgs(trailingOnly = TRUE)
source_dir = args[1]                                            # dir with .rds files
output_dir = if (length(args) >= 2 && nzchar(args[2])) args[2] else source_dir

if (is.na(source_dir) || !nzchar(source_dir)) {
  stop("microTrait_merger.R: a source directory containing .microtrait.rds files is required as the first argument.")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# define functions
load_microtrait = function(dir){
  # empty object
  mags = NULL

  # load data
  # Match only the per-genome result objects (<genome>.microtrait.rds); this
  # avoids picking up combined outputs the runner also writes here, e.g.
  # *_GenomeSet_Traits.rds, *_Genome_Features.rds, *_mingentime.rds.
  files = list.files(path = dir, pattern = "\\.microtrait\\.rds$", full.names = T)

  if(length(files) == 0){
    stop(paste0("microTrait_merger.R: no *.microtrait.rds files found in '", dir, "'."))
  }

  for(f in files){
    curr.file = readRDS(file = f)

    if(is.null(curr.file$growthrate_ConsistencyHE)){
      curr.file$growthrate_CUBHE = NA
      curr.file$growthrate_ConsistencyHE = NA
      curr.file$growthrate_CPB = NA
      curr.file$growthrate_d = NA
      curr.file$growthrate_LowerCI = NA
      curr.file$growthrate_UpperCI = NA
    }

    mags[[f]] = curr.file
  }

  # change names (basename() is robust to path separators / trailing slashes)
  names(mags) = gsub(".microtrait.rds", "", basename(names(mags)))

  # empty output
  output = NULL

  # parse traits
  output[["genome"]] = map_df(
    lapply(mags, function(x) x$allfeatures),
    ~as.data.frame(.x),
    .id="MAG"
  ) %>%
    select(-Genome)

  output[["growth"]] = map_df(
    lapply(mags, function(x) data.frame(growthrate_CUBHE =
                                          x$growthrate_CUBHE,
                                        growthrate_ConsistencyHE =
                                          x$growthrate_ConsistencyHE,
                                        growthrate_CPB =
                                          x$growthrate_CPB,
                                        growthrate_d =
                                          x$growthrate_d,
                                        growthrate_LowerCI =
                                          x$growthrate_LowerCI,
                                        growthrate_UpperCI =
                                          x$growthrate_UpperCI,
                                        ogt =
                                          x$ogt)),
    ~as.data.frame(.x),
    .id = "MAG"
  ) %>%
    gather(feature, value, -MAG)

  output[["functional"]] = map_df(
    lapply(mags, function(x) x$trait_counts_atgranularity3),
    ~as.data.frame(.x),
    .id="MAG"
  )

  # return
  return(output)
}

# run merger
microtrait_merge = load_microtrait(source_dir)

# write csvs
write.csv(microtrait_merge$genome,
          file.path(output_dir, "microtrait_genome_traits.csv"),
          quote = F,
          row.names = F)
write.csv(microtrait_merge$growth,
          file.path(output_dir, "microtrait_growth_traits.csv"),
          quote = F,
          row.names = F)
write.csv(microtrait_merge$functional,
          file.path(output_dir, "microtrait_functional_traits.csv"),
          quote = F,
          row.names = F)

cat("microTrait merge complete. Wrote 3 CSVs to", output_dir, "\n")
