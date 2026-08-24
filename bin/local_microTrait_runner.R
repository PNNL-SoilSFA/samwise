#!/usr/bin/env Rscript
# RED 2024; robert.danczak@pnnl.gov

# This script is going to run the microTrait package to identify some traits
# for genomes in the provided folder.

# For this function to be used appropriately, you should have your MAGs placed
# into a common folder which will be parsed into this script. I pulled much of
# this from the microTrait GitHub page, so I don't really take any credit.

# To-do:
# 1) Customize outputs

# parse arguments
args = commandArgs(trailingOnly = T)
fasta_dir = args[1] # arg1 = the directory containing fasta files
out_name = args[2] # arg2 = the output names
type = args[3] # arg3 = either protein or genome
               # if protein, make sure it's in the same directory as the .fa and
               # .fna files (OGT and mingen need these)

# setup
library(microtrait) # load in microtrait
library(dplyr) # for case_when because I'm lazy

# load in specific microtrait functions
# These helpers used to be source()d from a local checkout of the microtrait
# repo (e.g. "~/GitHub Packages/microtrait-022624/R/*.R"). Pulling them from
# the installed package namespace instead removes that hard-coded path and
# works wherever microtrait is installed. getFromNamespace() returns each
# function with its microtrait namespace intact, so the internal helpers they
# call (from protein.R, nucleotide.R, extern.R, utils.R, etc.) still resolve.
extract_features  <- utils::getFromNamespace("extract_features", "microtrait")
run_ogtmodel      <- utils::getFromNamespace("run_ogtmodel", "microtrait")
run.predictGrowth <- utils::getFromNamespace("run.predictGrowth", "microtrait")

# extracting traits based upon the provided type
if(type == "genomic"){
  # genomic traits can leverage the in-built parallel trait identification
  
  ###
  ### list files
  ###
  
  # list fasta files
  fasta.list = list.files(path = fasta_dir,
                          pattern = ".fa",
                          full.names = T)
  
  ###
  ### extracting traits
  ###
  
  # extract traits (spits out .rds in same directory of fasta)
  microtrait_results = extract.traits.parallel(fa_files = fasta.list)
  
  # extract rds file list
  rds_files = unlist(parallel::mclapply(microtrait_results, "[[", "rds_file", 
                                        mc.cores = 10)) # don't allocate too 
                                                        # many threads or else 
                                                        # you run into a massive 
                                                        # slowdown
  
  # combine traits
  genomeset_results = make.genomeset.results(rds_files = rds_files,
                                             ids = sub(".microtrait.rds", "", 
                                                       basename(rds_files)),
                                             ncores = 1)
  
  ###
  ### saving
  ###
  
  # write out combined .rds
  saveRDS(genomeset_results, paste0(unique(dirname(fasta.list)),  "/", out_name,
                                    "_GenomeSet_Traits.rds"))
  
  # finished
  cat("\n")
  cat(paste("The script has finished extracting traits from", 
            length(fasta.list), "genomes."))
  
} else if(type == "protein"){
  # protein traits require a bit of engineering for them to conform
  # the protein traits also need support from genomes if we want to predict OGT 
  # and mingen.
  
  ####
  #### list files ####
  ####
  
  # list fasta files
  fasta.list = list.files(path = fasta_dir,
                          pattern = ".fa$",
                          full.names = T)
  
  
  # list gene-called files
  called.files = list.files(path = fasta_dir,
                            pattern = ".fna$",
                            full.names = T)
  
  
  # list protein files
  protein.files = list.files(path = fasta_dir,
                             pattern = ".faa$",
                             full.names = T)
  
  
  ####
  #### extracting traits ####
  ####
  
  # logging progress
  tictoc::tic(paste0("Running microtrait for ", length(protein.files), 
                     " genomes"))
  
  # extracting traits, but not running OGT/mingen
  microtrait_results = parallel::mclapply(1:length(protein.files), function(i){
    returnList = extract.traits(protein.files[i], fasta_dir, type = "protein",
                                growthrate_predict = F, optimalT_predict = F)
    returnList}, 
  mc.cores = floor(parallel::detectCores()*0.7)) # pulled this from the 
                                                 # package itself, but 
                                                 # modified it for proteins
  
  # logging progress
  tictoc::toc(log = "TRUE")
  
  # extract rds file list
  rds_files = unlist(parallel::mclapply(microtrait_results, "[[", "rds_file", 
                                        mc.cores = 10)) # don't allocate too 
                                                        # many threads or else 
                                                        # you run into a massive 
                                                        # slowdown
  
  # combine traits
  genomeset_results = make.genomeset.results(rds_files = rds_files,
                                             ids = sub(".microtrait.rds", "", 
                                                       basename(rds_files)),
                                             growthrate = F, optimumT = F,
                                             ncores = 1)
  
  ####
  #### saving ####
  ####
  
  # write out combined .rds
  saveRDS(genomeset_results, paste0(unique(dirname(fasta.list)),  "/", out_name,
                                    "_GenomeSet_Traits.rds"))
  
  # finished
  cat("\n")
  cat(paste("The script has finished extracting traits from", 
            length(fasta.list), "genomes."))
  
  ####
  #### running OGT ####
  ####
  
  # logging progress
  tictoc::tic(paste0("Running OGT for ", length(protein.files), " genomes"))
  
  # extract genome features
  genome_features = parallel::mclapply(gsub(".fa", "", fasta.list), function(curr.bin){
    # extracting features used in OGT prediction
    temp = extract_features(genome_file = paste0(curr.bin, ".fa"),
                            cds_file = paste0(curr.bin, ".fna"),
                            proteins_file = paste0(curr.bin, ".faa"))
    
    # output
    temp}, 
    mc.cores = floor(parallel::detectCores()*0.7))
  
  names(genome_features) = gsub("/.*/", "", 
                                gsub(".fa", "", fasta.list))
  
  # logging progress
  tictoc::toc(log = "TRUE")
  
  # run the OGT model
  ogt_out = sapply(genome_features, run_ogtmodel)
  ogt_out = data.frame(user_genome = names(ogt_out),
                       OGT = ogt_out)
  
  # write out results
  saveRDS(genome_features, paste0(unique(dirname(fasta.list)),  "/", out_name,
                                  "_Genome_Features.rds"))
  write.csv(ogt_out,  paste0(unique(dirname(fasta.list)),  "/", out_name,
                             "_OGT_results.csv"), row.names = F)
  
  ####
  #### running minimum generation time ####
  ####
  
  # extract genome features
  mingentime = parallel::mclapply(gsub(".fa", "", fasta.list), function(curr.bin){
    # extracting features used in OGT prediction
    temp = run.predictGrowth(cds_file = paste0(curr.bin, ".fna"),
                            proteins_file = paste0(curr.bin, ".faa"))
    
    # output
    temp}, 
    mc.cores = floor(parallel::detectCores()*0.7))
  
  # write out results
  saveRDS(mingentime, paste0(unique(dirname(fasta.list)),  "/", out_name,
                                  "_mingentime.rds"))
  
}
