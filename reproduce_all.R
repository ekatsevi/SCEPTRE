####################################################
#
# Reproduce all numerical simulations and data 
# analysis from Katsevich and Roeder (2020).
# 
# Note: This script would take weeks to execute 
# sequentially. The most time-consuming parts of 
# this script have been broken up into tasks, which 
# if executed in parallel, would make the computation
# feasible. 
#
# Note: This script reproduces the entire analysis,
# starting from downloading the raw data and ending 
# the figures in the paper. However, parts of the 
# analysis can be reproduced instead by loading the 
# intermediate results files at bit.ly/SCEPTRE.
#
####################################################

# clear workspace
rm(list = ls())

# source files containing functions
source("analysis/run_one_experiment.R")
source("analysis/run_one_precomputation.R")
source("plotting/aux_plotting.R")

# load libraries
packages = c("R.utils", "reshape2", "MASS", "bigstatsr", 
             "VGAM", "sn", "GenomicRanges", "plyranges",
             "ggrepel", "readxl", "scales", "kableExtra", 
             "gridExtra", "tidyverse")
for(package in packages){
  suppressPackageStartupMessages(library(package, character.only = TRUE))    
}

# set base directory and create directories for output files, if necessary
base_dir = "../files" # set this to a directory where you would like all  
                      # the data and intermediate results to reside
directories_rel = c("", "data", "data/raw", "data/raw/CRISPR", "data/raw/ChIP-seq",
                "data/raw/HIC", "data/raw/GeneHancer", "data/processed", "precomp",
                "results", "results/pvalues_per_task", "results/resampled_zvalues",
                "figures")
directories_abs = sapply(directories_rel, function(dirname)(sprintf("%s/%s", base_dir, dirname)))
for(directory in directories_abs){
  if(!dir.exists(directory)){
    dir.create(directory)
  }
}

# download the raw data
source("analysis/download_raw_data.R")

# process the raw data into a form that's easier to compute with
source("analysis/process_raw_data.R")

# as a sanity check, reproduce p-values from Gasperini et al.
source("analysis/reproduce_Gasperini.R")

# run all precomputations per gene
input_mode = "count_tasks"
task = "precomputation_per_gene"
source("analysis/input_file.R")
for(task_index in 1:num_tasks){
  run_one_precomputation_per_gene(task_index, base_dir)
}

# run all precomputations per gRNA
input_mode = "count_tasks"
task = "precomputation_per_gRNA"
source("analysis/input_file.R")
for(task_index in 1:num_tasks){
  run_one_precomputation_per_grna(task_index, base_dir)
}

# run all gene-gRNA association tests
input_mode = "count_tasks"
task = "experiment"
source("analysis/input_file.R")
for(task_index in 1:num_tasks){
  run_one_experiment(task_index, base_dir)
}

# collate the results
source("analysis/collate_results.R")

# analyze the results
source("analysis/supplementary_analyses.R")

# numerical simulation
source("analysis/simulation.R")

# plot the results
source("plotting/load_results_for_plotting.R")
source("plotting/Figure1.R")
source("plotting/Figure2.R")
source("plotting/Figure3.R")
source("plotting/Figure4.R")
source("plotting/FigureS5.R")
source("plotting/FigureS6.R")
source("plotting/FigureS7.R")