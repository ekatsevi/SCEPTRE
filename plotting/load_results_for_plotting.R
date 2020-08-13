#########################################################
# 
# Load results from all analyses to prepare for plotting.
# 
#########################################################

# association results on real data
original_results = read_tsv(sprintf("%s/results/original_results.tsv", base_dir),
                            col_types = "cclcciiilccciicdddddd")
resampling_results = read_tsv(sprintf("%s/results/resampling_results.tsv", base_dir),
                              col_types = "cclcciiilccciiccddddddd")
likelihood_results = read_tsv(sprintf("%s/results/likelihood_results.tsv", base_dir),
                              col_types = "cccciiilccciiccccddd")

# results on simulated data
simulation_results = read_tsv(sprintf("%s/results/simulation_results.tsv", base_dir), col_types = "dcidddidd")

# supplementary analysis results
grna_confounding = read_tsv(sprintf("%s/results/gRNA_confounding.tsv", base_dir), col_types = "cd")
rejected_pairs_HIC = read_tsv(sprintf("%s/results/rejected_pairs_HIC.tsv", base_dir),
                              col_types = "cciiiciilliid")
TF_enrichments = read_tsv(sprintf("%s/results/TF_enrichments.tsv", base_dir), 
                          col_types = "ccd")
paired_fractions = read_tsv(sprintf("%s/results/TF_paired_enhancer_fractions.tsv", base_dir),
                            col_types = "cidd")

# details of methods implemented
input_mode = ""
task = ""
source("analysis/input_file.R")

# raw and processed data
source("analysis/read_data.R")