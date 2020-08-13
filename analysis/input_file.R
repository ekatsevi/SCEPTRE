#####################################################
#
# This script specifies the different methods to be
# applied to the data and how the data analysis is 
# divided across parallel computations.
#
#####################################################

# resampling methods
resampling_methods = tibble(resampling_type = character(), 
                            test_statistic = character(), 
                            dispersion_estimate = character(), 
                            confounder_adjustment_grna = character(),
                            confounder_adjustment_expression = character()) %>%
  add_row(resampling_type = "marginal_permutation", 
          test_statistic = "LOCO_NB",
          dispersion_estimate = "size_raw",
          confounder_adjustment_grna = NA,
          confounder_adjustment_expression = "percent.mito + prep_batch + log(total_umis) + log(guide_count) + log(gene_count)") %>%
  add_row(resampling_type = "conditional_randomization", 
          test_statistic = "LOCO_NB",
          dispersion_estimate = "size_raw",
          confounder_adjustment_grna =       "percent.mito + prep_batch + log(total_umis) + log(guide_count) + log(gene_count)",
          confounder_adjustment_expression = "percent.mito + prep_batch + log(total_umis) + log(guide_count) + log(gene_count)")

likelihood_methods = tibble(test = character(), dispersion_estimate = character(), confounder_adjustment = character()) %>%
  add_row(test = "z_test",
          dispersion_estimate = "size_shrunk",
          confounder_adjustment = "percent.mito + prep_batch + guide_count + offset(log(Size_Factor))") %>%
  add_row(test = "z_test",
          dispersion_estimate = "size_raw",
          confounder_adjustment = "percent.mito + prep_batch + guide_count + offset(log(Size_Factor))") %>%
  add_row(test = "z_test",
          dispersion_estimate = "size_shrunk",
          confounder_adjustment = "percent.mito + prep_batch + log(total_umis) + log(guide_count) + log(gene_count)") %>%
  add_row(test = "z_test",
          dispersion_estimate = "size_raw",
          confounder_adjustment = "percent.mito + prep_batch + log(total_umis) + log(guide_count) + log(gene_count)")

# number of CRT resamples
B = 500

# number of gene-enhancer pairs per experiment
pairs_per_task = 12591   
genes_per_task = 212     
grnas_per_task = 2200   

# test all pairs of gRNAs and expressed genes
all_deg_results = suppressWarnings(read_tsv(sprintf("%s/data/raw/CRISPR/GSE120861_all_deg_results.at_scale.txt", base_dir),
                           col_types = "cddddddccccciiciiccl"))

# divide all (pre)computations among a number of parallel scripts
parameters = all_deg_results %>% 
  rename(gene_id = ENSG, grna_group = gRNA_group) %>% 
  arrange(gene_id, grna_group) %>% 
  select(grna_group, gene_id) 

parameters_gene = parameters %>% select(gene_id) %>% unique()
genes_per_task = min(genes_per_task, nrow(parameters_gene))
num_gene_tasks = ceiling(nrow(parameters_gene)/genes_per_task)
parameters_gene = parameters_gene %>% 
  mutate(gene_precomp_index = 1 + ((row_number()-1) %% num_gene_tasks)) %>%
  group_by(gene_precomp_index) %>%
  mutate(gene_col_index = row_number(), genes_per_task = n()) %>%
  ungroup()

parameters_grna = parameters %>% select(grna_group) %>% unique()
grnas_per_task = min(grnas_per_task, nrow(parameters_grna))
num_grna_tasks = ceiling(nrow(parameters_grna)/grnas_per_task)
parameters_grna = parameters_grna %>% 
  mutate(grna_precomp_index = 1 + ((row_number()-1) %% num_grna_tasks)) %>%
  group_by(grna_precomp_index) %>%
  mutate(grna_col_index = row_number(), grnas_per_task = n()) %>%
  ungroup()

parameters = parameters %>%
  left_join(parameters_gene, by = "gene_id") %>% 
  left_join(parameters_grna, by = "grna_group") 

num_experiments = ceiling(nrow(parameters)/pairs_per_task)

parameters = parameters %>%
  arrange(gene_id) %>%
  mutate(experiment_index = 1 + ((row_number()-1) %% num_experiments))

######## SET INPUT MODE ###################
# Input mode is one of {"prepare_task", "count_tasks"} 
# and determines how this script behaves
if(input_mode == "prepare_task"){
  stopifnot(exists("task_index"))
  if(task == "experiment"){
    parameters = parameters %>% filter(experiment_index == task_index) %>% arrange(gene_id)
  } else if(task == "precomputation_per_gene" | task == "precomputation_per_gene_size_raw") {  
    # restrict to genes
    parameters = parameters %>% 
      filter(gene_precomp_index == task_index) %>%
      select(gene_id, gene_precomp_index, gene_col_index) %>% 
      unique()
  } else if(task == "precomputation_per_grna"){   
    # restrict to grnas
    parameters = parameters %>% 
      filter(grna_precomp_index == task_index) %>%
      select(grna_group, grna_precomp_index, grna_col_index) %>% 
      unique()
  } else{
    print("Invalid task specification!")
    stopifnot(FALSE)
  }
}
if(input_mode == "count_tasks"){
  if(task == "experiment"){
    num_tasks = num_experiments
  } else if(task == "precomputation_per_gene") {  
    num_tasks = num_gene_tasks
  } else if(task == "precomputation_per_grna"){   
    num_tasks = num_grna_tasks
  } else{
    print("Invalid task specification!")
    stopifnot(FALSE)
  }
}