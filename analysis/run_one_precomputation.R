####################################################
#
# This file defines functions to run precomputations
# per gene and per gRNA. 
#
####################################################

# precompute distillations for each gene
run_one_precomputation_per_gene = function(task_index, base_dir){
  # set up all analysis parameters
  task = "precomputation_per_gene"
  input_mode = "prepare_task"
  source("analysis/input_file.R", local = TRUE)
  
  # read in the data
  source("analysis/read_data.R", local = TRUE)
  
  # extract confounder adjustments to apply
  confounder_adjustments_expression = resampling_methods %>% 
    pull(confounder_adjustment_expression) %>% unique()
  
  # subset data to cells containing at least one guide RNA
  df_full = confounders_full
  df_full$cell_name = cell_names
  df = df_full %>% filter(guide_count > 0)
  cells_to_keep = df %>% pull(cell_name)
  num_cells_to_keep = length(cells_to_keep)
  
  # initialize file-backed matrices to store offsets
  genes_per_task = nrow(parameters)
  num_confounder_adjustments = length(confounder_adjustments_expression)
  offsets = vector("list", num_confounder_adjustments)
  names(offsets) = confounder_adjustments_expression
  for(confounder_adjustment_index in 1:length(confounder_adjustments_expression)){
    confounder_adjustment_expression = confounder_adjustments_expression[confounder_adjustment_index]
    offsets_filename = sprintf("%s/precomp/offsets_index_%d_confounders_%d", 
                               base_dir, task_index, confounder_adjustment_index)
    offsets[[confounder_adjustment_index]] = FBM(nrow = num_cells_to_keep, 
                                                 ncol = genes_per_task, 
                                                 type = "double",
                                                 backingfile = offsets_filename)
  }
  
  # run negative binomial regression of gene expression on confounders for each
  # gene and for each set of confounders
  for(index in 1:genes_per_task){
    gene_id = parameters$gene_id[index]
    gene_col_index = parameters$gene_col_index[index]
    gene_exp_full = expression_FBM[,match(gene_id, gene_ids)]
    df$gene_exp = gene_exp_full[match(cells_to_keep, cell_names)]
    for(confounder_adjustment_expression in confounder_adjustments_expression){
      cat(sprintf("Running precomputation for gene %s, adjusting for confounders %s...\n", 
                  gene_id, confounder_adjustment_expression))
      # run the negative binomial regression
      formula = as.formula(sprintf("gene_exp ~ %s", confounder_adjustment_expression))
      # note: raw dispersion estimate hardcoded here
      size_raw = 1/(disp_table %>% filter(gene_id == !!gene_id) %>% pull(disp))
      nb_fit = vglm(formula, data = df, family = negbinomial.size(size = size_raw))
      # get fitted values (on natural parameter scale)
      offsets[[confounder_adjustment_expression]][,gene_col_index] = log(nb_fit@fitted.values)
    }
  }
}

# precompute estimated probability of gRNA occurrence in each cell, per gRNA
run_one_precomputation_per_grna = function(task_index, base_dir){
  # set up all analysis parameters
  task = "precomputation_per_grna"
  input_mode = "prepare_task"
  source("analysis/input_file.R", local = TRUE)
  
  # read in the data
  source("analysis/read_data.R", local = TRUE)
  
  # extract confounder adjustments to apply
  confounder_adjustments_grna = resampling_methods %>% 
    filter(!is.na(confounder_adjustment_grna)) %>% 
    pull(confounder_adjustment_grna) %>%
    unique()
  
  # subset data to cells containing at least one guide RNA      
  df_full = confounders_full
  df_full$cell_name = cell_names
  df = df_full %>% filter(guide_count > 0)
  cells_to_keep = df %>% pull(cell_name)
  num_cells_to_keep = length(cells_to_keep)
  
  # file-backed matrix to store output
  grnas_per_task = nrow(parameters)
  num_confounder_adjustments = length(confounder_adjustments_grna)
  probabilities = vector("list", num_confounder_adjustments)
  names(probabilities) = confounder_adjustments_grna
  for(confounder_adjustment_index in 1:length(confounder_adjustments_grna)){
    confounder_adjustment_grna = confounder_adjustments_grna[confounder_adjustment_index]
    probabilities_filename = sprintf("%s/precomp/probabilities_index_%d_confounders_%d", 
                                     base_dir, task_index, confounder_adjustment_index)
    probabilities[[confounder_adjustment_index]] = FBM(nrow = num_cells_to_keep, ncol = grnas_per_task, 
                                                 type = "double", backingfile = probabilities_filename,
                                                 create_bk = TRUE)
  }

  # fit model grna ~ confounders for each grna and each set of confounders
  for(index in 1:grnas_per_task){
    grna_group = parameters$grna_group[index]
    grna_col_index = parameters$grna_col_index[index]
    
    grna_group_indicators_full = phenodata_FBM[,match(grna_group, pd_colnames)]
    df$grna_group_indicator = grna_group_indicators_full[match(cells_to_keep, cell_names)]
    for(confounder_adjustment_grna in confounder_adjustments_grna){
      cat(sprintf("Running precomputation for grna group %s, adjusting for confounders %s...\n", 
                  grna_group, confounder_adjustment_grna))
      formula = as.formula(sprintf("grna_group_indicator ~ %s", confounder_adjustment_grna))
      model_for_grna = glm(formula, family = binomial(), data = df)
      probabilities[[confounder_adjustment_grna]][,grna_col_index] = model_for_grna$fitted.values
    }
  }
}