####################################################
#
# Reproduce Gasperini et al p-values by re-running
# negative binomial GLM on expression/gRNA data.
#
####################################################

num_pairs_to_reproduce = 10
original_pvalues = numeric(num_pairs_to_reproduce)
reproduced_pvalues = numeric(num_pairs_to_reproduce)

# iterate over pairs
for(pair_idx in 1:num_pairs_to_reproduce){
  ###########################################
  # Set up for negative binomial GLM analysis
  ###########################################
  
  # randomly choose a gRNA/gene pair to analyze
  random_row = sample(1:nrow(all_deg_results), 1)
  grna_group = all_deg_results$gRNA_group[random_row]
  gene_id = all_deg_results$ENSG[random_row]
  cat(sprintf("Reproducing analysis for gRNA %s and gene %s.\n", grna_group, gene_id))
  
  # extract confounders
  confounder_names = c("guide_count", "prep_batch", "percent.mito")
  confounders_full = phenodata_FBM[,match(confounder_names, pd_colnames)]
  colnames(confounders_full) = confounder_names
  confounders_full = as.data.frame(confounders_full)
  
  # extract size factors
  Size_Factors_full = phenodata_FBM[,match("Size_Factor", pd_colnames)]
  
  # extra gRNA indicators for chosen gRNA
  grna_indicators_full = phenodata_FBM[,match(grna_group, pd_colnames)]
  
  # extract expression of chosen gene
  gene_exp_full = expression_FBM[,match(gene_id, gene_ids)]
  
  # keep cells that are either in reference or express the gRNA
  cells_to_keep = which((cell_names %in% reference_cells | (grna_indicators_full == 1)) & 
                          (confounders_full$guide_count > 0))
  
  # put all information in a data frame
  df = confounders_full[cells_to_keep,] 
  df$grna = grna_indicators_full[cells_to_keep]
  df$gene_exp = gene_exp_full[cells_to_keep]
  df$Size_Factor = Size_Factors_full[cells_to_keep]
  
  ###########################################
  # Carry out negative binomial GLM analysis
  ###########################################
  
  # extract the raw mean and dispersion of gene expression
  mean_estimate_raw = disp_table %>% filter(gene_id == !!gene_id) %>% pull(mu)
  disp_estimate_raw = disp_table %>% filter(gene_id == !!gene_id) %>% pull(disp)
  
  # calculate the shrunken dispersion estimate based on the DESeq2 formula (6)
  disp_estimate_shrunk = disp_coeffs[1] + disp_coeffs[2] / mean_estimate_raw
  
  # normalize expression by the size factor (I think the size factor would be better accounted for using offsets)
  df = df %>% mutate(normalized_gene_exp = round(gene_exp/Size_Factor))
  
  # fit the negative binomial model with grna
  full_model_fit <- VGAM::vglm(normalized_gene_exp ~ grna + guide_count + prep_batch + percent.mito, 
                               epsilon=1e-1, 
                               family=negbinomial.size(size=1/disp_estimate_shrunk), 
                               data = df)
  
  # fit the negative binomial model without grna
  reduced_model_fit <- VGAM::vglm(normalized_gene_exp ~ guide_count + prep_batch + percent.mito, 
                                  epsilon=1e-1, 
                                  family=negbinomial.size(size=1/disp_estimate_shrunk), 
                                  data = df)
  
  # obtain a p-value by comparing these two models via the likelihood ratio test
  reproduced_pvalues[pair_idx] = lrtest(full_model_fit, reduced_model_fit)@Body["Pr(>Chisq)"][2,]
  
  # compare to p-value downloaded from GEO
  original_pvalues[pair_idx] = all_deg_results$pvalue.raw[random_row]
}

print(cbind(original_pvalues, reproduced_pvalues))