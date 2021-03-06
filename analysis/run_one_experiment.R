#####################################################
#
# This is the primary function that carries out the 
# association analyses for several (gene,gRNA) pairs.
#
#####################################################

run_one_experiment = function(task_index, base_dir){
  # set up all analysis parameters
  task = "experiment"
  input_mode = "prepare_task"
  source("analysis/input_file.R", local = TRUE)
  
  # set up array to hold output from resampling methods
  num_pairs = nrow(parameters)
  num_resampling_methods = nrow(resampling_methods)
  resampling_output = array(NA,
                            dim = c(num_resampling_methods, num_pairs, 7),
                            dimnames = list(NULL, NULL,
                                            c("original_zvalue",
                                              "xi", "omega", "alpha", "nu",
                                              "corrected_pvalue_raw", "corrected_pvalue_st")))
  
  # set up file-backed matrices to hold resampled z-values from resampling methods
  zvalue_FBMs = vector("list", num_resampling_methods)
  for(resampling_method_index in 1:num_resampling_methods){
    zvalues_filename = sprintf("%s/results/resampled_zvalues/zvalues_index_%d_method_%d", 
                               base_dir, task_index, resampling_method_index)
    zvalue_FBMs[[resampling_method_index]] = FBM(nrow = B, 
                                                 ncol = num_pairs, 
                                                 type = "double",
                                                 backingfile = zvalues_filename, 
                                                 create_bk = TRUE)
  }
  
  # set up array to hold output from likelihood-based methods
  num_likelihood_methods = nrow(likelihood_methods)
  likelihood_output = array(NA,
                            dim = c(num_likelihood_methods, num_pairs, 3),
                            dimnames = list(NULL, NULL, c("beta", "pvalue", "zvalue")))
      
  # read in the data
  source("analysis/read_data.R", local = TRUE)

  # subset data to cells containing at least one guide RNA
  df_full = confounders_full
  df_full$cell_name = cell_names
  df = df_full %>% filter(guide_count > 0)
  cells_to_keep = df %>% pull(cell_name)
  num_cells_to_keep = length(cells_to_keep)
  
  # iterate over gene/gRNA pairs
  for(index in 1:num_pairs){
    gene_id =            parameters$gene_id[index]
    grna_group =         parameters$grna_group[index]
    grna_precomp_index = parameters$grna_precomp_index[index]
    gene_precomp_index = parameters$gene_precomp_index[index]
    grna_col_index =     parameters$grna_col_index[index]
    gene_col_index =     parameters$gene_col_index[index]
    grnas_per_task =     parameters$grnas_per_task[index]
    genes_per_task =     parameters$genes_per_task[index]

    cat(sprintf("Working on pair number %d: gRNA %s and gene %s...\n", index, grna_group, gene_id))

    # read gene expression
    gene_exp_full = expression_FBM[,match(gene_id, gene_ids)]
    df$gene_exp = gene_exp_full[match(cells_to_keep, cell_names)]
    
    # # read grna group indicators    
    grna_group_indicators_full = phenodata_FBM[,match(grna_group, pd_colnames)]
    df$grna_group_indicator = grna_group_indicators_full[match(cells_to_keep, cell_names)]
    
    # read precomputed probabilities for resampling
    confounder_adjustments_grna = resampling_methods %>%
      filter(!is.na(confounder_adjustment_grna)) %>%
      pull(confounder_adjustment_grna) %>%
      unique()
    probabilities_all = matrix(0, length(confounder_adjustments_grna), num_cells_to_keep)
    rownames(probabilities_all) = confounder_adjustments_grna
    for(confounder_adjustment_index in 1:length(confounder_adjustments_grna)){
      confounder_adjustment_grna = confounder_adjustments_grna[confounder_adjustment_index]
      probabilities_filename = sprintf("%s/precomp/probabilities_index_%d_confounders_%d",
                                       base_dir, grna_precomp_index, confounder_adjustment_index)
      probabilities_FBM = FBM(nrow = num_cells_to_keep, ncol = grnas_per_task,
                          type = "double", backingfile = probabilities_filename,
                          create_bk = FALSE)
      probabilities_all[confounder_adjustment_grna,] = probabilities_FBM[,grna_col_index]
    }

    # read precomputed offsets based on confounder_adjustment_expression
    confounder_adjustments_expression = resampling_methods %>% 
      pull(confounder_adjustment_expression) %>% unique()
    
    offsets = matrix(0, length(confounder_adjustments_expression), num_cells_to_keep)
    rownames(offsets) = confounder_adjustments_expression
    for(confounder_adjustment_index in 1:length(confounder_adjustments_expression)){
      confounder_adjustment_expression = confounder_adjustments_expression[confounder_adjustment_index]
      offsets_filename = sprintf("%s/precomp/offsets_index_%d_confounders_%d", 
                                 base_dir, gene_precomp_index, confounder_adjustment_index)
      offsets_FBM = FBM(nrow = num_cells_to_keep, ncol = genes_per_task, type = "double",
                        backingfile = offsets_filename, create_bk = FALSE)
      offsets[confounder_adjustment_expression,] = offsets_FBM[,gene_col_index]
    }

    # sizes
    size_raw    = 1/(disp_table %>% filter(gene_id == !!gene_id) %>% pull(disp))
    size_shrunk = 1/(disp_coeffs[1] + disp_coeffs[2] / mean(df$gene_exp))
    
    # run resampling methods
    if(num_resampling_methods > 0){
      for(resampling_method_idx in 1:num_resampling_methods){
        resampling_type = resampling_methods$resampling_type[resampling_method_idx]
        test_statistic = resampling_methods$test_statistic[resampling_method_idx]
        dispersion_estimate = resampling_methods$dispersion_estimate[resampling_method_idx]
        confounder_adjustment_grna = resampling_methods$confounder_adjustment_grna[resampling_method_idx]
        confounder_adjustment_expression = resampling_methods$confounder_adjustment_expression[resampling_method_idx]
        cat(sprintf("Running %s test with %s test statistics, %s dispersion estimate, and confounders %s and %s...\n", 
                    resampling_type, test_statistic, dispersion_estimate, confounder_adjustment_grna, 
                    confounder_adjustment_expression))
        
        # find size parameter for negative binomial regression
        if(dispersion_estimate == "size_raw"){
          size = size_raw
        } else if(dispersion_estimate == "size_shrunk"){
          size = size_shrunk
        } else{
          print("Invalid dispersion estimate specified!")
          stopifnot(FALSE)
        }
        
        ### COMPUTE ORIGINAL Z-VALUE ###
        cat(sprintf("Computing original z-value...\n"))
        if(test_statistic == "LOCO_NB"){
          # run NB regression of gene expression on gRNA
          formula = as.formula(sprintf("gene_exp ~ %s", confounder_adjustment_expression))
          nb_fit = vglm(formula, data = df, family = negbinomial.size(size = size_raw))
          NB_fit = vglm(gene_exp ~ 1, 
                        offset = offsets[confounder_adjustment_expression, df$grna_group_indicator == 1],
                        family=negbinomial.size(size = size), 
                        data = df[df$grna_group_indicator == 1,])
          original_zvalue = coefficients(summary(NB_fit))["(Intercept)", "z value"]
        }
        
        if(test_statistic == "full_NB"){
          # run NB regression of gene expression on gRNA and confounders
          formula = as.formula(sprintf("gene_exp ~ grna_group_indicator + %s", confounder_adjustment_expression))
          NB_fit = vglm(formula, family=negbinomial.size(size = size), data = df)
          original_zvalue = coefficients(summary(NB_fit))["grna_group_indicator", "z value"]
        }
        
        ### COMPUTE RESAMPLED Z-VALUES ###
        # extract probabilities for conditional resampling
        if(resampling_type == "conditional_randomization"){
          probabilities = probabilities_all[confounder_adjustment_grna,]
        }
        
        # set seed for reproducibility
        set.seed(1234) 

        resampled_zvalues = numeric(B)
        for(b in 1:B){
          if(b %% 10 == 0){
            cat(sprintf("Working on resample number %d...\n", b))
          }
          # resample 
          if(resampling_type == "marginal_permutation"){
            grna_indicators_null = numeric(num_cells_to_keep)
            grna_indicators_null[sample.int(num_cells_to_keep, size = sum(df$grna_group_indicator))] = 1
          } else if(resampling_type == "conditional_randomization"){
            grna_indicators_null = rbinom(n = num_cells_to_keep, size = 1, prob = probabilities)
          }
          # test
          if(test_statistic %in% c("LOCO_NB")){
            NB_fit = vglm(gene_exp ~ 1, 
                          offset = offsets[confounder_adjustment_expression, grna_indicators_null == 1],
                          family=negbinomial.size(size = size), 
                          data = df[grna_indicators_null == 1,])
            resampled_zvalues[b] = coefficients(summary(NB_fit))["(Intercept)", "z value"]
          }
          if(test_statistic == "full_NB"){
            formula = as.formula(sprintf("gene_exp ~ grna_indicators_null + %s", confounder_adjustment_expression))
            NB_fit = vglm(formula, family=negbinomial.size(size = size), data = cbind(df, grna_indicators_null))
            resampled_zvalues[b] = coefficients(summary(NB_fit))["grna_indicators_null", "z value"]      
          }
        }

        ### FIT SKEW-T DISTRIBUTION ###
        
        dp = tryCatch({
          selm(y ~ 1, family = "ST", data = tibble(y = resampled_zvalues))@param$dp
        }, error = function(error){
          cat(sprintf("skew-t fit did not work, returning NA for now...\n"))
          return(NA)
        })
        
        ### COMPUTE CORRECTED P-VALUE ###
        corrected_pvalue_raw = mean(c(-Inf, resampled_zvalues) <= original_zvalue)
        if(!any(is.na(dp))){
          corrected_pvalue_st = pst(x = original_zvalue, dp = dp) 
        } else{
          corrected_pvalue_st = corrected_pvalue_raw
        }

        ### RECORD THE RESULTS
        resampling_output[resampling_method_idx, index, "original_zvalue"] = original_zvalue
        zvalue_FBMs[[resampling_method_idx]][,index] = resampled_zvalues
        resampling_output[resampling_method_idx, index, c("xi", "omega", "alpha", "nu")] = dp
        resampling_output[resampling_method_idx, index, "corrected_pvalue_st"] = corrected_pvalue_st
        resampling_output[resampling_method_idx, index, "corrected_pvalue_raw"] = corrected_pvalue_raw
      }
      
      if(index %% 100 == 0 | index == num_pairs){
        if(index == num_pairs){
          cat(sprintf("Writing final output to file!\n")) 
        } else{
          cat(sprintf("Writing partial output to file!\n")) 
        }
        
        # massage resampling results into data frame
        df_resampling_output = melt(resampling_output)
        names(df_resampling_output) = c("method_idx", "pair_idx", "significance", "value")
        df_resampling_output = as_tibble(df_resampling_output)
        df_resampling_output = df_resampling_output %>% 
          mutate(gene_id = parameters$gene_id[pair_idx], grna_group = parameters$grna_group[pair_idx]) %>%
          mutate(resampling_type = resampling_methods$resampling_type[method_idx], 
                 test_statistic = resampling_methods$test_statistic[method_idx],
                 dispersion_estimate = resampling_methods$dispersion_estimate[method_idx],
                 confounder_adjustment = resampling_methods$confounder_adjustment_grna[method_idx]) %>%
          select(-pair_idx)
        
        resampling_output_filename = sprintf("%s/results/pvalues_per_task/resampling_output_%d.tsv", base_dir, task_index)
        write_tsv(df_resampling_output, path = resampling_output_filename)
      }
    }
    
    # run likelihood-based methods
    if(num_likelihood_methods > 0){
      for(likelihood_method_idx in 1:num_likelihood_methods){
        test = likelihood_methods$test[likelihood_method_idx]
        dispersion_estimate = likelihood_methods$dispersion_estimate[likelihood_method_idx]
        confounder_adjustment = likelihood_methods$confounder_adjustment[likelihood_method_idx]
        cat(sprintf("Running negative binomal test with %s statistics and %s dispersion estimate...\n", 
                    test, dispersion_estimate))

        # find size parameter for negative binomial regression
        if(dispersion_estimate == "size_raw"){
          size = size_raw
        } else if(dispersion_estimate == "size_shrunk"){
          size = size_shrunk
        } else{
          print("Invalid dispersion estimate specified!")
          stopifnot(FALSE)
        }
        
        formula = as.formula(sprintf("gene_exp ~ grna_group_indicator + %s", confounder_adjustment))
        NB_fit = vglm(formula, family=negbinomial.size(size = size), data = df)
        NB_fit_orig = vglm(formula, family=negbinomial.size(size = size), 
                           data = df %>% filter(grna_group_indicator == 1 | cell_name %in% reference_cells))
        beta_coef = unname(NB_fit@coefficients["grna_group_indicator"])
        if(test == "z_test"){
          zvalue = tryCatch({
            coefficients(summary(NB_fit))["grna_group_indicator", "z value"]
          }, error = function(error){
            cat(sprintf("NB regression did not work, returning NA for now...\n"))
            return(NA)
          })
          pvalue = pnorm(zvalue)
        }
        if(test == "lr_test"){
          formula = as.formula(sprintf("gene_exp ~ %s", confounder_adjustment))
          reduced_fit = vglm(formula, family=negbinomial.size(size = size), data = df)
          reduced_fit_orig = vglm(formula, family=negbinomial.size(size = size), 
                                  data = df %>% filter(grna_group_indicator == 1 | cell_name %in% reference_cells))
          pvalue = lrtest(NB_fit, reduced_fit)@Body["Pr(>Chisq)"][2,]
          pvalue_orig = lrtest(NB_fit_orig, reduced_fit_orig)@Body["Pr(>Chisq)"][2,]
          zvalue = NA
        }
        
        # record results
        likelihood_output[likelihood_method_idx, index, "beta"] = beta_coef
        likelihood_output[likelihood_method_idx, index, "pvalue"] = pvalue
        likelihood_output[likelihood_method_idx, index, "zvalue"] = zvalue
      }
      
      if(index %% 100 == 0 | index == num_pairs){
        if(index == num_pairs){
          cat(sprintf("Writing final output to file!\n")) 
        } else{
          cat(sprintf("Writing partial output to file!\n")) 
        }
        # write to file
        df_likelihood_output = melt(likelihood_output)
        names(df_likelihood_output) = c("method_idx", "pair_idx", "significance", "value")
        df_likelihood_output = as_tibble(df_likelihood_output)
        df_likelihood_output = df_likelihood_output %>%
          mutate(gene_id = parameters$gene_id[pair_idx], grna_group = parameters$grna_group[pair_idx]) %>%
          mutate(test = likelihood_methods$test[method_idx], 
                 dispersion_estimate = likelihood_methods$dispersion_estimate[method_idx],
                 confounder_adjustment = likelihood_methods$confounder_adjustment[method_idx]) %>%
          select(-pair_idx)
        
        likelihood_output_filename = sprintf("%s/results/pvalues_per_task/likelihood_output_%d.tsv", 
                                             base_dir, task_index)
        write_tsv(df_likelihood_output, path = likelihood_output_filename)
      }
    }
  }
}