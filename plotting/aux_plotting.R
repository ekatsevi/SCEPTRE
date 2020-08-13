################################
#
# Auxiliary plotting functions.
#
################################

# transformation to produce QQ plots on log scale
# code borrowed from https://gist.github.com/JoFrhwld/2266961
revlog_trans <- function(base = exp(1)){
  ## Define the desired transformation.
  trans <- function(x){
    -log(x, base)
  }
  ## Define the reverse of the desired transformation
  inv <- function(x){
    base^(-x)
  }
  ## Creates the transformation
  trans_new(paste("revlog-", base, sep = ""),
            trans, ## The transformation function (can be defined using anonymous functions)
            inv,  ## The reverse of the transformation
            log_breaks(base = base), ## default way to define the scale breaks
            domain = c(1e-100, Inf) ## The domain over which the transformation is valued
  )
}

# make the skew-t plot
plot_skew_t = function(gene_to_plot, grna_to_plot){
  B = 500
  experiment_index = parameters %>% filter(gene_id == gene_to_plot, grna_group == grna_to_plot) %>% pull(experiment_index)
  num_pairs = parameters %>% filter(experiment_index == !!experiment_index) %>% nrow()
  index = parameters %>% filter(experiment_index == !!experiment_index) %>% 
    mutate(index = row_number()) %>% filter(gene_id == gene_to_plot, grna_group == grna_to_plot) %>%
    pull(index)
  zvalues = FBM(nrow = B, ncol = num_pairs, type = "double", 
                backingfile = sprintf("%s/results/resampled_zvalues/zvalues_index_%d_method_1", 
                                      base_dir, experiment_index),
                create_bk = FALSE)
  resampled_zvalues = zvalues[,index]
  original_zvalue = resampling_results %>% 
    filter(method == "conditional_randomization",
           gene_id == gene_to_plot, grna_group == grna_to_plot) %>% pull(original_zvalue)
  pvalue_new = resampling_results %>% 
    filter(method == "conditional_randomization",
      gene_id == gene_to_plot, grna_group == grna_to_plot) %>% pull(corrected_pvalue_st)
  pvalue_old = likelihood_results %>% 
    filter(method == "raw_2", gene_id == gene_to_plot, grna_group == grna_to_plot) %>% pull(pvalue)
  plot_title = sprintf("Neg Binom p-val. = %0.0e\nSCEPTRE p-val. = %0.0e", pvalue_old, pvalue_new)
  dp = conditional_randomization_results %>% 
    filter(gene_id == gene_to_plot, grna_group == grna_to_plot) %>%
    select(xi, omega, alpha, nu) %>%
    as.numeric()
  z = seq(-4, 4, length.out = 1000)
  df_curves = tibble(z = z, fitted = dst(x = z, dp = dp), gaussian = dnorm(z)) %>%
    gather("curve", "y", fitted, gaussian) %>%
    mutate(curve = factor(curve, levels = c("fitted","gaussian"), labels = c("Conditional\nrandomization","Negative\nbinomial")))
  df_ribbon = df_curves %>% 
    filter(z <= original_zvalue, curve == "Conditional\nrandomization") %>% 
    mutate(lower = 0, upper = y) %>%
    select(z, lower, upper)
  p = ggplot() + 
    geom_histogram(aes(x = z, y = ..density..), 
                   data = tibble(z = resampled_zvalues), 
                   boundary = 0, colour = "black", fill = "lightgray", binwidth = 0.5) + 
    geom_line(aes(x = z, y = y, group = curve, colour = curve, linetype = curve), data = df_curves) + 
    geom_vline(xintercept = original_zvalue, colour = "firebrick3", linetype = "solid") + 
    geom_ribbon(aes(x = z, ymin = lower, ymax = upper), fill = "darkorchid2", alpha = 0.5, data = df_ribbon) +
    scale_colour_manual(values = c("darkorchid2", "black"), name = "Null distribution") + 
    scale_linetype_manual(values = c("solid", "dashed"), name = "Null distribution") + 
    scale_y_continuous(expand = c(0,0)) + 
    ggtitle(plot_title) +
    # ggtitle("Calibrating negative binomial z-value using conditional permutation") + 
    xlab("Negative binomial z-value") + theme_bw() +
    theme(legend.position = c(0.85,0.8), 
          legend.background = element_rect(fill = "transparent", colour = NA),
          plot.title = element_text(hjust = 0.5, size = 11),
          panel.grid = element_blank(), 
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  return(p)
}