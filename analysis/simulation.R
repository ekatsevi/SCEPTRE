##################################################
#
# Proof-of-concept simulation to show that SCEPTRE 
# gives calibrated p-values despite expression 
# model misspecification.
#
##################################################

# fixed simulation parameters
num_cells = 1000
grna_avg_prob = 0.05
mean_expression = 5
size = 1
reps = 500
B = 250

# variable simulation parameters
parameters = rbind(
  list(
    estimated_size = c(0.2, 1, 5),
    method = c("NB", "MP", "CR"),
    rep = 1:reps,
    size = 1,
    grna_effect = 0,
    zero_inflation = 0
  ) %>%
    cross_df(),
  list(
    estimated_size = c(1),
    method = c("NB", "MP", "CR"),
    rep = 1:reps,
    size = 1,
    grna_effect = 0,
    zero_inflation = 0.25
  ) %>%
    cross_df()
) %>%
  mutate(row_idx = row_number())

set.seed(1234)

results = tibble(pvalue = rep(NA, nrow(parameters)), zvalue = rep(NA, nrow(parameters)))

unique_simulation_params = parameters %>% select(size, grna_effect, zero_inflation) %>% unique()
for (sim_par_idx in 1:nrow(unique_simulation_params)) {
  size = unique_simulation_params$size[sim_par_idx]
  grna_effect = unique_simulation_params$grna_effect[sim_par_idx]
  zero_inflation = unique_simulation_params$zero_inflation[sim_par_idx]
  parameters_sim = parameters %>%
    filter(size == !!size,
           grna_effect == !!grna_effect,
           zero_inflation == !!zero_inflation)
  reps_sim = parameters_sim %>% pull(rep) %>% unique()
  for (rep in reps_sim) {
    cat(sprintf("Working on rep %d/%d...\n", rep, max(reps_sim)))
    # generate data
    Z = rnorm(num_cells)
    pi_true =  1 / (1 + exp(-(logit(grna_avg_prob) + 4*Z)))
    X = rbinom(n = num_cells, size = 1, prob = pi_true)
    Y = rnbinom(n = num_cells,
                size = size,
                mu = exp(log(mean_expression) + 4*Z - grna_effect * X))
    Y = Y * rbinom(n = num_cells,
                   size = 1,
                   prob = 1 - zero_inflation)
    
    methods_sim = parameters_sim %>% filter(rep == !!rep) %>% select(estimated_size, method, row_idx)
    
    # negative binomial
    estimated_sizes = methods_sim %>% filter(method == "NB") %>% pull(estimated_size)
    for (estimated_size in estimated_sizes) {
      # negative binomial regression
      cat(
        sprintf(
          "Running negative binomial regression for estimated size %0.1f...\n",
          estimated_size
        )
      )
      nb_fit = vglm(Y ~ X + Z,
                    family = negbinomial.size(size = estimated_size),
                    data = tibble(X, Y, Z))
      row_idx = methods_sim %>% filter(method == "NB", estimated_size == !!estimated_size) %>% pull(row_idx)
      results[[row_idx, "zvalue"]] = coefficients(summary(nb_fit))["X", "z value"]
      results[[row_idx, "pvalue"]] = pnorm(results[[row_idx, "zvalue"]])
    }
    
    # resampling methods
    estimated_sizes = methods_sim %>% filter(method != "NB") %>% pull(estimated_size) %>% unique()
    for (estimated_size in estimated_sizes) {
      cat(
        sprintf(
          "Working on resampling methods for estimated size %0.1f...\n",
          estimated_size
        )
      )
      # distilled NB test statistic
      nb_loco_fit = vglm(Y ~ Z,
                         family = negbinomial.size(size = estimated_size),
                         data = tibble(Y, Z))
      offsets = log(nb_loco_fit@fitted.values)
      nb_distilled_fit = vglm(
        Y ~ 1,
        offset = offsets[X == 1],
        family = negbinomial.size(size = estimated_size),
        data = tibble(Y = Y[X == 1])
      )
      row_idx = methods_sim %>% filter(method != "NB", estimated_size == !!estimated_size) %>% pull(row_idx)
      results[row_idx, "zvalue"] = coefficients(summary(nb_distilled_fit))["(Intercept)", "z value"]
      
      # logistic regression for conditional permutation
      log_reg_fit = glm(X ~ Z, family = "binomial", data = tibble(X, Z))
      pi_hat = log_reg_fit$fitted.values
      
      # resampling
      for (method in methods_sim %>% filter(method %in% c("CP", "MP", "CR"),
                                            estimated_size == !!estimated_size) %>% pull(method)) {
        cat(sprintf("Running %s...\n", method))
        resampled_zvalues = numeric(B)
        for (b in 1:B) {
          if (method == "MP") {
            cells_with_grna = sample.int(num_cells, size = sum(X))
          }
          if(method == "CR"){
            cells_with_grna = as.logical(rbinom(n = num_cells, size = 1, prob = pi_hat))
          }
          nb_distilled_fit = vglm(
            Y ~ 1,
            offset = offsets[cells_with_grna],
            family = negbinomial.size(size = estimated_size),
            data = tibble(Y = Y[cells_with_grna])
          )
          resampled_zvalues[b] = coefficients(summary(nb_distilled_fit))["(Intercept)", "z value"]
        }
        dp = selm(y ~ 1,
                  family = "ST",
                  data = tibble(y = resampled_zvalues))@param$dp
        row_idx = methods_sim %>% filter(method == !!method, estimated_size == !!estimated_size) %>% pull(row_idx)
        results[[row_idx, "pvalue"]] = pst(x = results[[row_idx, "zvalue"]], dp = dp)
      }
    }
  }
}

# write to file
df = cbind(parameters, results) %>% as_tibble()
write_tsv(df, sprintf("%s/results/simulation_results.tsv", base_dir))