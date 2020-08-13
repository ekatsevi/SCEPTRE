########################################################
#
# Read in results from parallel jobs and process them 
# into a more easily analyzable form.
#
########################################################

### post-process Gasperini et al results
original_results_raw = suppressWarnings(read_tsv(sprintf("%s/data/raw/CRISPR/GSE120861_all_deg_results.at_scale.txt", base_dir),
                                col_types = "cddddddccccciiciiccc"))
old_pairs = read_excel(path = sprintf("%s/data/raw/CRISPR/Gasperini_TableS2.xlsx", base_dir), sheet = 3)

get_target_site = function(grna_group){
  if(!grepl("_two", grna_group)){
    target_site = grna_group
  } else{
    target_site = strsplit(grna_group, "_")[[1]][1]
  }
}

# extract target sites, high confidence subset
original_results = original_results_raw %>% 
  group_by(gRNA_group) %>% 
  mutate(target_site = get_target_site(unique(gRNA_group))) %>% 
  ungroup() %>%
  left_join(old_pairs %>% 
              select(Target_Site, ENSG, high_confidence_subset) %>%
              mutate(quality_rank_grna = "top_two") %>%
              rename(target_site = Target_Site), 
            by = c("target_site", "ENSG", "quality_rank_grna")) %>%
  rename(grna_group = gRNA_group, gene_id = ENSG, pair_id = pairs4merge, chr = target_gene.chr) %>%
  mutate(TSS = ifelse(strand == "+", target_gene.start, target_gene.stop),
         rejected = ifelse(is.na(high_confidence_subset), FALSE, high_confidence_subset)) %>%
  select(chr, pair_id, rejected,
         gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene, 
         grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type, 
         beta, intercept, fold_change.transcript_remaining, pvalue.raw, pvalue.empirical, pvalue.empirical.adjusted)

# fill in missing empirical p-values
null_ecdf = ecdf(c(0,original_results %>% filter(site_type == "NTC") %>% pull(pvalue.raw)))
original_results = original_results %>%
  mutate(pvalue.empirical = ifelse(is.na(pvalue.empirical), null_ecdf(pvalue.raw), pvalue.empirical))

# write to file
write_tsv(original_results, path = sprintf("%s/results/original_results.tsv", base_dir))

results_files = list.files(sprintf("%s/results/pvalues_per_task", base_dir), full.names = TRUE)

### collate resampling results
resampling_output_files = results_files[grepl("resampling", results_files)]
resampling_output_list = vector("list", length(resampling_output_files))
for(index in 1:length(resampling_output_files)){
  cat(sprintf("Reading file %d out of %d...\n", index, length(resampling_output_files)))
  df = read_tsv(resampling_output_files[[index]], 
                col_types = "dcdcccccc", 
                comment = "#", progress = FALSE)
  df = df %>%
    group_by(gene_id, grna_group) %>%
    spread(significance, value) %>% ungroup()
  resampling_output_list[[index]] = df
}
resampling_output_df = do.call("rbind", resampling_output_list)

resampling_results = resampling_output_df %>% 
  select(gene_id, grna_group, resampling_type, original_zvalue, 
         corrected_pvalue_raw, corrected_pvalue_st, xi, omega, alpha, nu) %>%
  rename(method = resampling_type) %>%
  left_join(original_results %>% 
              select(chr, pair_id, 
                     gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene, 
                     grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type),
            by = c("gene_id", "grna_group")) %>%
  group_by(method, site_type, quality_rank_grna) %>%
  mutate(adjusted_pvalue = ifelse(site_type == "DHS" & quality_rank_grna == "top_two", 
                                  p.adjust(corrected_pvalue_st, "fdr"), NA),
         rejected = ifelse(is.na(adjusted_pvalue), FALSE, adjusted_pvalue <= 0.1)) %>%
  ungroup() %>%
  select(chr, pair_id, rejected,
         gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene, 
         grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type, 
         method, original_zvalue, corrected_pvalue_raw, corrected_pvalue_st, xi, omega, alpha, nu)

# write to file
write_tsv(resampling_results, path = sprintf("%s/results/resampling_results.tsv", base_dir))

### collate likelihood results

likelihood_output_files = results_files[grepl("likelihood", results_files)]
likelihood_output_list = vector("list", length(likelihood_output_files))
for(index in 1:length(likelihood_output_files)){
  cat(sprintf("Reading file %d out of %d...\n", index, length(likelihood_output_files)))
  df = read_tsv(likelihood_output_files[[index]], 
                col_types = "dcdccccc", 
                comment = "#", progress = FALSE)
  df = df %>%
    # filter(!is.na(value)) %>%
    group_by(gene_id, grna_group) %>%
    spread(significance, value) %>% ungroup()
  likelihood_output_list[[index]] = df
}
likelihood_output_df = do.call("rbind", likelihood_output_list)

likelihood_results = likelihood_output_df %>% 
  mutate(method = ifelse(method_idx %% 2 == 0, 
                         sprintf("raw_%d", round(method_idx/2)), 
                         sprintf("shrunk_%d", round((method_idx + 1)/2)))) %>%
  left_join(original_results %>% 
              select(chr, pair_id, 
                     gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene, 
                     grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type),
            by = c("gene_id", "grna_group")) %>%
  select(chr, pair_id, 
         gene_id, gene_short_name, target_gene.start, target_gene.stop, TSS, outlier_gene, 
         grna_group, quality_rank_grna, target_site, target_site.start, target_site.stop, site_type, 
         method, dispersion_estimate, confounder_adjustment, beta, pvalue, zvalue)

# write to file
write_tsv(likelihood_results, path = sprintf("%s/results/likelihood_results.tsv", base_dir))