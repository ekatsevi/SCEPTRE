#######################################################
#
# Reproduce Figure 4 from Katsevich and Roeder (2020).
#
#######################################################

# a: DHS and NTC p-values
ci = 0.95
subsampling_factor = 250
p = 
  resampling_results %>%
  filter(site_type %in% c("DHS", "NTC"), method == "conditional_randomization") %>%
  rename(pvalue = corrected_pvalue_st) %>%
  mutate(site_type = factor(site_type, levels = c("DHS", "NTC"), labels = c("Candidate enhancer", "Negative control"))) %>%
  group_by(site_type) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%
  mutate(clower = ifelse(site_type == "Candidate enhancer", NA, clower), 
         cupper = ifelse(site_type == "Candidate enhancer", NA, cupper)) %>%
  mutate(pvalue = ifelse(pvalue < 1e-20, 0, pvalue)) %>%
  filter(-log10(expected) > 2 | (site_type == "Candidate enhancer" & row_number() %% 20 == 0) | row_number() %% subsampling_factor == 0) %>%
  ggplot(aes(x = expected, y = pvalue, group = site_type, ymin = clower, ymax = cupper)) + 
  geom_point(aes(colour = site_type), size = 3, alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("firebrick3", "slategray"), name = "Perturbation target") + 
  scale_x_continuous(trans = revlog_trans(base = 10)) + scale_y_continuous(trans = revlog_trans(base = 10)) +
  xlab(expression(paste("Expected null p-value for gene-enhancer pair"))) + 
  ylab(expression(paste("SCEPTRE p-value"))) + 
  theme_bw() + theme(legend.position = c(0.25,0.8),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure4a.png", base_dir), plot = p, device = "png", 
       width = 5, height = 3)

# b: positive controls
combined_results = rbind(
  original_results %>% 
    mutate(pvalue.empirical = ifelse(beta > 0, 1, pvalue.empirical)) %>%
    select(gene_id, target_site, quality_rank_grna, site_type, pvalue.empirical) %>% 
    rename(pvalue = pvalue.empirical) %>%
    mutate(method = "Original"),
  resampling_results %>%
    filter(method == "conditional_randomization") %>%
    select(gene_id, target_site, quality_rank_grna, site_type, corrected_pvalue_st) %>% 
    rename(pvalue = corrected_pvalue_st) %>%
    mutate(method = "SCEPTRE")
)

df = combined_results %>% 
  filter(site_type %in% c("selfTSS", "positive_ctrl")) %>%
  mutate(pvalue = ifelse(pvalue == 0, 1e-17, pvalue),
         method = factor(method, levels = c("SCEPTRE", "Original")),
         site_type = factor(site_type, 
                            levels = c("selfTSS", "positive_ctrl"), 
                            labels = c("Positive controls (TSS)", "Positive controls (enhancer)"))) %>%
  group_by(site_type) %>% mutate(num_pairs = n()) %>% ungroup()

p = df %>%  
  ggplot(aes(x = site_type, y = pvalue, fill = method, colour = method)) + 
  geom_violin(scale = "width") +
  scale_fill_manual(values = c("blue", "gray")) + 
  scale_colour_manual(values = c("blue", "gray")) + 
  scale_y_log10() + 
  xlab("Perturbation target") + ylab("p-value for positive control pair") + 
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = c(0.85,0.2),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure4b.pdf", base_dir), plot = p, device = "pdf",
       width = 4.5, height = 3)

### c: old versus new DHS p-values ###
annotations = tibble(gene_short_name = c("B3GNT2", "PTPN1", "TOP1",  "AGFG1", "EIF1"),
                     eQTL = c(2.7e-26, NA, NA, 5.2e-8,NA),
                     eRNA = c(NA, 2e-18, 6.6e-5,NA, 1.2e-6))

promising_new_rejections = resampling_results %>%
  filter(gene_short_name %in% c("B3GNT2", "PTPN1", "TOP1",  "AGFG1", "EIF1"), rejected, method == "conditional_randomization") %>%
  inner_join(original_results %>% filter(pvalue.empirical.adjusted > 0.1) %>% 
               select(gene_id, grna_group, pvalue.empirical), by = c("gene_id", "grna_group")) %>%
  left_join(annotations, by = "gene_short_name") %>% 
  arrange(desc(quality_rank_grna), corrected_pvalue_st) %>%
  mutate(pair_number = row_number()) %>%
  select(pair_number, pair_id, gene_short_name, target_site, corrected_pvalue_st, 
         pvalue.empirical, eQTL, eRNA, quality_rank_grna)

thresh_new = resampling_results %>% filter(rejected, method == "conditional_randomization") %>% summarise(max(corrected_pvalue_st)) %>% pull()
thresh_old = original_results %>% filter(rejected) %>% summarise(max(pvalue.empirical)) %>% pull()
p = resampling_results %>% 
  filter(method == "conditional_randomization") %>%
  rename(new_pvalue = corrected_pvalue_st, new_rejected = rejected) %>%
  select(pair_id, gene_id, grna_group, new_pvalue, site_type, quality_rank_grna, new_rejected) %>%
  left_join(original_results %>% 
              rename(old_pvalue = pvalue.empirical, old_rejected = rejected) %>%
              select(gene_id, grna_group, old_pvalue, old_rejected, beta, outlier_gene), by = c("gene_id", "grna_group")) %>%
  filter(site_type == "DHS", quality_rank_grna == "top_two") %>%
  select(-c("site_type", "quality_rank_grna")) %>%
  mutate(reason_not_rejected = ifelse(outlier_gene, 
                                      "Outlier gene", 
                                      ifelse(beta > 0, 
                                             "Positive effect", 
                                             ifelse(!old_rejected & old_pvalue <= thresh_old, "Low confidence", "none"))),
         reason_not_rejected = factor(reason_not_rejected, levels = c("none", "Positive effect", "Outlier gene", "Low confidence")),
         rejected_by = ifelse(old_rejected, 
                              ifelse(new_rejected, "Both methods (368)", "Original only (102)"), 
                              ifelse(new_rejected, "SCEPTRE only (217)", "Neither method")),
         rejected_by = factor(rejected_by, levels = c("Both methods (368)", "SCEPTRE only (217)",
                                                      "Original only (102)", "Neither method"))) %>%
  left_join(promising_new_rejections %>% 
              select(pair_number, pair_id) %>% 
              mutate(pair_number = as.character(pair_number)), by = "pair_id") %>%
  mutate(pair_number = ifelse(is.na(pair_number), "", pair_number)) %>%
  arrange(old_rejected, new_rejected) %>%
  filter(old_pvalue < 1e-2 |
           new_pvalue < 1e-2 |
           reason_not_rejected != "none" |
           row_number() %% 10 == 0) %>%
  ggplot(aes(x = old_pvalue, y = new_pvalue, 
             colour = rejected_by, 
             shape = reason_not_rejected,
             size = reason_not_rejected,
             label = pair_number)) + geom_point() + 
  geom_text_repel(colour = "black", force = 0.8, point.padding = 0.1, box.padding = 0.4, size = 4, min.segment.length = 0) + 
  geom_vline(xintercept = thresh_old, linetype = "dashed") + geom_hline(yintercept = thresh_new, linetype = "dashed") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond"), 
                     breaks = c("Outlier gene", "Positive effect", "Low confidence"),
                     name = "Reason for original\nnon-discovery") + 
  scale_colour_manual(values = c("purple", "blue", "red", "gray40"), name = "Discovered by") + 
  scale_size_manual(values = c(1.5,1.5,1.5,2.5), 
                    breaks = c("Outlier gene", "Positive effect", "Low confidence"),
                    name = "Reason for original\nnon-discovery") + 
  guides(colour = guide_legend(order = 1)) +
  guides(shape = guide_legend(order = 2)) +
  guides(size = guide_legend(order = 2)) +
  scale_x_log10() + scale_y_log10() +
  xlab("Original empirical p-value of gene-enhancer pair") + 
  ylab("SCEPTRE p-value of gene-enhancer pair") + #ggtitle("Discovery of enhancer-gene relationships") + 
  theme_bw() + theme(
                     legend.position = c(0.85, 0.325),
                     legend.title=element_text(size=9), 
                     legend.text=element_text(size=8),
                     legend.key = element_rect(colour = "transparent", fill = "transparent"),
                     legend.margin = margin(-0.25,0,0,0, unit="cm"),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure4c.pdf", base_dir), plot = p, device = "pdf",
       width = 5, height = 4.2)

# d: New discoveries table
promising_new_rejections %>% 
  select(pair_number, gene_short_name, target_site, corrected_pvalue_st, pvalue.empirical, eQTL, eRNA) %>% 
  mutate_at(c("corrected_pvalue_st", "pvalue.empirical", "eQTL", "eRNA"), funs(format(signif(., 2), scientific = TRUE))) %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE, linesep = "",digits = 2,
        col.names = c("", "Gene", "Enhancer", "SCEPTRE", "Original", "eQTL", "eRNA")) %>%
  save_kable(sprintf("%s/figures/Figure4d.pdf", base_dir))


# e: HI-C interaction frequency comparison
old_scores = rejected_pairs_HIC %>% filter(rejected_old, !is.na(score_rank)) %>% pull(score_rank)
new_scores = rejected_pairs_HIC %>% filter(rejected_new, !is.na(score_rank)) %>% pull(score_rank)
x = seq(0, 1 + 1e-10, length.out = 1000)
old_ecdf = ecdf(old_scores)(x)
new_ecdf = ecdf(new_scores)(x)
p = tibble(x, old_ecdf, new_ecdf) %>% 
  gather(method, ecdf, -x) %>% 
  mutate(method = factor(method, levels = c("new_ecdf", "old_ecdf"), labels = c("SCEPTRE", "Original"))) %>%
  ggplot(aes(x = x, y = ecdf, group = method, colour = method)) + 
  geom_line() + geom_abline(slope = 1, linetype = "dashed") + 
  scale_colour_manual(values = c("blue", "red"), name = "Loci pairs") + 
  xlab("Rank of loci pair by HI-C interaction frequency within TAD") + 
  ylab("Cumulative fraction of loci pairs") +
  theme_bw() + theme(
                     legend.position = c(0.8,0.2),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure4e.pdf", base_dir), plot = p, device = "pdf",
       width = 4.5, height = 3)

# f: HIC table
rejected_pairs_HIC %>% 
  summarise(Original = sum(rejected_old & !is.na(TAD_left)), 
            SCEPTRE = sum(rejected_new & !is.na(TAD_left))) %>%
  gather(method, Hits) %>%
  left_join(rejected_pairs_HIC %>% 
              summarise(Original = sum(rejected_old), 
                        SCEPTRE = sum(rejected_new)) %>%
              gather(method, Total), 
            by = "method") %>%
  mutate(Fraction = Hits/Total) %>%
  kable(format = "latex", booktabs = TRUE, escape = FALSE, 
        linesep = "",digits = 2, col.names = c("", "Same TAD", "Total", "Fraction")) %>%
  add_header_above(c(" ", "Gene-enhancer pairs" = 3)) %>%
  save_kable(sprintf("%s/figures/Figure4f.pdf", base_dir))


# g: ChIP-seq odds ratios
p = TF_enrichments %>%
  mutate(method = factor(method, levels = c("Proposed", "Original"), labels =c("SCEPTRE", "Original"))) %>%
  ggplot(aes(x = TF, y = enrichment, fill = method)) + 
  geom_col(position = "dodge") + 
  xlab("Transcription Factor") + ylab("Enrichment (odds ratio)\namong paired enhancers") +
  scale_fill_manual(values = c("blue", "red"), name = "Paired enhancers") + 
  scale_y_continuous(expand = c(0,0)) + 
  theme_bw() + theme(
                     legend.position = c(0.825,0.8),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)
ggsave(plot = p, filename = sprintf("%s/figures/Figure4g.pdf", base_dir), width = 4.5, height = 3)