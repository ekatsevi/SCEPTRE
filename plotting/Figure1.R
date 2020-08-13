#######################################################
#
# Reproduce Figure 1 from Katsevich and Roeder (2020).
#
#######################################################

### a: DHS/NTC qq-plot ###
ci = 0.95
subsampling_factor = 250
p = original_results %>%
  rename(pvalue = pvalue.raw) %>%
  filter(site_type %in% c("DHS", "NTC"), beta < 0) %>%
  mutate(site_type = factor(site_type, levels = c("DHS", "NTC"), labels = c("Candidate enhancer", "Negative control"))) %>%
  group_by(site_type) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%
  mutate(clower = ifelse(site_type == "Candidate enhancer", NA, clower), 
         cupper = ifelse(site_type == "Candidate enhancer", NA, cupper)) %>%
  mutate(pvalue = ifelse(pvalue < 1e-50, 0, pvalue)) %>%
  filter(-log10(expected) > 2 | (site_type == "Candidate enhancer" & row_number() %% 20 == 0) | row_number() %% subsampling_factor == 0) %>%
  ggplot(aes(x = expected, y = pvalue, group = site_type, ymin = clower, ymax = cupper)) + 
  geom_point(aes(colour = site_type), size = 3, alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("firebrick3", "slategray"), name = "Perturbation target") + 
  xlab("Expected null p-value") + 
  ylab("Observed p-value for gene-gRNA pair") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(trans = revlog_trans(base = 10)) + 
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  theme_bw() + theme(legend.position = c(0.25,0.8), text = element_text(size = 14),
                     legend.background = element_rect(fill = "transparent", colour = NA),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/Figure1a.png", base_dir), plot = p, device = "png",
       width = 4.5, height = 4)

# b: dispersion estimation
KS_pvals = original_results %>% 
  filter(site_type == "NTC") %>% 
  group_by(gene_id, outlier_gene) %>% 
  summarise(KS_pval = ks.test(pvalue.empirical, "punif")$p.value) %>%
  ungroup()
df = disp_table %>% 
  rename(disp_raw = disp) %>% 
  mutate(disp_shrunk = disp_coeffs[1] + disp_coeffs[2]/mu) %>% 
  inner_join(KS_pvals, by = "gene_id") %>%
  mutate(KS_pval_adj = p.adjust(KS_pval, "fdr"))
problematic_genes = df %>% filter(KS_pval_adj < 0.1) %>% pull(gene_id)
p = df %>% arrange(desc(KS_pval)) %>% mutate(KS_pval = ifelse(KS_pval < 1e-5, 1e-5, KS_pval)) %>% 
  mutate(outlier_gene = factor(outlier_gene, levels = c(FALSE, TRUE), 
                               labels = c("Non-\"outlier\" gene", "\"Outlier\" gene"))) %>%
  ggplot() + geom_point(aes(x = mu, y = disp_raw, colour = KS_pval, shape = outlier_gene)) + 
  geom_point(aes(x = 1.29, y = 3.98), colour = "red") + 
  geom_point(aes(x = 0.964, y = 0.304), colour = "red") + 
  geom_point(aes(x = mu, y = disp_raw), size = 4, shape = 21, data = df %>% filter(KS_pval < 0.05/num_genes)) + 
  geom_line(aes(x = mu, y = disp_coeffs[1] + disp_coeffs[2]/mu), colour = "black", linetype = "dashed") + 
  scale_colour_continuous(trans = "log10", name = "Empirical p-value\nmiscalibration, per gene") +
  scale_x_log10() + scale_y_log10() + 
  xlab("Mean gene expression") + ylab("Dispersion of gene expression") + 
  guides(colour = guide_colorbar(title.position = "top", direction = "horizontal", order = 1)) +
  guides(shape = guide_legend(title = NULL, direction = "vertical", order = 2)) +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = c(0.65, 0.175),
                     panel.border = element_blank(), axis.line = element_line(),
                     legend.box = "horizontal",  legend.key = element_rect(colour = "transparent", fill = "white"),
                     text = element_text(size = 14),
                     legend.title = element_text(size = 12),
                     legend.background = element_rect(fill = "transparent", colour = NA))
plot(p)
ggsave(filename = sprintf("%s/figures/Figure1b.png", base_dir), plot = p, device = "png",
       width = 7, height = 5)

# c,d: undercorrection + overcorrection

thresh_old = original_results %>%
  filter(site_type == "DHS", pvalue.empirical.adjusted <= 0.1) %>%
  summarise(max(pvalue.empirical)) %>%
  pull()
undercorrection_gene = "ENSG00000124575"
overcorrection_gene = "ENSG00000146963"
df = original_results %>% 
  filter((gene_id == undercorrection_gene & (site_type == "NTC" | (site_type == "DHS" & beta < 0))) |
           (gene_id == overcorrection_gene & site_type == "NTC")) %>% 
  # mutate(pvalue.empirical = ifelse(site_type == "NTC", null_ecdf(pvalue.raw), pvalue.empirical)) %>%
  gather("pvalue_type", "pvalue", pvalue.raw, pvalue.empirical) %>%
  group_by(gene_id, site_type, pvalue_type) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%
  mutate(clower = ifelse(pvalue_type == "pvalue.empirical", NA, clower), 
         cupper = ifelse(pvalue_type == "pvalue.empirical", NA, cupper)) %>%
  mutate(site_type = factor(site_type, levels = c("NTC", "DHS"), labels = c("NTC gRNAs", "DHS gRNAs"))) %>%
  mutate(pvalue_type = factor(pvalue_type, levels = c("pvalue.raw", "pvalue.empirical"), 
                              labels = c("Raw p-value", "Empirical p-value"))) %>%
  group_by(grna_group) %>%
  mutate(grna_group_name = ifelse(site_type == "DHS gRNAs" & pvalue.empirical.adjusted < 0.1 & pvalue_type == "Raw p-value",
                                  strsplit(unique(grna_group), "_")[[1]][1], ""),
         gene_short_name = factor(ifelse(gene_short_name == "C7orf55-LUC7L2", "LUC7L2", "HIST1H1D"),
                                  levels = c("LUC7L2", "HIST1H1D"))) %>%
  ungroup() %>%
  select(gene_short_name, grna_group, grna_group_name, site_type, pvalue_type, pvalue, expected, clower, cupper)
hline_data = original_results %>% 
  filter(gene_id == overcorrection_gene, site_type == "selfTSS") %>% 
  select(pvalue.raw, pvalue.empirical, gene_short_name) %>% 
  gather(pvalue_type, pvalue, -gene_short_name) %>%
  mutate(site_type = factor("NTC", levels = c("NTC", "DHS"), labels = c("NTC gRNAs", "DHS gRNAs")),
         pvalue_type = factor(pvalue_type, levels = c("pvalue.raw", "pvalue.empirical"),
                              labels = c("Raw p-value", "Empirical p-value")),
         gene_short_name = factor(ifelse(gene_short_name == "C7orf55-LUC7L2", "LUC7L2", "HIST1H1D"),
                                  levels = c("LUC7L2", "HIST1H1D")))
arrow_data = df %>% 
  filter(gene_short_name == "HIST1H1D") %>% 
  select(grna_group, gene_short_name, expected, site_type, pvalue_type, pvalue) %>% 
  spread(pvalue_type, pvalue) %>%
  group_by(site_type) %>%
  arrange(`Raw p-value`) %>% 
  filter(row_number() <= 2) %>% 
  ungroup() %>% 
  mutate(y = `Raw p-value`, yend = `Empirical p-value`, x = expected, xend = expected) %>%
  select(gene_short_name, site_type, x,y,xend,yend) %>%
  add_row(gene_short_name = "LUC7L2", site_type = factor("NTC", levels = c("NTC", "DHS"), labels = c("NTC gRNAs", "DHS gRNAs")),
          x = 0.1, y = hline_data %>% filter(pvalue_type == "Raw p-value") %>% pull(pvalue), 
          xend = 0.1, yend = hline_data %>% filter(pvalue_type == "Empirical p-value") %>% pull(pvalue)) %>%
  mutate(dy = yend - y) 
p = df %>%
  ggplot(aes(x = expected, y = pvalue, ymin = clower, ymax = cupper,
             group = pvalue_type, label = grna_group_name)) + 
  geom_point(aes(colour = pvalue_type), size = 3, alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_text_repel() + 
  geom_text_repel(aes(x = x, y = y, label = label), 
                  data = tibble(x = 0.8, y = 10^(-3.6), 
                                label = "gRNA targeting   LUC7L2 TSS",
                                site_type = factor("NTC", levels = c("NTC", "DHS"), 
                                                   labels = c("NTC gRNAs", "DHS gRNAs")),
                                gene_short_name = factor("LUC7L2", levels = c("LUC7L2", "HIST1H1D"))), 
                  inherit.aes = FALSE) +
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(site_type ~ gene_short_name, scales = "free", nrow = 2, drop = FALSE, labeller = labeller(.multi_line = FALSE)) + 
  geom_hline(aes(yintercept = thresh_old), linetype = "dashed",
             data = tibble(thresh_old = thresh_old,
                           site_type = factor("DHS gRNAs",levels = c("NTC gRNAs", "DHS gRNAs")),
                           gene_short_name = factor("HIST1H1D", levels = c("LUC7L2", "HIST1H1D")))) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(angle = 20, type = "closed", length = unit(100*arrow_data$dy, "in")),
               data = arrow_data, inherit.aes = FALSE) +
  scale_colour_manual(values = c("darkorange1", "dodgerblue")) +
  guides(color = guide_legend(override.aes = list(linetype = 0))) + 
  geom_hline(aes(yintercept = pvalue, colour = pvalue_type), data = hline_data) +
  xlab("Expected null p-value") + 
  ylab("Observed p-value") + 
  scale_x_continuous(trans = revlog_trans(base = 10)) + 
  scale_y_continuous(trans = revlog_trans(base = 10), 
                     breaks = c(1,1e-2,1e-4,1e-6), 
                     labels = c("1", "1e-2", "1e-4", "1e-6")) +
  theme_bw() + theme(axis.title.x = element_text(hjust = 0.925),
                     axis.title.y = element_text(hjust = 0.9),
                     legend.title = element_blank(),
                     legend.position = c(0.7, 0.96),
                     legend.background = element_rect(fill = "transparent", colour = NA),
                     text = element_text(size = 14),
                     strip.background = element_blank(),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure1c-d.pdf", base_dir), plot = p, device = "pdf",
       width = 6, height = 6)


### e: gRNAs versus total UMIs
colnames(df) = c("guide_count", "total_umis")
p = confounders_full %>% 
  as_tibble() %>%
  filter(guide_count > 0, row_number() %% 100 == 0) %>%
  ggplot(aes(x = total_umis, y = guide_count)) + 
  stat_density2d(aes(fill = ..density..^0.25), adjust = 0.5, geom = "raster", contour = FALSE) +
  geom_smooth(method = "lm", colour = "black", se = FALSE, linetype = "dashed") + 
  scale_x_log10(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  xlab("Total UMIs per cell") + ylab("Total gRNAs per cell") + 
  scale_fill_gradient(low = "white", high = "dodgerblue4") + 
  theme_bw() + theme(legend.position = "none", text = element_text(size = 14), 
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure1e.png", base_dir), plot = p, device = "png", 
       width = 4.5, height = 4)

### f: fitted size factors do not cut it

conditional_effects_filename = sprintf("%s/figures/results/gRNA_confounding.tsv", base_dir)
grna_confounding = read_tsv(conditional_effects_filename, col_types = "cd")
grna_groups = c("scrambled_21", "scrambled_6", "random_9")
labels = sapply(grna_groups, function(grna_group)(sprintf("%s (z = %0.2f)", 
                                                          grna_group,
                                                          grna_confounding %>% 
                                                            filter(grna_group == !!grna_group) %>% 
                                                            pull(conditional_effect))))
ci = 0.95
p = likelihood_results %>% 
  filter(method == "raw_1", site_type == "NTC", 
         grna_group %in% grna_groups) %>% 
  mutate(grna_group = factor(grna_group, 
                             levels = grna_groups,
                             labels = labels)
                             # labels = c("scrambled_21 (z = -3.25)",
                             #            "scrambled_6 (z = -0.53)",
                             #            "random_9 (z = 3.30)"))
         ) %>%
  group_by(grna_group) %>%
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%
  mutate(clower = ifelse(grna_group == labels[1], clower, NA),
         cupper = ifelse(grna_group == labels[1], cupper, NA)) %>%
  ggplot(aes(x = expected, y = pvalue, group = grna_group, 
             ymin = clower, ymax = cupper)) + 
  geom_point(aes(colour = grna_group), size = 3, alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("skyblue1", "skyblue4", "slateblue4"), name = "Negative control gRNA") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  xlab("Expected null p-value") + 
  ylab("Observed gene-NTC p-value") + 
  scale_x_continuous(trans = revlog_trans(base = 10)) + 
  scale_y_continuous(trans = revlog_trans(base = 10)) +  
  theme_bw() + theme(legend.position = c(0.3, 0.8), text = element_text(size = 14),
                     legend.background = element_rect(fill = "transparent", colour = NA),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(), 
                     axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure1f.png", base_dir), plot = p, device = "png", 
       width = 4.5, height = 4)