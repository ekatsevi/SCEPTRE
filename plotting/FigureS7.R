#######################################################
#
# Reproduce Figure S7 from Katsevich and Roeder (2020).
#
#######################################################

# compute KS p-values
KS_pvals = original_results %>% 
  filter(site_type == "NTC") %>% 
  group_by(gene_id) %>% 
  summarise(KS_pval = ks.test(pvalue.empirical, "punif")$p.value) %>% ungroup() 

# collect several pieces of information about each gene-gRNA pair
thresh_old = original_results %>% filter(rejected) %>% summarise(max(pvalue.empirical)) %>% pull()
df = original_results %>% 
  rename(old_rejected = rejected, old_pvalue = pvalue.empirical) %>%
  left_join(grna_confounding, by = "grna_group") %>% 
  left_join(KS_pvals, by = "gene_id") %>% 
  filter(site_type == "DHS", quality_rank_grna == "top_two") %>%
  left_join(resampling_results %>% 
              filter(method == "conditional_randomization") %>%
              select(gene_id, grna_group, rejected) %>% 
              rename(new_rejected = rejected), 
            by = c("gene_id", "grna_group")) %>%
  left_join(rejected_pairs_HIC %>% 
              mutate(different_TADs = is.na(TAD_left)) %>% 
              select(gene_id, target_site, different_TADs), 
            by = c("gene_id", "target_site")) %>% 
  arrange((old_rejected + new_rejected) %% 2, old_rejected + new_rejected) %>%
  mutate(reason_not_rejected = ifelse(outlier_gene, 
                                      "Outlier gene", 
                                      ifelse(beta > 0, 
                                             "Positive effect", 
                                             ifelse(!old_rejected & old_pvalue <= thresh_old, "Low confidence", "none"))),
         reason_not_rejected = factor(reason_not_rejected, levels = c("none", "Positive effect", "Outlier gene", "Low confidence")),
         rejected_by = ifelse(old_rejected, 
                              ifelse(new_rejected, "Both methods", "Original only"), 
                              ifelse(new_rejected, "SCEPTRE only", "Neither method")),
         rejected_by = factor(rejected_by, levels = c("Both methods", "SCEPTRE only",
                                                      "Original only", "Neither method")))

# make the plot
p = df %>%  
  mutate(text_label = ifelse(gene_short_name == "HIST1H1D" & target_site == "chr6.1077", "HIST1H1D/\nchr6.1077","")) %>%
  ggplot(aes(x = conditional_effect, y = KS_pval, 
             colour = rejected_by, shape = reason_not_rejected, 
             size = reason_not_rejected)) + 
  geom_point() + 
  geom_point(aes(x = conditional_effect, y = KS_pval), size = 4, shape = 21, 
             data = df %>% 
               filter(old_rejected | new_rejected, 
                      different_TADs, 
                      conditional_effect < -2 | KS_pval < 0.05), show.legend = FALSE) + 
  geom_text_repel(aes(label = text_label),colour = "black", force = 1, box.padding = 0.4, size = 3, min.segment.length = 0) + 
  geom_vline(xintercept = -2, linetype = "dashed", colour = "black", alpha = 0.5) + 
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "black", alpha = 0.5) +
  scale_shape_manual(values = c("circle", "square", "triangle", "diamond"), 
                     breaks = c("Outlier gene", "Positive effect", "Low confidence"),
                     name = "Reason for original\nnon-discovery") + 
  scale_size_manual(values = c(1.5,1.5,1.5,2.5), 
                    breaks = c("Outlier gene", "Positive effect", "Low confidence"),
                    name = "Reason for original\nnon-discovery") + 
  scale_colour_manual(values = c("purple", "blue", "red", "gray40"), name = "Discovered by") + 
  guides(colour = guide_legend(order = 1)) +
  guides(shape = guide_legend(order = 2)) +
  guides(size = guide_legend(order = 2)) +
  scale_y_log10() + 
  xlab("Residual confounding (z)") + ylab("Empirical p-value miscalibration (KS p-value)") +
  theme_bw() + 
  theme(
    legend.position = c(0.15, 0.35),
    legend.key = element_rect(colour = "transparent", fill = "transparent"),
    panel.grid = element_blank(), 
    panel.border = element_blank(), 
    axis.line = element_line())
plot(p)
ggsave(filename = sprintf("%s/figures/FigureS7.png", base_dir), plot = p, device = "png", 
       width = 6, height = 5)
