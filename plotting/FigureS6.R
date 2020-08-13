#######################################################
#
# Reproduce Figure S6 from Katsevich and Roeder (2020).
#
#######################################################

# Collate candidate enhancer p-values from negative binomial and 
# conditional randomization approaches
df = likelihood_results %>% 
  filter(site_type == "DHS", method == "raw_2") %>% 
  select(gene_id, grna_group, pvalue, zvalue) %>% 
  dplyr::rename(old_pvalue = pvalue, old_zvalue = zvalue) %>%
  left_join(resampling_results %>% 
              filter(site_type == "DHS", method == "conditional_randomization") %>% 
              select(gene_id, grna_group, corrected_pvalue_st, original_zvalue, xi, alpha, omega, nu) %>% 
              dplyr::rename(new_pvalue = corrected_pvalue_st, new_zvalue = original_zvalue), 
            by = c("grna_group", "gene_id"))

# three gene/gRNA pairs to illustrate the differences in the two approaches
genes_to_plot = c("ENSG00000124657","ENSG00000145592","ENSG00000065268")
grnas_to_plot = c("chr6.1215_top_two", "chr5.1161_top_two","chr19.272_top_two")

# panel a: comparing the two ways of computing p-values
p0 = df %>% 
  mutate(highlighted = 
           (gene_id == genes_to_plot[1] & grna_group == grnas_to_plot[1]) |
           (gene_id == genes_to_plot[2] & grna_group == grnas_to_plot[2]) |
           (gene_id == genes_to_plot[3] & grna_group == grnas_to_plot[3])) %>%
  arrange(highlighted) %>%
  filter(old_pvalue > 1e-10, new_pvalue > 1e-10, !is.na(alpha), alpha < 4) %>% 
  ggplot(aes(x = old_pvalue, y = new_pvalue, colour = highlighted)) + 
  geom_point() + geom_abline(slope = 1, linetype = "dashed") +
  scale_x_log10() + scale_y_log10() + 
  xlab("Negative binomial p-value") + ylab("SCEPTRE p-value") +
  scale_colour_manual(values = c("dodgerblue", "red")) + 
  theme_bw() + theme(legend.position = "none",
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p0)

# panel b: NB p-value more significant
p1 = plot_skew_t(genes_to_plot[3], grnas_to_plot[3])  
plot(p1)

# panel c: SCEPTRE p-value more significant
p2 = plot_skew_t(genes_to_plot[2], grnas_to_plot[2]) + theme(legend.position = "none")
plot(p2)

# panel d: two p-values about equal
p3 = plot_skew_t(genes_to_plot[1], grnas_to_plot[1]) + theme(legend.position = "none")  
plot(p3)

# combine the panels
p = grid.arrange(p0, p1, p2, p3, nrow = 2)
plot(p)
ggsave(filename = sprintf("%s/figures/FigureS6.png", base_dir), plot = p, device = "png", 
       width = 6.5, height = 6.5)