#######################################################
#
# Reproduce Figure 2c from Katsevich and Roeder (2020).
#
#######################################################

gene_to_plot = "ENSG00000135390"
grna_to_plot = "chr12.1893_top_two"
p = plot_skew_t(gene_to_plot, grna_to_plot) + theme(plot.title = element_blank())
plot(p)
ggsave(filename = sprintf("%s/figures/Figure2c.png", base_dir), plot = p, device = "png", 
       width = 3.5, height = 2.5)
