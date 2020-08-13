#######################################################
#
# Reproduce Figure S8 from Katsevich and Roeder (2020).
#
#######################################################

# plot the fraction paired in each quintile
p = paired_fractions %>% 
  gather(method, paired_fraction, rejected_old, rejected_new) %>% 
  mutate(method = factor(method, levels = c("rejected_new", "rejected_old"), labels = c("SCEPTRE", "Original"))) %>%
  ggplot(aes(x = factor(quintile), y = paired_fraction, fill = method)) + 
  xlab("ChIP-seq quintiles of candidate enhancers") + ylab("Proportion enhancers paired with gene") +
  geom_col(position = "dodge") + scale_fill_manual(values = c("blue", "red")) + 
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(TF ~ ., nrow = 2) + theme_bw() + 
  theme(legend.title = element_blank(), 
        legend.position = c(0.09, 0.9), 
        strip.background = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line()
  )
plot(p)
ggsave(plot = p, filename = sprintf("%s/figures/FigureS8.pdf", base_dir), width = 7, height = 4)