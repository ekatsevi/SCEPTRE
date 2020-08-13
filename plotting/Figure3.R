#######################################################
#
# Reproduce Figure 3 from Katsevich and Roeder (2020).
#
#######################################################

# collate the NTC p-values of four methods
df_NTC = rbind(resampling_results %>% 
                 filter(site_type == "NTC", method == "conditional_randomization") %>%
                 select(gene_id, grna_group, corrected_pvalue_st) %>% 
                 rename(pvalue = corrected_pvalue_st) %>%
                 mutate(method = "Conditional randomization"),
               resampling_results %>% 
                 filter(site_type == "NTC", method == "marginal_permutation") %>%
                 select(gene_id, grna_group, corrected_pvalue_st) %>% 
                 rename(pvalue = corrected_pvalue_st) %>%
                 mutate(method = "Marginal permutation"),
               likelihood_results %>%
                 filter(site_type == "NTC", method == "raw_2") %>%
                 select(gene_id, grna_group, pvalue) %>%
                 mutate(method = "Improved NB"),
               likelihood_results %>%
                 filter(site_type == "NTC", method == "shrunk_1") %>%
                 select(gene_id, grna_group, pvalue) %>%
                 mutate(method = "Original NB")) %>%
  mutate(method = factor(method, levels = c("Original NB", "Improved NB", 
                                            "Marginal permutation", "Conditional randomization")))

# compute information for overall QQ plot
df1 = df_NTC %>% 
  group_by(method) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup()

# compute information for QQ plot by gRNA
df2 = df_NTC %>%   
  group_by(grna_group, method) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%  
  group_by(r, expected, clower, cupper, method) %>% 
  summarise(pmedian = median(pvalue), plower = quantile(pvalue, 0.025), pupper = quantile(pvalue, 0.975)) %>% 
  ungroup()

# compute information for QQ plot by gene
df3 = df_NTC %>%   
  group_by(gene_id, method) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%  
  group_by(r, expected, clower, cupper, method) %>% 
  summarise(pmedian = median(pvalue), plower = quantile(pvalue, 0.025), pupper = quantile(pvalue, 0.975)) %>% 
  ungroup()

# panel a: simulation results
p0 = simulation_results %>% 
  mutate(experiment = ifelse(estimated_size == 1 & zero_inflation == 0, 
                             "Correct expression model",
                             ifelse(estimated_size == 0.2 & zero_inflation == 0,
                                    "Overestimated dispersion",
                                    ifelse(estimated_size == 5 & zero_inflation == 0,
                                           "Underestimated dispersion", "Zero inflation"))),
         experiment = factor(experiment, 
                             levels = c("Correct expression model", 
                                        "Underestimated dispersion",
                                        "Overestimated dispersion",
                                        "Zero inflation"),
                             labels = c("Correct model", 
                                        "Dispersion too small",
                                        "Dispersion too large",
                                        "Zero inflation"))) %>%
  mutate(pvalue = ifelse(pvalue < 1e-5, 1e-5, pvalue)) %>%
  mutate(method = factor(method, levels = c("NB", "MP", "CR"))) %>%
  group_by(method, experiment) %>% 
  mutate(r = rank(pvalue, ties = "first"), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%  
  mutate(clower = ifelse(method == "CR", clower, NA), 
         cupper = ifelse(method == "CR", cupper, NA)) %>%
  ggplot(aes(x = expected, y = pvalue, ymin = clower, ymax = cupper, group = method)) + 
  geom_point(aes(colour = method), alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("red", "turquoise", "blue")) + 
  facet_wrap(. ~ experiment, nrow = 2) + 
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggtitle("Simulated gene-NTC pair") + 
  theme_bw() + theme(legend.position = "none",
                     panel.spacing.x = unit(1.25, "lines"),
                     plot.title = element_text(hjust = 0.5),
                     strip.background = element_blank(),
                     strip.text = element_text(size = 10),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.title = element_blank(),
                     axis.line = element_line())
plot(p0)

# panel b: overall calibration on real data
p1 = df1 %>% 
  filter(-log10(expected) > 2 | row_number() %% subsampling_factor == 0) %>%
  ungroup() %>%
  mutate(clower = ifelse(method == "Conditional randomization", clower, NA), 
         cupper = ifelse(method == "Conditional randomization", cupper, NA)) %>%
  mutate(method = factor(method, 
                         levels = c("Original NB", "Improved NB", 
                                    "Marginal permutation", "Conditional randomization"),
                         labels = c("Original NB", "Improved NB", 
                                    "Permutation", "SCEPTRE"))) %>%
  arrange(method) %>%
  mutate(pvalue = ifelse(pvalue < 1e-8, 1e-8, pvalue)) %>%
  ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) + 
  geom_point(aes(colour = method), size = 1, alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("red", "violet", "turquoise", "blue")) + 
  scale_x_continuous(trans = revlog_trans(base = 10)) + scale_y_continuous(trans = revlog_trans(base = 10)) +
  xlab(expression(paste("Expected null p-value"))) + 
  ylab(expression(paste("Observed p-value"))) + 
  ggtitle("All real gene-NTC pairs") + 
  theme_bw() + theme(legend.position = c(0.275,0.775),
                     legend.background = element_rect(fill = "transparent", colour = NA),
                     legend.text = element_text(size = 12),
                     legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5, size = 12),
                     panel.grid = element_blank(), 
                     strip.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p1)

# panel c: calibration per NTC on real data
p2 = df2 %>%
  mutate_at(c("expected", "pmedian", "plower", "pupper"), ~ifelse(.<1e-6, 1e-6,.)) %>%
  ggplot(aes(x = expected, y = pmedian)) + 
  geom_line(aes(colour = method), size = 1) + 
  geom_ribbon(aes(ymin = plower, ymax = pupper, fill = method), alpha = 0.4) +
  geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("red", "violet", "turquoise", "blue")) + 
  scale_fill_manual(values = c("red", "violet", "turquoise", "blue")) +
  facet_wrap(method~., scales = "fixed", nrow = 2) + 
  scale_x_continuous(trans = revlog_trans(base = 10), breaks = c(1,1e-2,1e-4)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggtitle("Real gene-NTC pairs for each NTC") + 
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(hjust = 0.5, size = 12),
                     panel.spacing.x = unit(1.25, "lines"),
                     axis.title = element_blank(),
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p2)

# panel d: calibration per gene on real data
p3 = df3 %>%
  mutate_at(c("expected", "pmedian", "plower", "pupper"), ~ifelse(.<5e-4, 5e-4,.)) %>%
  ggplot(aes(x = expected, y = pmedian, ymin = plower, ymax = pupper)) + 
  geom_line(aes(colour = method), size = 1) + 
  geom_ribbon(aes(ymin = plower, ymax = pupper, fill = method), alpha = 0.4) +
  geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("red", "violet", "turquoise", "blue")) + 
  scale_fill_manual(values = c("red", "violet", "turquoise", "blue")) +
  facet_wrap(method~., scales = "fixed", nrow = 2) + 
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggtitle("Real gene-NTC pairs for each gene") + 
  theme_bw() + theme(legend.position = "none",
                     panel.spacing.x = unit(1.25, "lines"),
                     plot.title = element_text(hjust = 0.5, size = 12),
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.title = element_blank(),
                     axis.line = element_line())
plot(p3)

# combine the panels
g <- ggplotGrob(p1)$grobs
x_axis_title <- g[[which(sapply(g, function(x) strsplit(x$name, split = "[.][.]")[[1]][1]) == "axis.title.x.bottom")]]
x_axis_title_height <- sum(x_axis_title$height)

y_axis_title <- g[[which(sapply(g, function(x) strsplit(x$name, split = "[.][.]")[[1]][1]) == "axis.title.y.left")]]
y_axis_title_width <- sum(y_axis_title$width)

main_plot = grid.arrange(
  arrangeGrob(p0, p1+theme(axis.title = element_blank()),p2,p3, nrow=2),
  x_axis_title,
  ncol = 1,
  heights = unit.c(unit(1, "npc") - x_axis_title_height, x_axis_title_height))

p = grid.arrange(
  y_axis_title, main_plot, ncol = 2,
  widths = unit.c(y_axis_title_width, unit(1, "npc") - y_axis_title_width))

ggsave(filename = sprintf("%s/figures/Figure3.png", base_dir), plot = p, device = "png",
       width = 7, height = 7)