#######################################################
#
# Reproduce Figure S5 from Katsevich and Roeder (2020).
#
#######################################################

# collate NTC p-values from four negative binomial variants
df_NTC = likelihood_results %>% 
  filter(site_type == "NTC") %>% 
  mutate(method = ifelse(method == "raw_1", 
                         "NB (disp)",
                         ifelse(method == "raw_2",
                                "NB (both)",
                                ifelse(method == "shrunk_1",
                                       "NB (orig)", "NB (conf)")))) %>%
  select(gene_id, grna_group, pvalue, method)

# compute information for overall QQ plot
ci = 0.95
df1 = df_NTC %>% 
  group_by(method) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>%
  ungroup()

# compute information for QQ plot by gene
df2 = df_NTC %>%   
  group_by(gene_id, method) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%  
  group_by(r, expected, clower, cupper, method) %>% 
  summarise(pmedian = median(pvalue), plower = quantile(pvalue, 0.025), pupper = quantile(pvalue, 0.975)) %>% 
  ungroup()

# compute information for QQ plot by gRNA
df3 = df_NTC %>%   
  group_by(grna_group, method) %>% 
  mutate(r = rank(pvalue), expected = ppoints(n())[r], 
         clower = qbeta(p=(1-ci)/2, shape1 = r, shape2 = n()+1-r),
         cupper = qbeta(p=(1+ci)/2, shape1 = r, shape2 = n()+1-r)) %>% 
  ungroup() %>%  
  group_by(r, expected, clower, cupper, method) %>% 
  summarise(pmedian = median(pvalue), plower = quantile(pvalue, 0.025), pupper = quantile(pvalue, 0.975)) %>% 
  ungroup()

# panel a: overall calibration on real data
subsampling_factor = 250
p1 = df1 %>% 
  mutate(method = factor(method, levels = c("NB (orig)", "NB (conf)", "NB (disp)", "NB (both)"),
                         labels = c("Original", "Confounder improved", "Dispersion improved", "Both improved"))) %>%
  filter(-log10(expected) > 2 | row_number() %% subsampling_factor == 0) %>%
  ungroup() %>%
  mutate(clower = ifelse(method == "Original", clower, NA), 
         cupper = ifelse(method == "Original", cupper, NA)) %>%
  arrange(method) %>%
  mutate(pvalue = ifelse(pvalue < 1e-8, 1e-8, pvalue)) %>%
  ggplot(aes(x = expected, y = pvalue, group = method, ymin = clower, ymax = cupper)) + 
  geom_point(aes(colour = method), size = 0.5, alpha = 0.5) + 
  geom_ribbon(alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(trans = revlog_trans(base = 10)) + scale_y_continuous(trans = revlog_trans(base = 10)) +
  xlab(expression(paste("Expected null p-value"))) + 
  ylab(expression(paste("Observed p-value"))) + 
  ggtitle("All gene-NTC pairs") + 
  theme_bw() + theme(legend.position = c(0.32,0.75),
                     legend.background = element_rect(fill = "transparent", colour = NA),
                     legend.title = element_blank(),
                     plot.title = element_text(hjust = 0.5, size = 12),
                     panel.grid = element_blank(), 
                     strip.background = element_blank(),
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p1)

# panel b: calibration by gene
p2 = df2 %>%
  mutate(method = factor(method, levels = c("NB (orig)", "NB (conf)", "NB (disp)", "NB (both)"),
                         labels = c("Original", "Confounder improved", "Dispersion improved", "Both improved"))) %>%
  mutate_at(c("expected", "pmedian", "plower", "pupper"), ~ifelse(.<5e-4, 5e-4,.)) %>%
  ggplot(aes(x = expected, y = pmedian, ymin = plower, ymax = pupper)) + 
  geom_line(aes(colour = method), size = 1) + 
  geom_ribbon(aes(ymin = plower, ymax = pupper, fill = method), alpha = 0.4) +
  geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(method~., scales = "fixed", nrow = 2) + 
  scale_x_continuous(trans = revlog_trans(base = 10)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggtitle("Gene-NTC pairs for each gene") + 
  theme_bw() + theme(legend.position = "none",
                     panel.spacing.x = unit(1.25, "lines"),
                     plot.title = element_text(hjust = 0.5, size = 12),
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.title = element_blank(),
                     axis.line = element_line())
plot(p2)

# panel c: calibration by NTC
p3 = df3 %>%
  mutate(method = factor(method, levels = c("NB (orig)", "NB (conf)", "NB (disp)", "NB (both)"),
                         labels = c("Original", "Confounder improved", "Dispersion improved", "Both improved"))) %>%
  mutate_at(c("expected", "pmedian", "plower", "pupper"), ~ifelse(.<1e-6, 1e-6,.)) %>%
  ggplot(aes(x = expected, y = pmedian)) + 
  geom_line(aes(colour = method), size = 1) + 
  geom_ribbon(aes(ymin = plower, ymax = pupper, fill = method), alpha = 0.4) +
  geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2) + 
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(method~., scales = "fixed", nrow = 1) + 
  scale_x_continuous(trans = revlog_trans(base = 10), breaks = c(1,1e-2,1e-4)) +
  scale_y_continuous(trans = revlog_trans(base = 10)) +
  ggtitle("Gene-NTC pairs for each NTC") + 
  theme_bw() + theme(legend.position = "none",
                     plot.title = element_text(hjust = 0.5, size = 12),
                     panel.spacing.x = unit(1.25, "lines"),
                     axis.title = element_blank(),
                     strip.background = element_blank(),
                     strip.text = element_blank(),
                     panel.grid = element_blank(), 
                     panel.border = element_blank(),
                     axis.line = element_line())
plot(p3)

# combine the three panels
g <- ggplotGrob(p1)$grobs
x_axis_title <- g[[which(sapply(g, function(x) strsplit(x$name, split = "[.][.]")[[1]][1]) == "axis.title.x.bottom")]]
x_axis_title_height <- sum(x_axis_title$height)

y_axis_title <- g[[which(sapply(g, function(x) strsplit(x$name, split = "[.][.]")[[1]][1]) == "axis.title.y.left")]]
y_axis_title_width <- sum(y_axis_title$width)

main_plot = grid.arrange(
  arrangeGrob(arrangeGrob(p1+theme(axis.title = element_blank()),p2,nrow=1), p3, nrow = 2, heights = c(1,0.6)),
  x_axis_title,
  ncol = 1,
  heights = unit.c(unit(1, "npc") - x_axis_title_height, x_axis_title_height))

p = grid.arrange(
  y_axis_title, main_plot, ncol = 2,
  widths = unit.c(y_axis_title_width, unit(1, "npc") - y_axis_title_width))

plot(p)

ggsave(filename = sprintf("%s/figures/FigureS5.png", base_dir), plot = p, device = "png",
       width = 6.5, height = 5)