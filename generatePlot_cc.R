pacman::p_load(tidyverse, viridis, hrbrthemes, patchwork)


# popu1
load("cc_pvalues.RData")

plotobj_cc <- ggplot(cc_mat_log, 
                          aes(x = expected, 
                              y = observed)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~ `Test`) + 
  labs() + 
  xlab("Expected -log10 p-values") +
  ylab("Observed -log10 p-values") +
  theme_bw() + 
  theme(legend.position = "none",
        line = element_blank(),
        text = element_text(size = 25)) + 
  scale_colour_manual(values = c("#CF5D6D", "#4896E0"))

ggsave(filename = "cc_plot.png", device = 'png', width = 15, height = 10)
