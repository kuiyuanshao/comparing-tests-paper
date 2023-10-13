pacman::p_load(tidyverse, viridis, hrbrthemes, patchwork)


# popu1
load("popu1_df25_pvalues.RData")
load("popu1_df35_pvalues.RData")
###########################################################


plotobj_25_pois <- ggplot(subset(df25_mat_log, 
                                 subset = `Outcome Distribution` == "Poisson"), 
                          aes(x = expected, 
                              y = observed,
                              colour = `Reference Distribution`)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(`Test` ~ q) + 
  labs() + 
  xlab("Expected -log10 p-values") +
  ylab("Observed -log10 p-values") +
  theme_bw() + 
  theme(legend.position = "none",
        line = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = 35)) + 
  scale_colour_manual(values = c("#CF5D6D", "#4896E0"))


plotobj_35_pois <- ggplot(subset(df35_mat_log, 
                                 subset = `Outcome Distribution` == "Poisson"), 
                          aes(x = expected, 
                              y = observed,
                              colour = `Reference Distribution`)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(`Test` ~ q) + 
  labs() + 
  xlab("Expected -log10 p-values") + 
  theme_bw() + 
  theme(legend.position = "none",
        line = element_blank(),
        axis.title.y = element_blank(), 
        strip.text.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 35)) + 
  scale_colour_manual(values = c("#CF5D6D", "#4896E0"))

plotobj_25_bino <- ggplot(subset(df25_mat_log, 
                                 subset = `Outcome Distribution` == "Binomial"), 
                          aes(x = expected, 
                              y = observed,
                              colour = `Reference Distribution`)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(`Test` ~ q) + 
  labs() + 
  xlab("Expected -log10 p-values") + 
  theme_bw() + 
  theme(legend.position = "none",
        line = element_blank(),
        axis.title.y = element_blank(), 
        strip.text.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 35)) + 
  scale_colour_manual(values = c("#CF5D6D", "#4896E0"))


plotobj_35_bino <- ggplot(subset(df35_mat_log, 
                                 subset = `Outcome Distribution` == "Binomial"), 
                          aes(x = expected, 
                              y = observed,
                              colour = `Reference Distribution`)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(`Test` ~ q) + 
  labs() + 
  xlab("Expected -log10 p-values") + 
  theme_bw() + 
  theme(legend.position = "none",
        line = element_blank(),
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 35)) + 
  scale_colour_manual(values = c("#CF5D6D", "#4896E0"))


finalplot <- plotobj_25_pois + plotobj_35_pois + plotobj_25_bino + plotobj_35_bino + 
  plot_layout(ncol = 4)

ggsave(filename = "popu1_plot.png", device = 'png', width = 30.5, height = 22)

##################################


# popu2
load("popu2_df25_pvalues.RData")
load("popu2_df35_pvalues.RData")


plotobj_25_pois <- ggplot(subset(df25_mat_log, 
                                 subset = `Outcome Distribution` == "Poisson"), 
                          aes(x = expected, 
                              y = observed,
                              colour = `Reference Distribution`)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(`Test` ~ q) + 
  labs() + 
  xlab("Expected -log10 p-values") +
  ylab("Observed -log10 p-values") +
  theme_bw() + 
  theme(legend.position = "none",
        line = element_blank(),
        strip.text.y = element_blank(),
        text = element_text(size = 35)) + 
  scale_colour_manual(values = c("#CF5D6D", "#4896E0"))


plotobj_35_pois <- ggplot(subset(df35_mat_log, 
                                 subset = `Outcome Distribution` == "Poisson"), 
                          aes(x = expected, 
                              y = observed,
                              colour = `Reference Distribution`)) + 
  geom_point(alpha = 0.3) + 
  geom_abline(slope = 1, intercept = 0) + 
  facet_grid(`Test` ~ q) + 
  labs() + 
  xlab("Expected -log10 p-values") + 
  theme_bw() + 
  theme(legend.position = "none",
        line = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        text = element_text(size = 35)) + 
  scale_colour_manual(values = c("#CF5D6D", "#4896E0"))


finalplot <- plotobj_25_pois + plotobj_35_pois+
  plot_layout(ncol = 2)

ggsave(filename = "popu2_plot.png", device = 'png', width = 15, height = 21)
