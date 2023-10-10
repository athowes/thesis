library(tidyverse)
library(patchwork)

library(INLA)
data("Epil")

Epil %>%
  ggplot(aes(x =  as.factor(Trt), y = y)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(col = as.factor(V4)), width = 0.2, alpha = 0.7, size = 2, shape = 1) +
  scale_colour_manual(values = c("grey", "#56B4E9")) +
  coord_flip() +
  labs(x = "Treatment", y = "Number of seizures", col = "Final visit?") +
  theme_minimal()

ggsave("figures/naomi-aghq/epil.png", h = 3.5, w = 6.25, bg = "white")
