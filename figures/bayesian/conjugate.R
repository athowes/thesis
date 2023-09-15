library(ggplot2)

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

a <- 3
b <- 1

set.seed(1)

truth <- 2.5
y <- rpois(3, lambda = truth)

ggplot(data = data.frame(x = c(0, 8)), aes(x)) +
  geom_function(aes(col = "Prior"), fun = dgamma, n = 500, args = list(shape = a, rate = b)) +
  geom_function(aes(col = "Posterior"), fun = dgamma, n = 500, args = list(shape = a + sum(y), rate = b + length(y))) +
  geom_vline(aes(col = "Truth", xintercept = truth)) +
  scale_color_manual(values = cbpalette) +
  geom_point(data = data.frame(x = y, y = 0), aes(x = x, y = y), inherit.aes = FALSE, alpha = 0.7, size = 2, shape = 1) +
  geom_col(data = data.frame(x = c(0, 0, 0), y = c(0, 0, 0), type = c("Posterior", "Prior", "Truth")), aes(x = x, y = y, fill = type)) +
  scale_fill_manual(values = cbpalette) +
  labs(x = "", y = "", col = "", fill = "") +
  guides(col = "none", fill = guide_legend(override.aes = list(alpha = 1, shape = 15))) +
  theme_minimal()

ggsave("figures/bayesian/conjugate.png", h = 3, w = 6.25, bg = "white")

3 + sum(y)