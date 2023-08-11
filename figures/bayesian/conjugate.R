cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

a <- 3
b <- 1

set.seed(1)

truth <- 2.5
y <- rpois(3, lambda = truth)

prior <- ggplot(data = data.frame(x = c(0, 8)), aes(x)) +
  stat_function(fun = dgamma, n = 500, args = list(shape = a, rate = b), col = "#56B4E9") +
  annotate("text", x = 5.5, y = 0.2, label = "Prior: Gamma(3, 1)", col = "#56B4E9", size = 4) +
  labs(x = "", y = "") +
  theme_minimal()

3 + sum(y)


prior +
  geom_point(data = data.frame(x = y, y = 0), aes(x = x, y = y), inherit.aes = FALSE, alpha = 0.7, size = 2) +
  stat_function(data = data.frame(x = c(0, 8)), aes(x), fun = dgamma, n = 500, args = list(shape = a + sum(y), rate = b + length(y)), col = "#009E73") +
  annotate("text", x = 5, y = 0.3, label = "Posterior: Gamma(9, 4)", col = "#009E73", size = 4) +
  geom_vline(xintercept = truth, col = "#E69F00", linetype = "dashed") +
  annotate("text", x = 3.5, y = 0.4, label = "Truth", col = "#E69F00", size = 4)

ggsave("figures/bayesian/conjugate.png", h = 3, w = 6.25, bg = "white")
