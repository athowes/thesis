library(ggplot2)
library(rstan)
library(patchwork)

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

a <- 3
b <- 1

set.seed(1)

truth <- 2.5
y <- rpois(3, lambda = truth)

model <- rstan::stan_model("figures/bayesian/conjugate.stan")
fit <- rstan::sampling(model, data = list(N = 3, y = y, a = a, b = b), chains = 4, iter = 500, warmup = 250)
samples <- as.data.frame(rstan::extract(fit)) %>%
  mutate(
    iter = row_number(),
    mean = cumsum(phi) / iter
  )

posterior <- ggplot(data = data.frame(x = c(0, 8)), aes(x)) +
  geom_histogram(data = samples, aes(x = phi, y = stat(density)), bins = 20, fill = "#CC79A7", alpha = 0.9) +
  geom_function(aes(col = "Posterior"), fun = dgamma, n = 500, args = list(shape = a + sum(y), rate = b + length(y))) +
  scale_color_manual(values = c("#56B4E9", "#CC79A7")) +
  geom_col(data = data.frame(x = c(0, 0), y = c(0, 0), type = c("Exact posterior", "MCMC posterior")), aes(x = x, y = y, fill = type)) +
  scale_fill_manual(values = c("#56B4E9", "#CC79A7")) +
  labs(x = "", y = "", col = "", fill = "", tag = "A") +
  guides(col = "none", fill = guide_legend(override.aes = list(alpha = 1, shape = 15))) +
  theme_minimal() +
  theme(legend.position = c(0.75, 0.9))

trace <- ggplot(samples, aes(x = iter, y = phi)) +
  geom_line(col = "#CC79A7", alpha = 0.9) +
  theme_minimal() +
  labs(x = "", y = "", tag = "B")

mean <- ggplot(samples, aes(x = iter, y = mean)) +
  geom_line(col = "#CC79A7") +
  geom_segment(aes(x = 0, xend = 1000, y = 9 / 4, yend = 9 / 4), col = "#56B4E9") +
  labs(x = "Iteration", y = "", tag = "C") +
  theme_minimal()

posterior | (trace / mean)

ggsave("figures/bayesian/stan.png", h = 3.5, w = 6.25, bg = "white")
