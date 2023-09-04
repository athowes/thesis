library(ggplot2)

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

a <- 3
b <- 1

set.seed(1)

truth <- 2.5
y <- rpois(3, lambda = truth)

fn <- function(x) dgamma(x, a + sum(y), b + length(y), log = TRUE)

# Here we are using numerical derivatives
ff <- list(
  fn = fn,
  gr = function(x) numDeriv::grad(fn, x),
  he = function(x) numDeriv::hessian(fn, x)
)

opt_bfgs <- aghq::optimize_theta(
  ff, 1, control = aghq::default_control(method = "BFGS")
)

ggplot(data = data.frame(x = c(0, 8)), aes(x)) +
  geom_function(aes(col = "Prior"), fun = dgamma, n = 500, args = list(shape = a, rate = b)) +
  geom_function(aes(col = "Posterior"), fun = dgamma, n = 500, args = list(shape = a + sum(y), rate = b + length(y))) +
  stat_function(aes(col = "Laplace"), fun = dnorm, n = 500, args = list(mean = opt_bfgs$mode, sd = sqrt(1 / opt_bfgs$hessian))) +
  stat_function(fun = unnormalised_posterior, n = 500) +
  geom_vline(aes(col = "Truth", xintercept = truth)) +
  scale_color_manual(values = c("#CC79A7", "#56B4E9","#009E73", "#E69F00")) +
  geom_point(data = data.frame(x = y, y = 0), aes(x = x, y = y), inherit.aes = FALSE, alpha = 0.7, size = 2) +
  geom_col(data = data.frame(x = c(0, 0, 0, 0), y = c(0, 0, 0, 0), type = as.factor(c("Laplace", "Posterior", "Prior", "Truth"))), aes(x = x, y = y, fill = type)) +
  scale_fill_manual(values = c("#CC79A7", "#56B4E9","#009E73", "#E69F00")) +
  labs(x = "", y = "", col = "", fill = "") +
  guides(col = "none", fill = guide_legend(override.aes = list(alpha = 1, shape = 15))) +
  theme_minimal()

ggsave("figures/naomi-aghq/laplace.png", h = 3, w = 6.25, bg = "white")

opt_bfgs$mode
sqrt(1 / opt_bfgs$hessian)

unnormalised_posterior <- function(phi) {
  phi^8 * exp(-4 * phi)
}

lognormconst <- log(gamma(9)) - 9 * log(4)
lognormconst
laplace_lognormconst <- log(unnormalised_posterior(2)) - dnorm(2, mean = opt_bfgs$mode, sd = sqrt(1 / opt_bfgs$hessian), log = TRUE)
laplace_lognormconst
