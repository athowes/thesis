library(tidyverse)
library(patchwork)
library(INLA)
library(TMB)
library(gt)

data("Epil")

Epil %>%
  ggplot(aes(x =  as.factor(Trt), y = y)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 2, shape = 1) +
  # scale_colour_manual(values = c("grey", "#56B4E9")) +
  coord_flip() +
  labs(x = "Treatment", y = "Number of seizures") +
  theme_minimal()

ggsave("figures/naomi-aghq/epil.png", h = 3.5, w = 6.25, bg = "white")

centre <- function(x) (x - mean(x))

Epil <- Epil %>%
  mutate(CTrt    = centre(Trt),
         ClBase4 = centre(log(Base/4)),
         CV4     = centre(V4),
         ClAge   = centre(log(Age)),
         CBT     = centre(Trt * log(Base/4)))

N <- 59
J <- 4
K <- 6
X <- model.matrix(formula(~ 1 + CTrt + ClBase4 + CV4 + ClAge + CBT), data = Epil)
y <- Epil$y

make_epsilon_matrix <- function(N, J) {
  t(outer(1:N, 1:(N * J), function(r, c) as.numeric((J*(r - 1) < c) & (c <= J*r))))
}

dat <- list(N = N, J = J, K = K, X = X, y = y, E = make_epsilon_matrix(N, J))

tau_prior <- list(prec = list(
  prior = "loggamma",
  param = c(0.001, 0.001),
  initial = 1,
  fixed = FALSE)
)

formula <- y ~ 1 + CTrt + ClBase4 + CV4 + ClAge + CBT +
  f(rand, model = "iid", hyper = tau_prior) +
  f(Ind,  model = "iid", hyper = tau_prior)

beta_prior <- list(mean = 0, prec = 1 / 100^2)

epil_inla <- function(strat, int_strat) {
  inla(
    formula,
    control.fixed = beta_prior,
    family = "poisson",
    data = Epil,
    control.inla = list(strategy = strat, int.strategy = int_strat),
    control.predictor = list(compute = TRUE)
  )
}

fit <- epil_inla(strat = "gaussian", int_strat = "eb")

compile("resources/naomi-aghq/epil.cpp")
dyn.load(dynlib("resources/naomi-aghq/epil"))

param <- list(
  beta = rep(0, K),
  epsilon = rep(0, N),
  nu = rep(0, N * J),
  l_tau_epsilon = 0,
  l_tau_nu = 0
)

obj <- MakeADFun(
  data = dat,
  parameters = param,
  random = c("beta", "epsilon", "nu"),
  DLL = "epil",
  silent = TRUE
)

its <- 1000

opt <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(iter.max = its, trace = 0)
)

sd_out <- sdreport(
  obj,
  par.fixed = opt$par,
  getJointPrecision = TRUE
)

stan <- tmbstan::tmbstan(obj = obj, chains = 4, refresh = 0)

df <- cbind(
  "R-INLA" = as.vector(t(fit$summary.fixed[1:6, 1:2])),
  "TMB" = as.vector(t(data.frame(sd_out$par.random[1:6], sqrt(sd_out$diag.cov.random[1:6])))),
  "NUTS" = as.vector(t(summary(stan)$summary[1:6, c(1, 3)]))
) %>% as.data.frame()

beta_i <- function(i) { c(paste0("beta[", i, "]"), paste0("sd(beta[", i, "])")) }
rownames(df) <- c(sapply(0:5, beta_i))

saveRDS(df, "resources/naomi-aghq/epil.rds")
