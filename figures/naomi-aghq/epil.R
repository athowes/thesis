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

ggsave("figures/naomi-aghq/epil.png", h = 3.5, w = 6.25)

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

start <- Sys.time() 
inla_g_eb <- epil_inla(strat = "gaussian", int_strat = "eb")
end <- Sys.time()
inla_g_eb_time <- end - start

start <- Sys.time() 
inla_sl_eb <- epil_inla(strat = "simplified.laplace", int_strat = "eb") 
end <- Sys.time()
inla_sl_eb_time <- end - start

start <- Sys.time() 
inla_l_eb <- epil_inla(strat = "laplace", int_strat = "eb") 
end <- Sys.time()
inla_l_eb_time <- end - start

start <- Sys.time() 
inla_g_grid <- epil_inla(strat = "gaussian", int_strat = "grid")
end <- Sys.time()
inla_g_grid_time <- end - start

start <- Sys.time() 
inla_sl_grid <- epil_inla(strat = "simplified.laplace", int_strat = "grid") 
end <- Sys.time()
inla_sl_grid_time <- end - start

start <- Sys.time() 
inla_l_grid <- epil_inla(strat = "laplace", int_strat = "grid") 
end <- Sys.time()
inla_l_grid_time <- end - start

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

start <- c(param$l_tau_epsilon, param$l_tau_nu)
aghq <- aghq::marginal_laplace_tmb(obj, k = 3, startingvalue = start)

aghq_samples <- aghq::sample_marginal(aghq, M = 1000)$samps %>%
  t() %>%
  as.data.frame() %>%
  inf.utils::replace_duplicate_colnames()

stan <- tmbstan::tmbstan(obj = obj, chains = 4, refresh = 0)

df <- cbind(
  "R-INLA" = as.vector(t(fit$summary.fixed[1:6, 1:2])),
  "TMB" = as.vector(t(data.frame(sd_out$par.random[1:6], sqrt(sd_out$diag.cov.random[1:6])))),
  "AGHQ" = as.vector(t(data.frame(mean = apply(aghq_samples[, 1:6], 2, mean), sd = apply(aghq_samples[, 1:6], 2, sd)))),
  "NUTS" = as.vector(t(summary(stan)$summary[1:6, c(1, 3)]))
) %>% as.data.frame()

beta_i <- function(i) { c(paste0("beta[", i, "]"), paste0("sd(beta[", i, "])")) }
rownames(df) <- c(sapply(0:5, beta_i))

saveRDS(df, "resources/naomi-aghq/epil.rds")

time_df <- data.frame(
  time = c(inla_g_eb_time, inla_sl_eb_time, inla_l_eb_time, inla_g_grid_time, inla_sl_grid_time, inla_l_grid_time),
  method = c("Gaussian, EB", "Simplified Laplace, EB", "Laplace, EB", "Gaussian, grid", "Simplified Laplace, grid", "Laplace, grid"),
  software = rep("R-INLA", times = 6)
)

time_df %>% ggplot(aes(x = forcats::fct_reorder(method, time), y = time, fill = software)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values = c("#56B4E9","#009E73", "#E69F00")) +
  labs(x = "Method", y = "Time taken (s)", fill = "Software") +
  coord_flip()

ggsave("figures/naomi-aghq/epil-time.png", h = 3.5, w = 6.25)
