library(tidyverse)
library(patchwork)
library(INLA)
library(TMB)
library(gt)

set.seed(1)

source("figures/naomi-aghq/functions.R")

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
  mutate(
    CTrt    = centre(Trt),
    ClBase4 = centre(log(Base/4)),
    CV4     = centre(V4),
    ClAge   = centre(log(Age)),
    CBT     = centre(Trt * log(Base/4))
  )

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
    control.predictor = list(compute = TRUE),
    control.compute = list(config = TRUE)
  )
}

tictoc::tic()
inla_g_eb <- epil_inla(strat = "gaussian", int_strat = "eb")
inla_g_eb_time <- tictoc::toc()

tictoc::tic()
inla_sl_eb <- epil_inla(strat = "simplified.laplace", int_strat = "eb") 
inla_sl_eb_time <- tictoc::toc()

tictoc::tic()
inla_l_eb <- epil_inla(strat = "laplace", int_strat = "eb") 
inla_l_eb_time <- tictoc::toc()

tictoc::tic()
inla_g_grid <- epil_inla(strat = "gaussian", int_strat = "grid")
inla_g_grid_time <- tictoc::toc()

tictoc::tic()
inla_sl_grid <- epil_inla(strat = "simplified.laplace", int_strat = "grid") 
inla_sl_grid_time <- tictoc::toc()

tictoc::tic()
inla_l_grid <- epil_inla(strat = "laplace", int_strat = "grid") 
inla_l_grid_time <- tictoc::toc()

compile("figures/naomi-aghq/TMB/epil.cpp")
dyn.load(dynlib("figures/naomi-aghq/TMB/epil"))

tictoc::tic()
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

tmb_g_eb_time <- tictoc::toc()

sd_summary <- summary(sd_out)
tab <- table(rownames(sd_summary))

tmb_parameter_names <- sapply(split(tab, names(tab)), function(x) {
  if(x > 1) paste0(names(x), "[", 1:x, "]")
  else(names(x))
}) %>%
  unlist() %>%
  as.vector()

H <- as.matrix(sd_out$jointPrecision)
rownames(H) <- tmb_parameter_names[1:nrow(H)]
colnames(H) <- tmb_parameter_names[1:ncol(H)]

indices_small <- which(startsWith(tmb_parameter_names, "beta") | startsWith(tmb_parameter_names, "l_"))
H_small <- H[indices_small, indices_small]

reshape2::melt(H_small) %>%
  ggplot(aes(x = Var1, y = Var2, fill = log(value))) +
    geom_tile() +
    coord_fixed(ratio = 1) +
    scale_fill_viridis_c(na.value = "grey90") +
    theme_void() +
    labs(x = "", y = "", fill = "") +
    theme(
      axis.text.y = element_text(color = "grey30", hjust = 1, vjust = 0.5, angle = 0, lineheight = 0.9, size = 8.8),
    )

ggsave("figures/naomi-aghq/hessian-matrix.png", h = 1.5, w = 6.25)

tictoc::tic()
init <- c(param$l_tau_epsilon, param$l_tau_nu)
aghq <- aghq::marginal_laplace_tmb(obj, k = 3, startingvalue = init)
tmb_g_aghq_time <- tictoc::toc()

aghq_samples <- aghq::sample_marginal(aghq, M = 4000)$samps %>%
  t() %>%
  as.data.frame() %>%
  inf.utils::replace_duplicate_colnames()

# Laplace marginals

random <- obj$env$random
N <- length(random)
x_names <- names(obj$env$par[random])
x_lengths <- lengths(param[unique(x_names)])
x_starts <- cumsum(x_lengths) - x_lengths

dat$x_starts <- as.numeric(x_starts)
dat$x_lengths <- as.numeric(x_lengths)

param[unique(x_names)] <- NULL
param$x_minus_i <- rep(0, sum(x_lengths) - 1)
param$x_i <- 0

compile("figures/naomi-aghq/TMB/epil_modified.cpp")
dyn.load(dynlib("figures/naomi-aghq/TMB/epil_modified"))

compute_laplace_marginal <- function(i, quad) {
  dat$i <- i
  
  obj_i <- TMB::MakeADFun(
    data = dat,
    parameters = param,
    random = "x_minus_i",
    DLL = "epil_modified",
    silent = TRUE,
  )
  
  random_i <- obj_i$env$random
  mode_i <- quad$modesandhessians[["mode"]][[1]][-i]
  gg <- create_approx_grid(quad$modesandhessians, i = i, k = 5)
  out <- data.frame(index = i, par = x_names[i], x = mvQuad::getNodes(gg), w = mvQuad::getWeights(gg))
  
  theta_names <- make.unique(names(obj$par), sep = "")
  
  .g <- function(x) {
    lp <- vector(mode = "numeric", length = nrow(quad$modesandhessians))
    
    for(z in 1:nrow(quad$modesandhessians)) {
      theta <- as.numeric(quad$modesandhessians[z, theta_names])
      obj_i$env$last.par[random_i] <- quad$modesandhessians[z, "mode"][[1]][-dat$i]
      lp[z] <- as.numeric(- obj_i$fn(c(x, theta)))
    }
    
    return(logSumExpWeights(lp, w = quad$normalized_posterior$nodesandweights$weights))
  }
  
  out$lp <- purrr::map_dbl(out$x, .g)
  lognormconst <- logSumExpWeights(out$lp, out$w)
  out$lp_normalised <- out$lp - lognormconst
  
  return(out)
}

#' Laplace marginals and EB with TMB

tictoc::tic()
init <- c(param$l_tau_epsilon, param$l_tau_nu)
eb <- aghq::marginal_laplace_tmb(obj, k = 1, startingvalue = init)
eb_laplace_marginals <- purrr::map(.x = 1:N, .f = compute_laplace_marginal, quad = eb, .progress = TRUE)
eb_laplace_marginals <- dplyr::bind_rows(eb_laplace_marginals)
tmb_l_eb_time <- tictoc::toc()

#' Laplace marginals and AGHQ with TMB

tictoc::tic()
aghq <- aghq::marginal_laplace_tmb(obj, k = 3, startingvalue = init)
quad_laplace_marginals <- purrr::map(.x = 1:N, .f = compute_laplace_marginal, quad = aghq, .progress = TRUE)
tmb_l_aghq_time <- tictoc::toc()

#' NUTS

tictoc::tic()
tmbstan <- tmbstan::tmbstan(obj = obj, chains = 4, refresh = 0)
tmbstan_time <- tictoc::toc()

tictoc::tic()
stan <- rstan::stan(file = "figures/naomi-aghq/epil.stan", data = dat, chains = 4, refresh = 0)
stan_time <- tictoc::toc()

#' Means and SD
sample_adam <- function(i, M) {
  q <- runif(M)
  pdf_and_cdf <- compute_pdf_and_cdf(nodes = quad_laplace_marginals[[i]]$x, quad_laplace_marginals[[i]]$lp_normalised)
  s <- numeric(length(q))
  for(j in 1:length(q)) s[j] <- pdf_and_cdf$x[max(which(pdf_and_cdf$cdf < q[j]))]
  return(s)
}

samples_adam <- lapply(1:6, sample_adam, M = 4000)

beta_df <- function(value, method, software) {
  data.frame(
    index = rep(0:5, each = 2),
    type = rep(c("Posterior mean", "Posterior SD"), 6),
    value = value,
    method = method,
    software = software
  )
}

df <- bind_rows(
  beta_df(as.vector(t(data.frame(mean = apply(aghq_samples[, 1:6], 2, mean), sd = apply(aghq_samples[, 1:6], 2, sd)))), "Gaussian", "TMB"),
  beta_df(as.vector(sapply(samples_adam, function(x) c(mean(x), sd(x)))), "Laplace", "TMB"),
  beta_df(as.vector(t(inla_g_grid$summary.fixed[1:6, 1:2])), "Gaussian", "R-INLA"),
  beta_df(as.vector(t(inla_l_grid$summary.fixed[1:6, 1:2])), "Laplace", "R-INLA")
) %>%
  left_join(
    data.frame(
      index = rep(0:5, each = 2),
      type = rep(c("Posterior mean", "Posterior SD"), 6),
      gold = c(t(summary(tmbstan)$summary[1:6, c(1, 3)]))
    ),
    by = c("index", "type")
  )
  
ggplot(df, aes(x = as.factor(index), y = (gold - value) / gold, color = method)) +
  geom_point(size = 1.5, shape = 1) +
  facet_grid(type ~ software) +
  scale_color_manual(values = c("#56B4E9","#009E73")) +
  geom_abline(intercept = 0, slope = 0, col = "#E69F00", linetype = "dashed") +
  labs(x = "Regression parameter", y = "Percentage difference to NUTS", col = "Method") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.spacing = unit(1.5, "lines"))

ggsave("figures/naomi-aghq/beta-mean-sd.png", h = 3.5, w = 6.25, bg = "white")

intercept_samples_inla <- function(fit) {
  samples <- INLA::inla.posterior.sample(n = 1000, result = fit)
  intercepts <- INLA::inla.posterior.sample.eval("Intercept", samples)
  as.numeric(intercepts)
}

grid <- seq(from = round(min(aghq_samples[[1]]) - 0.1, digits = 2), to = round(max(aghq_samples[[1]]) + 0.1, digits = 2), length.out = 1000)

inla_g_grid_ecdf <- stats::ecdf(intercept_samples_inla(inla_g_grid))
inla_g_grid_ecdf_df <- data.frame(x = grid, ecdf = inla_g_grid_ecdf(grid), method = "Gaussian", software = "R-INLA")

inla_l_grid_ecdf <- stats::ecdf(intercept_samples_inla(inla_l_grid))
inla_l_grid_ecdf_df <- data.frame(x = grid, ecdf = inla_l_grid_ecdf(grid), method = "Laplace", software = "R-INLA")

aghq_ecdf <- stats::ecdf(aghq_samples$`beta[1]`)
aghq_ecdf_df <- data.frame(x = grid, ecdf = aghq_ecdf(grid), method = "Gaussian", software = "TMB")

adam_ecdf_df <- compute_pdf_and_cdf(nodes = quad_laplace_marginals[[1]]$x, lps = quad_laplace_marginals[[1]]$lp_normalised, finegrid = grid) %>%
  select(x = x, ecdf = cdf) %>%
  mutate(method = "Laplace", software = "TMB")

tmbstan_ecdf <- stats::ecdf(rstan::extract(tmbstan, pars = "beta[1]")[[1]])
tmbstan_ecdf_df <- data.frame(x = grid, ecdf = tmbstan_ecdf(grid), method = "NUTS", software = "tmbstan")

inla_g_grid_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - inla_g_grid_ecdf_df$ecdf
inla_l_grid_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - inla_l_grid_ecdf_df$ecdf
aghq_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - aghq_ecdf_df$ecdf
adam_ecdf_df$ecdf_diff <- tmbstan_ecdf_df$ecdf - adam_ecdf_df$ecdf
tmbstan_ecdf_df$ecdf_diff <- 0

dummyA <- tmbstan_ecdf_df
dummyB <- tmbstan_ecdf_df
dummyA$software <- "TMB"
dummyB$software <- "R-INLA"

ecdf_df <- bind_rows(inla_g_grid_ecdf_df, inla_l_grid_ecdf_df, aghq_ecdf_df, adam_ecdf_df, dummyA, dummyB)

ecdf_df %>%
  rename("ECDF" = "ecdf", "ECDF difference to NUTS" = "ecdf_diff") %>%
  pivot_longer(cols = c("ECDF", "ECDF difference to NUTS"), names_to = "indicator", values_to = "value") %>%
  ggplot(aes(x = x, y = value, col = method)) +
  geom_line() +
  facet_wrap(indicator ~ software, scales = "free") +
  labs(x = "", y = "", col = "Method") +
  scale_color_manual(values = c("#56B4E9","#009E73", "#E69F00")) +
  theme_minimal()

ggsave("figures/naomi-aghq/intercept-comparison.png", h = 5, w = 6.25, bg = "white")

cbpalette <- c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#E69F00", "#E69F00")
bayesplot::color_scheme_set(rev(cbpalette))

tmbstan_summary <- summary(tmbstan)$summary
tmbstan_summary <- tmbstan_summary[1:(nrow(tmbstan_summary) - 1), ]
tmbstan_rhats <- bayesplot::rhat(tmbstan)
tmbstan_rhats <- tmbstan_rhats[1:(nrow(tmbstan_summary) - 1)]

round(min(tmbstan_summary[, "n_eff"])) #' 377
names(which.min(tmbstan_summary[, "n_eff"])) #' l_tau_nu

max(tmbstan_rhats) #' 1.006
names(which.max(tmbstan_rhats)) #' beta[3]

bayesplot::mcmc_trace(tmbstan, pars = c(names(which.min(tmbstan_summary[, "n_eff"])), names(which.max(tmbstan_rhats)))) +
  theme_minimal()

ggsave("figures/naomi-aghq/tmbstan-epil.png", h = 3, w = 6.25)

stan_summary <- summary(stan)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column("par")

stan_summary <- stan_summary %>%
  filter(par != "lp__")

stan_rhats <- bayesplot::rhat(stan)
stan_rhats <- stan_rhats[1:(nrow(stan_summary) - 1)]

round(min(stan_summary[, "n_eff"])) #' 437
stan_summary$par[which.min(stan_summary[, "n_eff"])] #' tau_nu

stan_summary$Rhat[order(stan_summary$Rhat, decreasing = TRUE)[2]] #' 1.008
stan_summary$par[order(stan_summary$Rhat, decreasing = TRUE)[2]] #' epsilon[18]

#' tau_nu had Rhat 1.009

bayesplot::mcmc_trace(stan, pars = c(stan_summary$par[which.min(stan_summary[, "n_eff"])], stan_summary$par[order(stan_summary$Rhat, decreasing = TRUE)[2]])) +
  theme_minimal()

ggsave("figures/naomi-aghq/stan-epil.png", h = 3, w = 6.25)


get_time <- function(t) {
  t$toc - t$tic
}

times <- list(tmb_g_eb_time, tmb_g_aghq_time, tmb_l_eb_time, tmb_l_aghq_time, inla_g_eb_time, inla_g_grid_time, inla_l_eb_time, inla_l_grid_time, tmbstan_time, stan_time)
times <- sapply(times, get_time)
methods <- c("Gaussian, EB", "Gaussian, AGHQ", "Laplace, EB", "Laplace, AGHQ", " Gaussian, EB", "Gaussian, grid", " Laplace, EB", "Laplace, grid", "NUTS", " NUTS")
softwares <- c(rep("TMB", times = 4), rep("R-INLA", times = 4), "tmbstan", "rstan")

time_df <- data.frame(time = times, method = methods, software = softwares)

ggplot(time_df, aes(x = forcats::fct_reorder(method, time), y = time, fill = software)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values = c("#0072B2", "#D55E00", "#E69F00", "#F0E442")) +
  labs(x = "", y = "Time taken (s)", fill = "Software") +
  coord_flip()

ggsave("figures/naomi-aghq/epil-time.png", h = 3.5, w = 6.25)

write_csv(time_df, "figures/naomi-aghq/epil-time.csv")
