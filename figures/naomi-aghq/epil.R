library(tidyverse)
library(patchwork)
library(INLA)
library(TMB)
library(gt)

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

start <- Sys.time() 
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
end <- Sys.time()
tmb_g_eb_time <- end - start

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

ggsave("figures/naomi-aghq/hessian-matrix.png", h = 2.5, w = 6.25)

init <- c(param$l_tau_epsilon, param$l_tau_nu)
aghq <- aghq::marginal_laplace_tmb(obj, k = 3, startingvalue = init)

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

compile("resources/naomi-aghq/epil_modified.cpp")
dyn.load(dynlib("resources/naomi-aghq/epil_modified"))

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

# Laplace marginals and EB with TMB

init <- c(param$l_tau_epsilon, param$l_tau_nu)
eb <- aghq::marginal_laplace_tmb(obj, k = 1, startingvalue = init)
eb_laplace_marginals <- purrr::map(.x = 1:N, .f = compute_laplace_marginal, quad = eb, .progress = TRUE)
eb_laplace_marginals <- dplyr::bind_rows(eb_laplace_marginals)

# Laplace marginals and AGHQ with TMB

quad_laplace_marginals <- purrr::map(.x = 1:N, .f = compute_laplace_marginal, quad = aghq, .progress = TRUE)
quad_laplace_marginals <- dplyr::bind_rows(quad_laplace_marginals)

start <- Sys.time() 
tmbstan <- tmbstan::tmbstan(obj = obj, chains = 4, refresh = 0)
end <- Sys.time()
tmbstan_time <- end - start

df <- cbind(
  "Gaussian, EB (TMB)" = as.vector(t(data.frame(sd_out$par.random[1:6], sqrt(sd_out$diag.cov.random[1:6])))),
  "Gaussian, EB (R-INLA)" = as.vector(t(inla_g_eb$summary.fixed[1:6, 1:2])),
  "SL, EB (R-INLA)" = as.vector(t(inla_sl_eb$summary.fixed[1:6, 1:2])),
  "Laplace, EB (R-INLA)" = as.vector(t(inla_l_eb$summary.fixed[1:6, 1:2])),
  "Gaussian, grid (R-INLA)" = as.vector(t(inla_g_grid$summary.fixed[1:6, 1:2])),
  "SL, grid (R-INLA)" = as.vector(t(inla_sl_grid$summary.fixed[1:6, 1:2])),
  "Laplace, grid (R-INLA)" = as.vector(t(inla_l_grid$summary.fixed[1:6, 1:2])),
  "Gaussian, EB (TMB)" = as.vector(t(data.frame(sd_out$par.random[1:6], sqrt(sd_out$diag.cov.random[1:6])))),
  "Gaussian, AGHQ (TMB)" = as.vector(t(data.frame(mean = apply(aghq_samples[, 1:6], 2, mean), sd = apply(aghq_samples[, 1:6], 2, sd)))),
  "NUTS (tmbstan)" = as.vector(t(summary(tmbstan)$summary[1:6, c(1, 3)]))
) %>%
  as.data.frame() %>%
  round(digits = 3)

beta_i <- function(i) { c(paste0("beta[", i, "]"), paste0("sd(beta[", i, "])")) }
rownames(df) <- c(sapply(0:5, beta_i))

saveRDS(df, "resources/naomi-aghq/epil.rds")

time_df <- data.frame(
  time = c(tmb_g_eb_time, inla_g_eb_time, inla_sl_eb_time, inla_l_eb_time, inla_g_grid_time, inla_sl_grid_time, inla_l_grid_time, tmbstan_time),
  method = c(" Gaussian, EB", "Gaussian, EB", "SL, EB", "Laplace, EB", "Gaussian, grid", "SL, grid", "Laplace, grid", "NUTS"),
  software = c("TMB", rep("R-INLA", times = 6), "tmbstan")
)

ggplot(time_df, aes(x = forcats::fct_reorder(method, time), y = time, fill = software)) +
  geom_col(width = 0.7) +
  theme_minimal() +
  scale_fill_manual(values = c("#56B4E9","#009E73", "#E69F00")) +
  labs(x = "", y = "Time taken (s)", fill = "Software") +
  coord_flip()

ggsave("figures/naomi-aghq/epil-time.png", h = 3.5, w = 6.25)

beta0_inla_marginals <- bind_rows(
  data.frame(inla_g_grid$marginals.fixed$`(Intercept)`) %>%
    mutate(method = "Gaussian"),
  data.frame(inla_l_grid$marginals.fixed$`(Intercept)`) %>%
    mutate(method = "Laplace")
)

plot0 <- ggplot(data.frame(x = rstan::extract(tmbstan, pars = "beta[1]")[[1]]), aes(x = x)) +
  geom_histogram(aes(y = ..density..), alpha = 0.6, fill = "#E69F00", bins = 40) +
  theme_minimal() +
  lims(x = c(1.2, 2), y = c(0, 5.5)) +
  labs(x = "", y = "PDF")

inla_comparison <- plot0 +
  geom_line(data = beta0_inla_marginals, aes(x = x, y = y, col = method)) +
  scale_color_manual(values = c("#56B4E9", "#009E73")) +
  geom_col(data = data.frame(x = rep(inla_g_grid$summary.fixed["(Intercept)", ]$mean, 3), y = c(0, 0, 0), type = c("Gaussian (R-INLA)", "Laplace (R-INLA)", "NUTS")), aes(x = x, y = y, fill = type)) +
  scale_fill_manual(values = c("#56B4E9","#009E73", "#E69F00")) +
  labs(col = "", title = "", fill = "", tag = "A") +
  guides(col = "none", fill = guide_legend(override.aes = list(alpha = 1, shape = 15)))

beta0_quad_laplace_marginal <- filter(quad_laplace_marginals, index == 1)
beta0_quad_laplace_marginal <- compute_pdf_and_cdf(nodes = beta0_quad_laplace_marginal$x, lps = beta0_quad_laplace_marginal$lp_normalised)

tmb_comparison <- plot0 +
  geom_histogram(data = data.frame(x = aghq_samples$`beta[1]`), aes(y = ..density..), alpha = 0.6, fill = "#56B4E9", bins = 30) +
  geom_line(data = beta0_quad_laplace_marginal, aes(x = x, y = pdf), col = "#009E73") +
  geom_col(data = data.frame(x = rep(inla_g_grid$summary.fixed["(Intercept)", ]$mean, 3), y = c(0, 0, 0), type = c("Gaussian (TMB)", "Laplace (TMB)", "NUTS")), aes(x = x, y = y, fill = type)) +
  scale_fill_manual(values = c("#56B4E9","#009E73", "#E69F00")) +
  labs(col = "", title = "", fill = "", tag = "B") +
  guides(col = "none", fill = guide_legend(override.aes = list(alpha = 1, shape = 15)))

inla_comparison / tmb_comparison

ggsave("figures/naomi-aghq/intercept-comparison.png", h = 4, w = 6.25)
