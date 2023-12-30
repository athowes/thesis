library(TMB)
library(geostatsp)
library(aghq)
library(tmbstan)
library(stars)
library(tidyverse)
library(patchwork)
library(pals)

source("figures/naomi-aghq/functions.R")

#' The Loa loa data from geostatsp package
data(loaloa, package = "geostatsp")

#' It is useful to have it as sf for future
loaloa <- terra::vect(loaloa)
loaloa_sf <- sf::st_as_sf(loaloa)

#' Area files for the countries containing the villages sampled
cmr <- sf::st_read("figures/naomi-aghq/gadm41_CMR_2.json")
nga <- sf::st_read("figures/naomi-aghq/gadm41_NGA_2.json")
areas <- rbind(cmr, nga)
areas <- st_transform(areas, crs = st_crs(loaloa_sf))

#' Add column for the direct prevalence, and whether or not it's a zero
loaloa_sf <- mutate(loaloa_sf, p = y / N, zero = p == 0)

#' Initial EDA figure
figA <- ggplot() +
  geom_sf(data = sf::st_crop(areas, sf::st_bbox(loaloa_sf)), col = "grey20") +
  geom_sf(data = loaloa_sf, aes(col = p, size = N, shape = zero), alpha = 0.7) +
  scale_color_viridis_c() +
  scale_size(range = c(1, 4)) +
  theme_void() +
  labs(x = "", y = "", col = "Prevalence", size = "Sample size", shape = "Zero", tag = "A") +
  guides(
    col = guide_colourbar(order = 1),
    shape = guide_legend(override.aes = list(size = 2.5, col = "grey20"), order = 2),
    size = guide_legend(override.aes = list(shape = 16, col = "grey20"), order = 3)
  ) +
  theme(
    legend.direction = "vertical", 
    legend.box = "vertical",
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.width = unit(1, "line"),
    legend.key.height = unit(1, "line")
  )

figB <- ggplot(loaloa_sf, aes(x = p)) +
  geom_histogram(col = "grey60", fill = "grey80") +
  labs(x = "", y = "", tag = "B") +
  coord_fixed(ratio = 0.01) +
  theme_minimal()

figA / figB + plot_layout(heights = c(1.5, 1))

ggsave("figures/naomi-aghq/loa-loa-data.png", h = 5, w = 6.25, bg = "white")

#' TMB model
compile("figures/naomi-aghq/TMB/loaloazip.cpp")
dyn.load(dynlib("figures/naomi-aghq/TMB/loaloazip"))

Amat <- Diagonal(nrow(loaloa_sf))
Xmat <- cbind(rep(1, nrow(Amat)))
design <- bdiag(cbind(Amat, Xmat), cbind(Amat, Xmat))
y <- loaloa_sf$y
N <- loaloa_sf$N
n <- nrow(Xmat)
p <- ncol(Xmat) * 2
m <- ncol(Amat) * 2
Wd <- ncol(design)
stopifnot(Wd == m + p)

sigma_u <- 1
sigma_alpha <- 0.025
rho_u <- 2e05
rho_alpha <- 0.975

matern <- list()
matern$d <- 2
matern$nu <- 1

get_kappa <- function(sigma, rho) {
  sqrt(8 * matern$nu) / rho
}

get_tau <- function(sigma, rho) {
  sigma * get_kappa(sigma, rho)^(matern$nu) * sqrt(gamma(matern$nu + matern$d / 2) * (4 * pi)^(matern$d / 2) / gamma(matern$nu))
}

get_sigma <- function(kappa, tau) {
  tau / (kappa^(matern$nu) * sqrt(gamma(matern$nu + matern$d / 2) * (4 * pi)^(matern$d / 2) / gamma(matern$nu)))
}

get_rho <- function(kappa, tau) {
  sqrt(8 * matern$nu) / kappa
}

beta_prec <- 0.001

dat <- list(
  y = y,
  N = N,
  design = design,
  nu = matern$nu,
  rho_u = rho_u,
  rho_alpha = rho_alpha,
  sigma_u = sigma_u,
  sigma_alpha = sigma_alpha,
  D = raster::pointDistance(loaloa_sf, lonlat = FALSE),
  betaprec = beta_prec
)

starting_sig <- 1
starting_rho <- 4.22 * 1e04

set.seed(4564)

param <- list(
  W = rnorm(ncol(design)),
  logkappa = log(get_kappa(starting_sig, starting_rho)),
  logtau = log(get_tau(starting_sig, starting_rho))
)

obj <- MakeADFun(
  data = dat,
  parameters = param,
  random = "W",
  DLL = "loaloazip",
  ADreport = FALSE,
  silent = TRUE
)

#' Run AGHQ
quad <- aghq::marginal_laplace_tmb(obj, 3, startingvalue = c(param$logkappa, param$logtau))

# Large number of samples here to estimate beta well for fixing later
beta_samples <- sample_marginal(quad, 5000)
beta_fixed <- round(rowMeans(beta_samples$samps[c(191, 382), ]), 2) # 2.95 -1.99

#' NUTS can't run the whole model, so the approach is to run AGHQ, find beta, fix beta,
#' then run all inference procedures (Gaussian, Laplace, NUTS) again with fixed beta

#' The fixed beta version has a new TMB template
compile("figures/naomi-aghq/TMB/loaloazip_fixed.cpp")
dyn.load(dynlib("figures/naomi-aghq/TMB/loaloazip_fixed"))

#' The parameters and data need to be changed for the fixed template
param_new <- with(param, list(
  Urisk = W[192:381],
  Uzi = W[1:190],
  betarisk = beta_fixed[2],
  betazi = beta_fixed[1],
  logkappa = log(get_kappa(starting_sig, starting_rho)),
  logtau = log(get_tau(starting_sig, starting_rho))
))

dat_new <- within(dat,{
  A <- as(as(Amat, "dgCMatrix"), "dgTMatrix")
  X <- Xmat
})

dat_new$design <- NULL

#' The map option here fixes betazi and betarisk to their value in param_new
obj_fixed <- MakeADFun(
  data = dat_new,
  parameters = param_new,
  random = c("Urisk", "Uzi"),
  map = list(betarisk = factor(NA), betazi = factor(NA)),
  DLL = "loaloazip_fixed",
  ADreport = FALSE,
  silent = TRUE
)

tictoc::tic()
quad_fixed <- aghq::marginal_laplace_tmb(obj_fixed, 3, startingvalue = c(param$logkappa, param$logtau))
time_aghq <- tictoc::toc()

#' Small number of samples as illustration
aghq_samples <- sample_marginal(quad, 100)

#' Producing the output plot for AGHQ with no fixed beta
random_field_samples <- random_field_simulation(
  u_samples = aghq_samples$samps[c(1:190), ],
  v_samples = aghq_samples$samps[c(192:381), ],
  beta_samples = aghq_samples$samps[c(191, 382), ],
  theta_samples = aghq_samples$theta
)

phi_samples <- random_field_samples$phi_samples
rho_samples <- random_field_samples$rho_samples
phi_sf <- st_sf("phi" = rowMeans(phi_samples), "geometry" = random_field_samples$grid$geometry)
rho_sf <- st_sf("rho" = rowMeans(rho_samples), "geometry" = random_field_samples$grid$geometry)

plot_suitability <- function(phi_sf) {
  ggplot() +
    geom_sf(data = sf::st_intersection(phi_sf, areas), aes(fill = phi), alpha = 0.8, col = NA) +
    geom_sf(data = sf::st_crop(areas, sf::st_bbox(phi_sf)), fill = NA, col = "grey20", alpha = 0.8) +
    geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
    scale_fill_viridis_c(option = "A", direction = 1, labels = scales::percent) +
    theme_void() +
    labs(fill = "Suitability", tag = "A")
}

plot_prevalence <- function(rho_sf) {
  ggplot() +
    geom_sf(data = sf::st_intersection(rho_sf, areas), aes(fill = rho), alpha = 0.8, col = NA) +
    geom_sf(data = sf::st_crop(areas, sf::st_bbox(rho_sf)), fill = NA, col = "grey20", alpha = 0.8) +
    geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
    scale_fill_viridis_c(option = "D", direction = 1, labels = scales::percent) +
    theme_void() +
    labs(fill = "Prevalence", tag = "B")
}

plot_suitability(phi_sf) / plot_prevalence(rho_sf)

ggsave("figures/naomi-aghq/conditional-simulation.png", h = 5, w = 6.25, bg = "white")

#' Doing the same as above for AGHQ with a fixed beta
#' Use a larger number of samples, because this figure will be included in thesis
#' It's relatively slow to do this!

n_medium <- 500
aghq_fixed_samples <- sample_marginal(quad_fixed, n_medium)

random_field_samples_fixed <- random_field_simulation(
  u_samples = aghq_fixed_samples$samps[which(rownames(aghq_fixed_samples$samps) == "Uzi"), ],
  v_samples = aghq_fixed_samples$samps[which(rownames(aghq_fixed_samples$samps) == "Urisk"), ],
  beta_samples = t(data.frame(rep(param_new$betazi, n_medium), rep(param_new$betarisk, n_medium))),
  theta_samples = aghq_fixed_samples$theta,
  nsim = n_medium
)

phi_fixed_samples <- random_field_samples_fixed$phi_samples
rho_fixed_samples <- random_field_samples_fixed$rho_samples
phi_fixed_sf <- st_sf("phi" = rowMeans(phi_fixed_samples), "geometry" = random_field_samples_fixed$grid$geometry)
rho_fixed_sf <- st_sf("rho" = rowMeans(rho_fixed_samples), "geometry" = random_field_samples_fixed$grid$geometry)

plot_suitability(phi_fixed_sf) / plot_prevalence(rho_fixed_sf)

ggsave("figures/naomi-aghq/conditional-simulation-fixed.png", h = 5, w = 6.25, bg = "white")

# Laplace marginals with fixed beta
set.seed(4564)

W_names <- names(param_new)[1:2]
W_lengths <- lengths(param_new[unique(W_names)])
N <- sum(W_lengths)
W_starts <- cumsum(W_lengths) - W_lengths

dat_new$W_starts <- as.numeric(W_starts)
dat_new$W_lengths <- as.numeric(W_lengths)

W_init <- c(param_new$Urisk, param_new$Uzi)
param_new[unique(W_names)] <- NULL

compile("figures/naomi-aghq/TMB/loaloazip_fixed_modified.cpp")
dyn.load(dynlib("figures/naomi-aghq/TMB/loaloazip_fixed_modified"))

compute_laplace_marginal <- function(i, quad) {
  dat_new$i <- i
  
  param_new$W_minus_i <- W_init[-i]
  param_new$W_i <- W_init[i]
  
  obj_fixed_i <- TMB::MakeADFun(
    data = dat_new,
    parameters = param_new,
    random = "W_minus_i",
    map = list(betazi = factor(NA), betarisk = factor(NA)),
    DLL = "loaloazip_fixed_modified",
    ADreport = FALSE,
    silent = FALSE
  )
  
  random_i <- obj_fixed_i$env$random
  mode_i <- quad$modesandhessians[["mode"]][[1]][-i]
  
  gg <- create_approx_grid(quad$modesandhessians, i = i, k = 5)
  out <- data.frame(index = i, par = "W", x = mvQuad::getNodes(gg), w = mvQuad::getWeights(gg))
  
  theta_names <- make.unique(names(obj$par), sep = "")
  
  .g <- function(x) {
    lp <- vector(mode = "numeric", length = nrow(quad$modesandhessians))
    
    for(z in 1:nrow(quad$modesandhessians)) {
      theta <- as.numeric(quad$modesandhessians[z, theta_names])
      obj_fixed_i$env$last.par[random_i] <- quad$modesandhessians[z, "mode"][[1]][-dat_new$i]
      lp[z] <- as.numeric(- obj_fixed_i$fn(c(theta, x))) # Look at order here, usually it's theta, x (?)
    }
    
    return(logSumExpWeights(lp, w = quad$normalized_posterior$nodesandweights$weights))
  }
  
  out$lp <- purrr::map_dbl(out$x, .g)
  lognormconst <- logSumExpWeights(out$lp, out$w)
  out$lp_normalised <- out$lp - lognormconst
  
  return(out)
}

tictoc::tic()
test_laplace <- compute_laplace_marginal(i = 1, quad = quad_fixed)
time_test <- tictoc::toc()

#' This would take around 3 hours to run for 1:N
# tictoc::tic()
# quad_fixed_laplace_marginals <- purrr::map(.x = 1:N, .f = compute_laplace_marginal, quad = quad_fixed, .progress = TRUE)
# time_laplace <- tictoc::toc()
#
# saveRDS(quad_fixed_laplace_marginals, "figures/naomi-aghq/loa-loa_laplace.rds")

quad_fixed_laplace_marginals <- readRDS("figures/naomi-aghq/loa-loa-laplace.rds")

#' Function to sample from the Laplace marginals using CDF inversion
sample_adam <- function(i, quad, M) {
  q <- runif(M)
  pdf_and_cdf <- compute_pdf_and_cdf(nodes = quad[[i]]$x, quad[[i]]$lp_normalised)
  s <- numeric(length(q))
  for(j in 1:length(q)) s[j] <- pdf_and_cdf$x[max(which(pdf_and_cdf$cdf < q[j]))]
  return(s)
}

#' Laplace marginals maps
#' Also use n_medium simulations here 

samples_adam <- lapply(1:N, sample_adam, M = n_medium, quad = quad_fixed_laplace_marginals)
samples_adam <- do.call(rbind, samples_adam)

random_field_samples_fixed_adam <- random_field_simulation(
  u_samples = samples_adam[c(191:380), ],
  v_samples = samples_adam[c(1:190), ],
  beta_samples = t(data.frame(rep(param_new$betazi, n_medium), rep(param_new$betarisk, n_medium))),
  theta_samples = aghq_fixed_samples$theta,
  nsim = n_medium
)

phi_fixed_adam_samples <- random_field_samples_fixed_adam$phi_samples
rho_fixed_adam_samples <- random_field_samples_fixed_adam$rho_samples
phi_fixed_adam_sf <- st_sf("phi" = rowMeans(phi_fixed_adam_samples), "geometry" = random_field_samples_fixed_adam$grid$geometry)
rho_fixed_adam_sf <- st_sf("rho" = rowMeans(rho_fixed_adam_samples), "geometry" = random_field_samples_fixed_adam$grid$geometry)

plot_suitability(phi_fixed_adam_sf) / plot_prevalence(rho_fixed_adam_sf)

ggsave("figures/naomi-aghq/conditional-simulation-adam-fixed.png", h = 5, w = 6.25, bg = "white")

#' Run tmbstan
#' Alex writes that "the sampler converged in just over 19 hours, for 10,000 iterations"
#' I found that getting 10,000 iterations (as below) required 23 hours, though I did use my laptop

# tictoc::tic()
# nuts <- tmbstan::tmbstan(obj_fixed, chains = 4, iter = 5000)
# time_nuts <- tictoc::toc()
# 
# (time$toc - time$tic) / 60 / 60
# 
# saveRDS(time_nuts, file = "figures/naomi-aghq/nuts-big-time.rds")
# saveRDS(nuts, file = "figures/naomi-aghq/nuts-big.rds")

nuts <- readRDS("figures/naomi-aghq/nuts-big.rds")

#' Diagnostics look good
bayesplot::color_scheme_set(rev(c("#56B4E9", "#009E73", "#E69F00", "#F0E442", "#E69F00", "#E69F00")))

nuts_summary <- head(summary(nuts)$summary, -1)
nuts_rhats <- head(bayesplot::rhat(nuts), -1)

round(min(nuts_summary[, "n_eff"])) #' 677
names(which.min(nuts_summary[, "n_eff"])) #' Uzi[174]

max(nuts_rhats) #' 1.004941
names(which.max(nuts_rhats)) #' Uzi[88]

bayesplot::mcmc_trace(nuts, pars = c(names(which.min(nuts_summary[, "n_eff"])), names(which.max(nuts_rhats)))) +
  theme_minimal()

ggsave("figures/naomi-aghq/nuts-loa-loa.png", h = 3, w = 6.25)

nuts_df <- as.data.frame(nuts)

#' I am going to thin the nuts_df so that there are 500 samples to pass into conditional simulation
#' To do this, I'll take the final 5000 iterations, and keep every 10th one
nuts_thin_df <- nuts_df[5001:10000, ]
nuts_thin_df <- nuts_thin_df[seq(1, 5000, by = 10), ]

nuts_random_field_samples <- random_field_simulation(
  u_samples = t(nuts_thin_df[, c(191:380)]),
  v_samples = t(nuts_thin_df[, c(1:190)]),
  beta_samples = t(data.frame(rep(param_new$betazi, n_medium), rep(param_new$betarisk, n_medium))),
  theta_samples = nuts_thin_df[, c(381, 382)],
  nsim = n_medium
)

phi_nuts_samples <- nuts_random_field_samples$phi_samples
rho_nuts_samples <- nuts_random_field_samples$rho_samples
phi_nuts_sf <- st_sf("phi" = rowMeans(phi_nuts_samples), "geometry" = nuts_random_field_samples$grid$geometry)
rho_nuts_sf <- st_sf("rho" = rowMeans(rho_nuts_samples), "geometry" = nuts_random_field_samples$grid$geometry)

plot_suitability(phi_nuts_sf) / plot_prevalence(rho_nuts_sf)

phi_diff_sf <- st_sf(
  "Gaussian" = rowMeans(phi_fixed_samples) - rowMeans(phi_nuts_samples),
  "Laplace" = rowMeans(phi_fixed_adam_samples) - rowMeans(phi_nuts_samples),
  "geometry" = nuts_random_field_samples$grid$geometry
) %>%
  pivot_longer(cols = c("Gaussian", "Laplace"), names_to = "method", values_to = "value")

rho_diff_sf <- st_sf(
  "Gaussian" = rowMeans(rho_fixed_samples) - rowMeans(rho_nuts_samples),
  "Laplace" = rowMeans(rho_fixed_adam_samples) - rowMeans(rho_nuts_samples),
  "geometry" = nuts_random_field_samples$grid$geometry
) %>%
  pivot_longer(cols = c("Gaussian", "Laplace"), names_to = "method", values_to = "value")

abs(phi_diff_sf$value) %>% max()

fig_phi_diff <- ggplot() +
  geom_sf(data = sf::st_intersection(phi_diff_sf, areas), aes(fill = value), alpha = 0.8, col = NA) +
  geom_sf(data = sf::st_crop(areas, sf::st_bbox(phi_diff_sf)), fill = NA, col = "grey20", alpha = 0.8) +
  geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
  facet_wrap(~ method, ncol = 1) +
  scale_fill_gradientn(colours = pals::ocean.balance(100), labels = scales::percent, limits = c(-0.086, 0.086)) +
  theme_void() +
  labs(fill = "Suitability\ndifference\nto NUTS")

fig_phi_diff

ggsave("figures/naomi-aghq/conditional-simulation-phi-diff-fixed.png", h = 4, w = 6.25, bg = "white")

abs(rho_diff_sf$value) %>% max()

fig_rho_diff <- ggplot() +
  geom_sf(data = sf::st_intersection(rho_diff_sf, areas), aes(fill = value), alpha = 0.8, col = NA) +
  geom_sf(data = sf::st_crop(areas, sf::st_bbox(rho_diff_sf)), fill = NA, col = "grey20", alpha = 0.8) +
  geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
  facet_wrap(~ method, ncol = 1) +
  scale_fill_gradientn(colours = pals::ocean.curl(100), labels = scales::percent, limits = c(-0.066, 0.066)) +
  theme_void() +
  labs(fill = "Prevalence\ndifference\nto NUTS")

fig_rho_diff

ggsave("figures/naomi-aghq/conditional-simulation-rho-diff-fixed.png", h = 4, w = 6.25, bg = "white")

#' Point estimate differences
samples_adam <- lapply(1:N, sample_adam, M = 5000, quad = quad_fixed_laplace_marginals)
aghq_fixed_samples <- sample_marginal(quad_fixed, 5000)

df <- bind_rows(
  data.frame(method = "Laplace", mean = sapply(samples_adam, mean), sd = sapply(samples_adam, sd), index = 1:N),
  data.frame(method = "Gaussian", mean = apply(aghq_fixed_samples$samps[1:N, ], 1, mean), sd = apply(aghq_fixed_samples$samps[1:N, ], 1, sd), index = 1:N),
  data.frame(method = "NUTS", mean = nuts_summary[1:380, "mean"], sd = nuts_summary[1:380, "sd"], index = 1:N) %>% `rownames<-`(NULL)
)

df_plot <- df %>%
  pivot_longer(cols = c("mean", "sd"), names_to = "indicator", values_to = "value") %>%
  pivot_wider(id_cols = c("index", "indicator"), names_from = "method") %>%
  pivot_longer(cols = c("Laplace", "Gaussian"), names_to = "method", values_to = "estimate") %>%
  mutate(
    indicator = fct_recode(indicator, "Posterior mean" = "mean", "Posterior SD" = "sd"),
    diff_abs = estimate - NUTS,
    diff_pct = diff_abs / (sign(NUTS) * pmax(0.25, abs(NUTS)))
  ) %>%
  pivot_longer(cols = c("diff_abs", "diff_pct"), names_to = "metric", values_to = "value")

node_diff_df <- df %>%
  pivot_longer(cols = c("mean", "sd"), names_to = "indicator", values_to = "value") %>%
  pivot_wider(id_cols = c("index", "indicator"), names_from = "method") %>%
  mutate(
    indicator = fct_recode(indicator, "Posterior mean" = "mean", "Posterior SD" = "sd"),
    diff_laplace = abs(Laplace - NUTS),
    diff_gaussian = abs(Gaussian - NUTS),
    diff = diff_gaussian - diff_laplace,
  )

fig_abs <- df_plot %>%
  filter(metric == "diff_abs") %>%
  ggplot(aes(x = estimate, y = value)) +
  geom_point(size = 1.5, shape = 1, alpha = 0.6) +
  facet_grid(indicator ~ method) +
  geom_abline(intercept = 0, slope = 0, col = "#E69F00", linetype = "dashed") +
  labs(y = "Absolute difference to NUTS", x = "Estimate", tag = "A") +
  theme_minimal()

fig_abs_hist <- ggplot(node_diff_df, aes(x = diff)) +
  geom_histogram(col = "grey60", fill = "grey80", bins = 30) +
  facet_wrap(~ indicator, scales = "free") +
  labs(x = "Difference between Laplace and Gaussian", y = "", tag = "B") +
  theme_minimal()

fig_abs / fig_abs_hist + plot_layout(heights = c(1, 0.5))
  
ggsave("figures/naomi-aghq/loa-loa-mean-sd-abs.png", h = 7, w = 6.25, bg = "white")

fig_pct <- df_plot %>%
  filter(metric == "diff_pct") %>%
  ggplot(aes(x = estimate, y = value)) +
  geom_point(size = 1.5, shape = 1, alpha = 0.6) +
  facet_grid(indicator ~ method) +
  geom_abline(intercept = 0, slope = 0, col = "#E69F00", linetype = "dashed") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percent difference to NUTS", x = "Estimate") +
  theme_minimal()

fig_pct

ggsave("figures/naomi-aghq/loa-loa-mean-sd-pct.png", h = 6, w = 6.25, bg = "white")

#' Index 184 has the largest difference in mean comparing Laplace and Gaussian marginals
node_diff_df %>%
  arrange(-diff)

grid <- seq(
  from = round(min(rstan::extract(nuts, pars = "Urisk[184]")[[1]]) - 0.1, digits = 2),
  to = round(max(rstan::extract(nuts, pars = "Urisk[184]")[[1]]) + 0.1, digits = 2),
  length.out = 1000
)

gaussian_ecdf <- stats::ecdf(aghq_fixed_samples$samps[184, ])
gaussian_ecdf_df <- data.frame(x = grid, ecdf = aghq_ecdf(grid), method = "Gaussian")

laplace_ecdf <- stats::ecdf(samples_adam[[184]])
laplace_ecdf_df <- data.frame(x = grid, ecdf = laplace_ecdf(grid), method = "Laplace")

nuts_ecdf <- stats::ecdf(rstan::extract(nuts, pars = "Urisk[184]")[[1]])
nuts_ecdf_df <- data.frame(x = grid, ecdf = nuts_ecdf(grid), method = "NUTS")

gaussian_ecdf_df$ecdf_diff <- nuts_ecdf_df$ecdf - gaussian_ecdf_df$ecdf
laplace_ecdf_df$ecdf_diff <- nuts_ecdf_df$ecdf - laplace_ecdf_df$ecdf
nuts_ecdf_df$ecdf_diff <- 0

ecdf_df <- bind_rows(gaussian_ecdf_df, laplace_ecdf_df, nuts_ecdf_df)

ecdf_df %>%
  rename("ECDF" = "ecdf", "ECDF difference to NUTS" = "ecdf_diff") %>%
  pivot_longer(cols = c("ECDF", "ECDF difference to NUTS"), names_to = "indicator", values_to = "value") %>%
  ggplot(aes(x = x, y = value, col = method)) +
  geom_line() +
  facet_wrap(~ indicator, scales = "free", ncol = 2) +
  labs(x = "", y = "", col = "Method") +
  scale_color_manual(values = c("#56B4E9","#009E73", "#E69F00")) +
  theme_minimal()

ggsave("figures/naomi-aghq/loa-loa-worst-node.png", h = 3.5, w = 6.25, bg = "white")

time_nuts <- readRDS(file = "figures/naomi-aghq/nuts-big-time.rds")

time_df <- data.frame(
  "time" = c(time_aghq$toc - time_aghq$tic, (time_test$toc - time_test$tic) * N, time_nuts$toc - time_nuts$tic),
  "method" = c("Gaussian, AGHQ", "Laplace, AGHQ", "NUTS"),
  "software" = c("TMB", "TMB", "tmbstan")
  )

readr::write_csv(time_df, "figures/naomi-aghq/loa-loa-time.csv")

time_df %>%
  mutate(mins = time / 60) %>%
  ggplot(aes(x = method, y = mins, fill = software)) +
  geom_col() +
  theme_minimal() +
  scale_fill_manual(values = c("#D55E00", "#CC79A7")) +
  labs(x = "", y = "Time taken (m)", fill = "Software") +
  coord_flip()

ggsave("figures/naomi-aghq/loa-loa-time.png", h = 3, w = 6.25)

#' Section in conclusions about using a bigger grid being possibly better than more accurate marginals
#' Comparison between k = 3 AGHQ with Gaussian marginals and k = 7 AGHQ with Gaussian marginals

tictoc::tic()
big_quad_fixed <- aghq::marginal_laplace_tmb(obj_fixed, 7, startingvalue = c(param$logkappa, param$logtau))
time_big_aghq <- tictoc::toc()

aghq_big_fixed_samples <- sample_marginal(big_quad_fixed, n_medium)

random_field_samples_big_fixed <- random_field_simulation(
  u_samples = aghq_big_fixed_samples$samps[which(rownames(aghq_big_fixed_samples$samps) == "Uzi"), ],
  v_samples = aghq_big_fixed_samples$samps[which(rownames(aghq_big_fixed_samples$samps) == "Urisk"), ],
  beta_samples = t(data.frame(rep(param_new$betazi, n_medium), rep(param_new$betarisk, n_medium))),
  theta_samples = aghq_big_fixed_samples$theta,
  nsim = n_medium
)

phi_big_fixed_samples <- random_field_samples_big_fixed$phi_samples
rho_big_fixed_samples <- random_field_samples_big_fixed$rho_samples
phi_big_fixed_sf <- st_sf("phi" = rowMeans(phi_big_fixed_samples), "geometry" = random_field_samples_big_fixed$grid$geometry)
rho_big_fixed_sf <- st_sf("rho" = rowMeans(rho_big_fixed_samples), "geometry" = random_field_samples_big_fixed$grid$geometry)

phi_diff_sf <- st_sf(
  "3" = rowMeans(phi_fixed_samples) - rowMeans(phi_nuts_samples),
  "7" = rowMeans(phi_big_fixed_samples) - rowMeans(phi_nuts_samples),
  "geometry" = nuts_random_field_samples$grid$geometry
) %>%
  pivot_longer(cols = c("3", "7"), names_to = "method", values_to = "value")

rho_diff_sf <- st_sf(
  "3" = rowMeans(rho_fixed_samples) - rowMeans(rho_nuts_samples),
  "7" = rowMeans(rho_big_fixed_samples) - rowMeans(rho_nuts_samples),
  "geometry" = nuts_random_field_samples$grid$geometry
) %>%
  pivot_longer(cols = c("3", "7"), names_to = "method", values_to = "value")

abs(phi_diff_sf$value) %>% max()

fig_phi_diff_k <- ggplot() +
  geom_sf(data = sf::st_intersection(phi_diff_sf, areas), aes(fill = value), alpha = 0.8, col = NA) +
  geom_sf(data = sf::st_crop(areas, sf::st_bbox(phi_diff_sf)), fill = NA, col = "grey20", alpha = 0.8) +
  geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
  facet_wrap(~ method, ncol = 1) +
  scale_fill_gradientn(colours = pals::ocean.balance(100), labels = scales::percent, limits = c(-0.096, 0.096)) +
  theme_void() +
  labs(fill = "Suitability\ndifference\nto NUTS", tag = "A") +
  theme(
    legend.position = "bottom"
  )

fig_rho_diff_k <- ggplot() +
  geom_sf(data = sf::st_intersection(rho_diff_sf, areas), aes(fill = value), alpha = 0.8, col = NA) +
  geom_sf(data = sf::st_crop(areas, sf::st_bbox(rho_diff_sf)), fill = NA, col = "grey20", alpha = 0.8) +
  geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
  facet_wrap(~ method, ncol = 1) +
  scale_fill_gradientn(colours = pals::ocean.curl(100), labels = scales::percent, limits = c(-0.053, 0.053), breaks = c(-0.05, 0, 0.05)) +
  theme_void() +
  labs(fill = "Prevalence\ndifference\nto NUTS", tag = "B") +
  theme(
    legend.position = "bottom"
  )

fig_phi_diff_k + fig_rho_diff_k

ggsave("figures/naomi-aghq/conditional-simulation-diff-k-fixed.png", h = 4, w = 6.25, bg = "white")
