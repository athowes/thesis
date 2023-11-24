library(TMB)
library(geostatsp)
library(aghq)
library(tmbstan)
library(stars)
library(tidyverse)
library(patchwork)

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

quad_fixed <- aghq::marginal_laplace_tmb(obj_fixed, 3, startingvalue = c(param$logkappa, param$logtau))

#' Small number of samples
aghq_samples <- sample_marginal(quad, 100)
aghq_fixed_samples <- sample_marginal(quad_fixed, 100)

#' Function to run conditional simulation of a random field. It is all done using
#' samples, so each simulation uses one sample of u, v, beta and theta. This requires
#' lapply over 1:sim and running gstat::krige for each iteration
random_field_simulation <- function(u_samples, v_samples, beta_samples, theta_samples, nsim = 100) {
  covariance_samples <- cbind(
    var = get_sigma(exp(theta_samples$logkappa), exp(theta_samples$logtau))^2,
    range = get_rho(exp(theta_samples$logkappa), exp(theta_samples$logtau)),
    shape = matern$nu
  )
  
  extent <- st_bbox(loaloa_sf) + 10000 * c(-1, -1, 1, 1)
  grid <- st_as_stars(extent, dx = 10000)

  phi_samples <- lapply(1:nsim, function(i) {
    vgm <- gstat::vgm(model = "Mat", range = covariance_samples[i, "range"], shape = 1, psill = covariance_samples[i, "var"])
    u_sf <- st_sf("u" = u_samples[, i], "geometry" = loaloa_sf$geometry)
    u_kriging_stars <- gstat::krige(u ~ 1, u_sf, grid, nmax = 30, model = vgm, nsim = 1)
    u_kriging_samples_sf <- st_as_sf(u_kriging_stars)
    u_kriging_samples <- st_drop_geometry(u_kriging_samples_sf)
    eta_phi_samples <- u_kriging_samples + rep(beta_samples[1, i], each = nrow(u_kriging_samples))
    apply(eta_phi_samples, 2, plogis)
  })
  
  phi_samples <- data.frame(dplyr::bind_cols(phi_samples))
  names(phi_samples) <- paste0("sim", 1:nsim)

  rho_samples <- lapply(1:nsim, function(i) {
    vgm <- gstat::vgm(model = "Mat", range = covariance_samples[i, "range"], shape = 1, psill = covariance_samples[i, "var"])
    v_sf <- st_sf("v" = v_samples[, i], "geometry" = loaloa_sf$geometry)
    v_kriging_stars <- gstat::krige(v ~ 1, v_sf, grid, nmax = 30, model = vgm, nsim = 1)
    v_kriging_samples_sf <- st_as_sf(v_kriging_stars)
    v_kriging_samples <- st_drop_geometry(v_kriging_samples_sf)
    eta_rho_samples <- v_kriging_samples + rep(beta_samples[2, i], each = nrow(v_kriging_samples))
    apply(eta_rho_samples, 2, plogis)
  })
  
  rho_samples <- data.frame(dplyr::bind_cols(rho_samples))
  names(rho_samples) <- paste0("sim", 1:nsim)
  
  # Samples from the phi and rho random fields on a grid
  return(list("phi_samples" = phi_samples, "rho_samples" = rho_samples, "grid" = st_as_sf(grid)))
}

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
    scale_fill_viridis_c(option = "A", direction = 1) +
    theme_void() +
    labs(fill = "Suitability", tag = "A")
}

plot_prevalence <- function(rho_sf) {
  ggplot() +
    geom_sf(data = sf::st_intersection(rho_sf, areas), aes(fill = rho), alpha = 0.8, col = NA) +
    geom_sf(data = sf::st_crop(areas, sf::st_bbox(rho_sf)), fill = NA, col = "grey20", alpha = 0.8) +
    geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
    scale_fill_viridis_c(option = "D", direction = 1) +
    theme_void() +
    labs(fill = "Prevalence", tag = "B")
}

plot_suitability(phi_sf) / plot_prevalence(rho_sf)

ggsave("figures/naomi-aghq/conditional-simulation.png", h = 5, w = 6.25, bg = "white")

#' Doing the same as above for AGHQ with a fixed beta
random_field_samples_fixed <- random_field_simulation(
  u_samples = aghq_fixed_samples$samps[which(rownames(aghq_fixed_samples$samps) == "Uzi"), ],
  v_samples = aghq_fixed_samples$samps[which(rownames(aghq_fixed_samples$samps) == "Urisk"), ],
  beta_samples = t(data.frame(rep(param_new$betarisk, 100), rep(param_new$betazi, 100))),
  theta_samples = aghq_fixed_samples$theta
)

phi_fixed_samples <- random_field_samples_fixed$phi_samples
rho_fixed_samples <- random_field_samples_fixed$rho_samples
phi_fixed_sf <- st_sf("phi" = rowMeans(phi_samples), "geometry" = random_field_samples_fixed$grid$geometry)
rho_fixed_sf <- st_sf("rho" = rowMeans(rho_samples), "geometry" = random_field_samples_fixed$grid$geometry)

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
    silent = TRUE
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
      lp[z] <- as.numeric(- obj_fixed_i$fn(c(x, theta)))
    }
    
    return(logSumExpWeights(lp, w = quad$normalized_posterior$nodesandweights$weights))
  }
  
  out$lp <- purrr::map_dbl(out$x, .g)
  lognormconst <- logSumExpWeights(out$lp, out$w)
  out$lp_normalised <- out$lp - lognormconst
  
  return(out)
}

# i <- 1
# dat_new$i <- i
# 
# param_new$W_minus_i <- W_init[-i]
# param_new$W_i <- W_init[i]
# 
# obj_fixed_i <- TMB::MakeADFun(
#   data = dat_new,
#   parameters = param_new,
#   random = "W_minus_i",
#   map = list(betazi = factor(NA), betarisk = factor(NA)),
#   DLL = "loaloazip_fixed_modified",
#   ADreport = FALSE,
#   silent = FALSE
# )
# 
# random_i <- obj_fixed_i$env$random
# mode_i <- quad_fixed$modesandhessians[["mode"]][[1]][-i]
# gg <- create_approx_grid(quad_fixed$modesandhessians, i = i, k = 5)
# out <- data.frame(index = i, par = "W", x = mvQuad::getNodes(gg), w = mvQuad::getWeights(gg))
# 
# theta_names <- make.unique(names(obj$par), sep = "")
# 
# .g <- function(x) {
#   lp <- vector(mode = "numeric", length = nrow(quad_fixed$modesandhessians))
#   
#   for(z in 1:nrow(quad_fixed$modesandhessians)) {
#     theta <- as.numeric(quad_fixed$modesandhessians[z, theta_names])
#     obj_fixed_i$env$last.par[random_i] <- quad_fixed$modesandhessians[z, "mode"][[1]][-dat_new$i]
#     lp[z] <- as.numeric(- obj_fixed_i$fn(c(x, theta)))
#   }
#   
#   return(logSumExpWeights(lp, w = quad$normalized_posterior$nodesandweights$weights))
# }

#' It's saying it's doing the marginal for u[1] not v[1] (!?)
#' #' v[1] is -3.21
# .g(-0.242215130)
# .g(-3)

tictoc::tic()
test_laplace <- compute_laplace_marginal(i = 1, quad = quad_fixed)
time <- tictoc::toc()

(time$toc - time$tic) * N / 60 / 60

#' This would take around 3 hours to run
quad_fixed_laplace_marginals <- purrr::map(.x = 1:5, .f = compute_laplace_marginal, quad = quad_fixed, .progress = TRUE)

#' Function to sample from the Laplace marginals using CDF inversion
sample_adam <- function(i, M) {
  q <- runif(M)
  pdf_and_cdf <- compute_pdf_and_cdf(nodes = quad_laplace_marginals[[i]]$x, quad_laplace_marginals[[i]]$lp_normalised)
  s <- numeric(length(q))
  for(j in 1:length(q)) s[j] <- pdf_and_cdf$x[max(which(pdf_and_cdf$cdf < q[j]))]
  return(s)
}

samples_adam <- lapply(1:5, sample_adam, M = 1000)
aghq_fixed_samples <- sample_marginal(quad_fixed, 1000)

df <- bind_rows(
  data.frame(method = "laplace", mean = sapply(samples_adam, mean), sd = sapply(samples_adam, sd), index = 1:5),
  data.frame(method = "gaussian", mean = apply(aghq_fixed_samples$samps[1:5, ], 1, mean), sd = apply(aghq_fixed_samples$samps[1:5, ], 1, sd), index = 1:5)
)

df %>%
  pivot_longer(cols = c("mean", "sd"), names_to = "indicator", values_to = "value") %>%
  pivot_wider(id_cols = c("index", "indicator"), names_from = "method") %>%
  mutate(
    indicator = fct_recode(indicator, "Posterior mean" = "mean", "Posterior SD" = "sd")
  ) %>%
  ggplot(aes(x = gaussian, y = (laplace - gaussian) / laplace)) +
    geom_point(size = 1.5, shape = 1) +
    facet_wrap(~ indicator, scales = "free_x") +
    geom_abline(intercept = 0, slope = 0, col = "#009E73", linetype = "dashed") +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Gaussian estimate", y = "Percentage difference to Laplace") +
    theme_minimal()

ggsave("figures/naomi-aghq/loa-loa-mean-sd.png", h = 3.5, w = 6.25, bg = "white")

#' Try running tmbstan
nuts <- tmbstan::tmbstan(obj, chains = 1, warmup = 50, iter = 100)

nuts_random_field_samples <- random_field_simulation(t(as.data.frame(nuts)), as.data.frame(nuts)[, c(383, 384)], nsim = 50)
nuts_phi_samples <- nuts_random_field_samples$phi_samples
nuts_rho_samples <- nuts_random_field_samples$rho_samples
nuts_phi_sf <- st_sf("phi" = rowMeans(nuts_phi_samples), "geometry" = nuts_random_field_samples$grid$geometry)
nuts_rho_sf <- st_sf("rho" = rowMeans(nuts_rho_samples), "geometry" = nuts_random_field_samples$grid$geometry)

plot_suitability(nuts_phi_sf)
plot_prevalence(nuts_rho_sf)
