library(TMB)
library(geostatsp)
library(aghq)
library(tmbstan)
library(stars)
library(tidyverse)
library(patchwork)

source("figures/naomi-aghq/functions.R")

data(loaloa, package = "geostatsp")
loaloa <- terra::vect(loaloa)
loaloa_sf <- sf::st_as_sf(loaloa)

cmr <- sf::st_read("figures/naomi-aghq/gadm41_CMR_2.json")
nga <- sf::st_read("figures/naomi-aghq/gadm41_NGA_2.json")
areas <- rbind(cmr, nga)
areas <- st_transform(areas, crs = st_crs(loaloa_sf))

compile("resources/naomi-aghq/loaloazip.cpp")
dyn.load(dynlib("resources/naomi-aghq/loaloazip"))

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

quad <- aghq::marginal_laplace_tmb(obj, 3, startingvalue = c(param$logkappa, param$logtau))

loaloa_sf <- mutate(loaloa_sf, p = y / N, zero = p == 0)

figA <- ggplot() +
  geom_sf(data = areas, col = "grey20") +
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

aghq_samples <- sample_marginal(quad, 100)
u_samples <- aghq_samples$samps[c(1:190), ]
v_samples <- aghq_samples$samps[c(192:381), ]
beta_samples <- aghq_samples$samps[c(191, 382), ]

theta_samples <- aghq_samples$theta

covariance_samples <- cbind(
  var = get_sigma(exp(theta_samples$logkappa), exp(theta_samples$logtau))^2,
  range = get_rho(exp(theta_samples$logkappa), exp(theta_samples$logtau)),
  shape = matern$nu
)

extent <- st_bbox(loaloa_sf) + 10000 * c(-1, -1, 1, 1)
grid <- st_as_stars(extent, dx = 10000)

vgm <- gstat::vgm(model = "Mat", range = 60000, shape = 1, psill = 1)
kriging_results <- gstat::krige(p ~ 1, loaloa_sf, grid, model = vgm)

ggplot() +
  geom_sf(data = areas, fill = NA) +
  geom_stars(data = kriging_results, aes(fill = var1.pred, x = x, y = y)) +
  geom_sf(data = loaloa_sf, shape = 4) +
  scale_fill_viridis_c() +
  theme_void() +
  labs(fill = "Prevalence")

nsim <- 100
phi_samples <- list()

for(i in 1:nsim) {
  vgm <- gstat::vgm(model = "Mat", range = covariance_samples[i, "range"], shape = 1, psill = covariance_samples[i, "var"])
  u_sf <- st_sf("u" = u_samples[, i], "geometry" = loaloa_sf$geometry)
  u_kriging_stars <- gstat::krige(u ~ 1, u_sf, grid, nmax = 30, model = vgm, nsim = 1)
  u_kriging_samples_sf <- st_as_sf(u_kriging_stars)
  u_kriging_samples <- st_drop_geometry(u_kriging_samples_sf)
  eta_phi_samples <- u_kriging_samples + rep(beta_samples[1, i], each = nrow(u_kriging_samples))
  phi_samples[[i]] <- apply(eta_phi_samples, 2, plogis)
}

phi_samples <- data.frame(dplyr::bind_cols(phi_samples))
phi_sf <- st_sf("phi" = rowMeans(phi_samples), "geometry" = u_kriging_samples_sf$geometry)

fig_suitability <- ggplot() +
  geom_sf(data = sf::st_intersection(phi_sf, areas), aes(fill = phi), alpha = 0.8, col = NA) +
  geom_sf(data = sf::st_crop(areas, sf::st_bbox(phi_sf)), fill = NA, col = "grey20", alpha = 0.8) +
  geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
  scale_fill_viridis_c(option = "A", direction = 1) +
  theme_void() +
  labs(fill = "Suitability", tag = "A")

rho_samples <- list()

for(i in 1:nsim) {
  vgm <- gstat::vgm(model = "Mat", range = covariance_samples[i, "range"], shape = 1, psill = covariance_samples[i, "var"])
  v_sf <- st_sf("v" = v_samples[, i], "geometry" = loaloa_sf$geometry)
  v_kriging_stars <- gstat::krige(v ~ 1, v_sf, grid, nmax = 30, model = vgm, nsim = 1)
  v_kriging_samples_sf <- st_as_sf(v_kriging_stars)
  v_kriging_samples <- st_drop_geometry(v_kriging_samples_sf)
  eta_rho_samples <- v_kriging_samples + rep(beta_samples[2, i], each = nrow(v_kriging_samples))
  rho_samples[[i]] <- apply(eta_rho_samples, 2, plogis)
}

rho_samples <- data.frame(dplyr::bind_cols(rho_samples))
rho_sf <- st_sf("rho" = rowMeans(rho_samples), "geometry" = v_kriging_samples_sf$geometry)

fig_prevalence <- ggplot() +
  geom_sf(data = sf::st_intersection(rho_sf, areas), aes(fill = rho), alpha = 0.8, col = NA) +
  geom_sf(data = sf::st_crop(areas, sf::st_bbox(rho_sf)), fill = NA, col = "grey20", alpha = 0.8) +
  geom_sf(data = loaloa_sf, shape = 4, size = 0.5, col = "grey30") +
  scale_fill_viridis_c(option = "D", direction = 1) +
  theme_void() +
  labs(fill = "Prevalence", tag = "B")

fig_suitability / fig_prevalence

ggsave("figures/naomi-aghq/conditional-simulation.png", h = 5, w = 6.25, bg = "white")

# Laplace marginals
random <- obj$env$random
N <- length(random)
W_names <- names(obj$env$par[random])

param[unique(W_names)] <- NULL
param$W_minus_i <- rnorm(N - 1)
param$W_i <- rnorm(1)

compile("resources/naomi-aghq/loaloazip_modified.cpp")
dyn.load(dynlib("resources/naomi-aghq/loaloazip_modified"))

compute_laplace_marginal <- function(i, quad) {
  dat$i <- i
  
  obj_i <- TMB::MakeADFun(
    data = dat,
    parameters = param,
    random = "W_minus_i",
    DLL = "loaloazip_modified",
    ADreport = FALSE,
    silent = TRUE
  )
  
  random_i <- obj_i$env$random
  mode_i <- quad$modesandhessians[["mode"]][[1]][-i]
  gg <- create_approx_grid(quad$modesandhessians, i = i, k = 5)
  out <- data.frame(index = i, par = W_names[i], x = mvQuad::getNodes(gg), w = mvQuad::getWeights(gg))
  
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

tictoc::tic()
temp <- compute_laplace_marginal(i = 1, quad)
time <- tictoc::toc()

(time$toc - time$tic) * N / 60 / 60

sample_adam <- function(M) {
  q <- runif(M)
  pdf_and_cdf <- compute_pdf_and_cdf(nodes = temp$x, temp$lp_normalised)
  s <- numeric(length(q))
  for(j in 1:length(q)) s[j] <- pdf_and_cdf$x[max(which(pdf_and_cdf$cdf < q[j]))]
  return(s)
}

samples <- sample_adam(1000)

ggplot(data.frame(samples), aes(x = samples)) +
  geom_histogram()

#' Try running tmbstan
nuts <- tmbstan::tmbstan(obj, chains = 1, warmup = 50, iter = 100)
nuts_samples <- as.data.frame(nuts)
u_samples <- t(nuts_samples[, c(1:190)])
v_samples <- t(nuts_samples[, c(192:381)])
beta_samples <- t(nuts_samples[, c(191, 382)])
theta_samples <- nuts_samples[, c(383, 384)]
