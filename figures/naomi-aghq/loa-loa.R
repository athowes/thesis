library(TMB)
library(geostatsp)
library(aghq)
library(tmbstan)
library(stars)

data(loaloa, package = "geostatsp")
loaloa <- terra::vect(loaloa)
loaloa_sf <- sf::st_as_sf(loaloa)

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

# Sampling from AGHQ fitted model

cmr <- sf::st_read("figures/naomi-aghq/gadm41_CMR_2.json")
nga <- sf::st_read("figures/naomi-aghq/gadm41_NGA_2.json")
areas <- rbind(cmr, nga)

loaloa_sf <- mutate(loaloa_sf, p = y / N)

grid <- raster::raster(loaloa_sf, nrows = 30, ncols = 30)
grid_stars <- stars::st_as_stars(grid)
b <- gstat::krige(p ~ 1, loaloa_sf, grid_stars)
b_sf <- st_as_sf(b)

ggplot() +
  geom_sf(data = areas, col = "grey70") +
  geom_sf(data = loaloa_sf, aes(col = p)) +
  # geom_sf(data = b_sf, aes(fill = var1.pred)) +
  scale_color_viridis_c() +
  theme_void() +
  lims(x = c(7.5, 16), y = c(3, 7)) +
  labs(x = "", y = "", col = "")

ggsave("figures/naomi-aghq/loa-loa-data.png", h = 4, w = 6.25)

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

# Resolution for spatial interpolations
resolution <- list(nrow = 50, ncol = 50)

options(useRandomFields = TRUE)

loaloa_u <- vect(loaloa_sf$geometry)
values(loaloa_u) <- as.data.frame(u_samples)

u_brick <- geostatsp::RFsimulate(
  model = covariance_samples,
  x = terra::rast(loaloa, nrows = resolution$nrow, ncols = resolution$ncol),
  data = loaloa_u
)

nrow(covariance_samples)

# Error in UnifyXT(x, y, z, T, grid = grid, distances = distances, dim = dim) : 
#   is.numeric(x) is not TRUE
# 'matrix,SpatRaster'

# Happens with no data so must be either model or x that is the issue!

loaloa_v <- vect(loaloa_sf$geometry)
values(loaloa_v) <- as.data.frame(v_samples)

v_brick <- geostatsp::RFsimulate(
  model = covariance_samples,
  x = terra::rast(loaloa, nrows = resolution$nrow, ncols = resolution$ncol),
  data = loaloa_v
)

eta_phi_brick <- beta_samples[1, ] + u_brick
eta_rho_brick <- beta_samples[2, ] + v_brick
