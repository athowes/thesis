library(TMB)
library(geostatsp)
library(aghq)
library(tmbstan)
library(stars)

data(loaloa, package = "geostatsp")
loaloa <- terra::vect(loaloa)
loaloa_sf <- sf::st_as_sf(loaloa)

# Set the resolution for the spatial interpolations
reslist <- list(nrow = 50, ncol = 100)

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

maternconstants <- list()
maternconstants$d <- 2
maternconstants$nu <- 1

get_kappa <- function(sigma, rho) {
  sqrt(8 * maternconstants$nu) / rho
}

get_tau <- function(sigma, rho) {
  sigma * get_kappa(sigma, rho)^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d / 2) * (4 * pi)^(maternconstants$d / 2) / gamma(maternconstants$nu))
}

get_sigma <- function(kappa, tau) {
  tau / (kappa^(maternconstants$nu) * sqrt(gamma(maternconstants$nu + maternconstants$d / 2) * (4 * pi)^(maternconstants$d / 2) / gamma(maternconstants$nu)))
}

get_rho <- function(kappa, tau) {
  sqrt(8 * maternconstants$nu) / kappa
}

beta_prec <- 0.001

dat <- list(
  y = y,
  N = N,
  design = design,
  nu = maternconstants$nu,
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

quad <- aghq::marginal_laplace_tmb(obj, 1, startingvalue = c(param$logkappa, param$logtau))

# Sampling from AGHQ fitted model

cmr <- sf::st_read("figures/naomi-aghq/gadm41_CMR_2.json")
nga <- sf::st_read("figures/naomi-aghq/gadm41_NGA_2.json")
areas <- rbind(cmr, nga)

loaloa_sf <- mutate(loaloa_sf, p = y / N)

grid <- raster::raster(loaloa_sf, nrows = 30, ncols = 30)
grid_stars <- stars::st_as_stars(grid)
b <- gstat::krige(p ~ 1, loaloa_sf, grid_stars)
b <- gstat::krige(p ~ 1, loaloa_sf, grid_stars, nsim = 1)
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

# loazippostsamples <- sample_marginal(quad, 100)
# postU <- loazippostsamples$samps[c(1:190), ]
# postV <- loazippostsamples$samps[c(192:381), ]
# postBeta <- loazippostsamples$samps[c(191, 382), ]
# 
# U <- postU  
# theta <- loazippostsamples$theta
# pointsdata <- loaloa_sf
# resolution <- reslist
# 
# modpar <- cbind(
#   var = get_sigma(exp(theta$logkappa), exp(theta$logtau))^2,
#   range = get_rho(exp(theta$logkappa), exp(theta$logtau)),
#   shape = maternconstants$nu
# )
# 
# fielddat <- pointsdata
# fielddat@data <- as.data.frame(U)
# 
# geostatsp::RFsimulate(
#   model = modpar,
#   data = fielddat,
#   x = raster::raster(fielddat, nrows = resolution$nrow, ncols = resolution$ncol)
# )
