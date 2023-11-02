library(TMB)
library(geostatsp)
data(loaloa, package = "geostatsp")
loaloa <- terra::vect(loaloa)
loaloa <- sf::st_as_sf(loaloa)
library(aghq)
library(tmbstan)

# Set the resolution for the spatial interpolations
reslist <- list(nrow = 50, ncol = 100)

compile("resources/naomi-aghq/loaloazip.cpp")
dyn.load(dynlib("resources/naomi-aghq/loaloazip"))

Amat <- Diagonal(nrow(loaloa))
Xmat <- cbind(rep(1, nrow(Amat)))
design <- bdiag(cbind(Amat, Xmat), cbind(Amat, Xmat))
y <- loaloa$y
N <- loaloa$N
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
  D = raster::pointDistance(loaloa, lonlat = FALSE),
  betaprec = beta_prec
)

startingsig <- 1
startingrho <- 4.22*1e04

set.seed(4564)

param <- list(
  W = rnorm(ncol(design)),
  logkappa = log(get_kappa(startingsig, startingrho)),
  logtau = log(get_tau(startingsig, startingrho))
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

simulate_spatial_fields <- function(U, theta, pointsdata, resolution = list(nrow = 100, ncol = 100)) {
  modpar <- cbind(
    var = get_sigma(exp(theta$theta1), exp(theta$theta2))^2,
    range = get_rho(exp(theta$theta1), exp(theta$theta2)),
    shape = maternconstants$nu
  )
  
  fielddat <- pointsdata
  fielddat@data <- as.data.frame(U)
  
  geostatsp::RFsimulate(
    model = modpar,
    data = fielddat,
    x = raster(fielddat, nrows = resolution$nrow,ncols = resolution$ncol)
  )
}

