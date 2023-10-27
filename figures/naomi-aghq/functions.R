#' Return a quadrature rule at which to evaluate the Laplace marginal
#'
#' The grid is approximate in that it uses the Gaussian approximation to the mode
#' and standard deviation rather than the "true" mode and standard deviation which
#' one might calculate directly using the Laplace marginal. Most of the time these
#' should be very similar (Gaussian approximations tend to be good at capturing the
#' first two moments) so it should not matter.
#'
#' @param modeandhessian The row of `modesandhessians` containing the node which
#' is the mode of the Laplace approximation, or alternatively just the node which
#' has the highest log posterior evaluation.
#' @param i The index of the latent field to choose
#' @param k The number of AGHQ grid points to choose
create_approx_grid <- function(modeandhessian, i, k = 5) {
  mode <- modeandhessian[["mode"]][[1]]
  mode_i <- mode[i]
  H <- modeandhessian[["H"]][[1]]
  var_i <- diag(solve(H))[i]
  # LL <- Cholesky(H, LDL = FALSE)
  # var_i <- (colSums(solve(expand(LL)$L)^2))[i]
  gg <- mvQuad::createNIGrid(dim = 1, type = "GHe", level = k) # Create Gauss-Hermite quadrature
  mvQuad::rescale(gg, m = mode_i, C = var_i) # Adapt to mode_i and sd_i
  return(gg)
}

#' Calculate weighted sum of probabilities using `matrixStats::logSumExp`
logSumExpWeights <- function(lp, w) {
  matrixStats::logSumExp(log(w) + lp)
}

#' Given a small number of log function evaluations `lps` at points `nodes`
#' calculate PDF and CDF using spline or polynomial interpolation
#'
#' @param nodes Locations at which the function has been evaluated
#' @param lps Log-posterior function evaluations
#' @param method Not in use currently
#' @param normalise Use the trapezoid rule to normalise the posterior at the finegrid stage?
compute_pdf_and_cdf <- function(nodes, lps, method = "auto", normalise = FALSE) {
  k <- length(nodes)
  if(k >= 4) method <- "spline"
  if(k < 4) method <- "polynomial"
  
  rn <- range(nodes)
  rnl <- diff(rn)
  min <- min(rn) - rnl / 2
  max <- max(rn) + rnl / 2
  
  if(method == "spline") {
    ss <- splines::interpSpline(nodes, lps, bSpline = TRUE, sparse = TRUE)
    interpolant <- function(x) { as.numeric(stats::predict(ss, x)$y) }
  }
  
  if(method == "polynomial") {
    interpolant <- as.function(polynom::poly.calc(x = nodes, y = lps))
  }
  
  finegrid <- seq(min, max, length.out = 1000)
  lps <- interpolant(finegrid)
  
  if(normalise) {
    #' Make sure that the log probabilities produce a normalised PDF
    logC <- trapezoid_rule_log(lps, spacing = finegrid[2] - finegrid[1])
    lps <- lps - logC
  }
  
  df <- data.frame(
    x = finegrid,
    pdf = exp(lps),
    cdf = cumsum(exp(lps)) * c(0, diff(finegrid))
  )
  
  return(df)
}

#' Trapezoid integration rule on the log scale
#'
#' @param x Log value of function evaluated on a regular grid
#' @param spacing The distance between grid points
#' @return Integral of the function
trapezoid_rule_log <- function(x, spacing) {
  w <- rep(spacing, length(x))
  w[1] <- w[1] / 2
  w[length(x)] <- w[length(x)] / 2
  matrixStats::logSumExp(log(w) + x)
}