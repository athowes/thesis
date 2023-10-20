library(tidyverse)
library(patchwork)
library(TMB)

TMB::compile("figures/naomi-aghq/2d.cpp")
dyn.load(TMB::dynlib("figures/naomi-aghq/2d"))

obj <- TMB::MakeADFun(data = list(), parameters = list(theta1 = 0, theta2 = 0), DLL = "2d")

box_lower <- -5
box_upper <- 10
box_size <- box_upper - box_lower

grid <- expand.grid(
  theta1 = seq(box_lower, box_upper, length.out = box_size * 50),
  theta2 = seq(box_lower, box_upper, length.out = box_size * 50)
)

ground_truth <- cbind(grid, pdf = apply(grid, 1, function(x) exp(-1 * obj$fn(x))))

opt <- nlminb(
  start = obj$par,
  objective = obj$fn,
  gradient = obj$gr,
  control = list(iter.max = 1000, trace = 0)
)

sd_out <- TMB::sdreport(
  obj,
  par.fixed = opt$par,
  getJointPrecision = TRUE
)

mu <- opt$par
cov <- sd_out$cov.fixed

fig0 <- ggplot(ground_truth, aes(x = theta1, y = theta2, z = pdf)) +
  geom_contour(col = "lightgrey") +
  coord_fixed(xlim = c(box_lower, box_upper), ylim = c(box_lower, box_upper), ratio = 1) +
  labs(x = "", y = "") +
  theme_minimal() +
  guides(size = "none") +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

#' R-INLA grid method
Lambda <- diag(eigen(cov)$values)
E <- eigen(cov)$vectors
E %*% Lambda %*% t(E)

z_to_theta <- function(z) {
  as.vector(opt$par + E %*% sqrt(Lambda) %*% z)
}

theta_mode <- z_to_theta(c(0, 0))
theta_mode

test_statistic <- function(theta_proposal) {
  obj$fn(theta_proposal) - obj$fn(theta_mode)
}

#' @param j is the index of exploration direction
#' @param m is dim(theta)
#' @param delta_z is the step size
#' @param delta_pi is the acceptable drop-off
explore_direction <- function(j, m, delta_z, delta_pi) {
  z_mode <- rep(0, m)
  z_names <- paste0("z", 1:m)
  names(z_mode) <- z_names
  
  unit_vector <- rep(0, m)
  unit_vector[j] <- 1
  
  points <- rbind(z_mode)
  
  # Increasing
  i <- 0
  condition <- TRUE
  while(condition) {
    i <- i + 1
    proposal <- c(z_mode + i * delta_z %*% unit_vector)
    names(proposal) <- z_names
    statistic <- test_statistic(z_to_theta(proposal))
    condition <- (statistic < delta_pi)
    if(condition){
      points <- rbind(points, proposal)
    }
  }
  
  # Decreasing
  i <- 0
  condition <- TRUE
  while(condition) {
    i <- i + 1
    proposal <- c(z_mode - i * delta_z %*% unit_vector)
    names(proposal) <- z_names
    statistic <- test_statistic(z_to_theta(proposal))
    condition <- (statistic < delta_pi)
    if(condition){
      points <- rbind(points, proposal)
    }
  }
  as.data.frame(points)
}

d_z <- 0.75
d_pi <- 2

z_grid <- expand.grid(
  z1 = explore_direction(1, 2, delta_z = d_z, delta_pi = d_pi)$z1,
  z2 = explore_direction(2, 2, delta_z = d_z, delta_pi = d_pi)$z2
)

theta_grid_full <- t(apply(z_grid, 1, z_to_theta)) %>%
  as.data.frame() %>%
  rename(theta1 = V1, theta2 = V2)

theta_grid <- theta_grid_full %>%
  mutate(
    statistic = apply(theta_grid_full, 1, test_statistic),
    condition = statistic < d_pi
  ) %>%
  filter(condition == TRUE)

fig1 <- fig0 +
  geom_point(
    data = theta_grid,
    aes(x = theta1, y = theta2, size = 1),
    alpha = 0.8,
    col = "#009E73",
    inherit.aes = FALSE
  ) +
  scale_size_continuous(range = c(1, 2)) +
  labs(tag = "A") +
  theme(
    plot.caption = element_text(hjust = 0.5, vjust = 1)
  )

ccd <- rsm::ccd(basis = 2, n0 = 1)
ccd <- as.data.frame(ccd)

z_to_theta <- function(z1, z2) {
  as.vector(opt$par + E %*% sqrt(Lambda) %*% c(z1, z2))
}

ccd <- t(mapply(z_to_theta, z1 = ccd$x1, z2 = ccd$x2)) %>%
  as.data.frame()
  
fig2 <- fig0 +
  geom_point(
    data = ccd,
    aes(x = V1, y = V2, size = 1),
    alpha = 0.8,
    col = "#009E73",
    inherit.aes = FALSE
  ) +
  scale_size_continuous(range = c(1, 2)) +
  labs(tag = "B") +
  theme(
    plot.caption = element_text(hjust = 0.5, vjust = 1)
  )

fig1 + fig2

ggsave("figures/naomi-aghq/inla-grid-demo.png", h = 3.5, w = 6.25)
