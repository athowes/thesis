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

gg <- mvQuad::createNIGrid(2, "GHe", 3)

add_points <- function(fig0, gg) {
  
  points <- mvQuad::getNodes(gg) %>%
    as.data.frame() %>%
    mutate(weights = mvQuad::getWeights(gg))
  
  colnames(points) <- c("theta1", "theta2", "weights")
  
  fig0 +
    geom_point(
      data = points,
      aes(x = theta1, y = theta2, size = weights),
      alpha = 0.8,
      col = "#009E73",
      inherit.aes = FALSE
    ) +
    scale_size_continuous(range = c(1, 2)) +
    theme(
      plot.caption = element_text(hjust = 0.5, vjust = 1)
    )
}

fig1 <- add_points(fig0, gg) +
  labs(tag = "A", size = "", caption = "GHQ")

#' Adapt by the mean
gg2 <- gg
mvQuad::rescale(gg2, m = mu, C = diag(c(1, 1)), dec.type = 1)

fig2 <- add_points(fig0, gg2) +
  labs(tag = "B", size = "", caption = "Shifted")

#' Adapt by the lower Cholesky
gg3 <- gg
mvQuad::rescale(gg3, m = mu, C = cov, dec.type = 2)

fig3 <- add_points(fig0, gg3) +
  labs(tag = "C", size = "", caption = "AGHQ (Cholesky)")

#' Adapt by the spectral
gg4 <- gg
mvQuad::rescale(gg4, m = mu, C = cov, dec.type = 1)

fig4 <- add_points(fig0, gg4) +
  labs(tag = "D", size = "", caption = "AGHQ (spectral)")

fig1 + fig2 + fig3 + fig4

ggsave("figures/naomi-aghq/aghq-demo.png", h = 7, w = 6.25)

#' PCA-AGHQ
gg5 <- mvQuad::createNIGrid(2, "GHe", level = c(3, 1))
mvQuad::rescale(gg5, m = mu, C = cov, dec.type = 1)

lambda <- eigen(cov)$values
cumsum(lambda) / sum(lambda)

xstart <- 6.2
ystart <- -2.3

x1end <- xstart + 4 * eigen(cov)$vectors[1, 1]
y1end <- ystart + 4 * eigen(cov)$vectors[2, 1]

x2end <- xstart + 1 * eigen(cov)$vectors[1, 2]
y2end <- ystart + 1 * eigen(cov)$vectors[2, 2]

fig5 <- add_points(fig0, gg5) +
  geom_segment(aes(x = xstart, y = ystart, xend = x1end, yend = y1end), arrow = arrow(length = unit(0.25, "cm")), col = "darkgrey") +
  annotate("text", x = x1end + 1, y = y1end - 3, label = "95%", col = "darkgrey") +
  geom_segment(aes(x = xstart, y = ystart, xend = x2end, yend = y2end), arrow = arrow(length = unit(0.25, "cm")), col = "darkgrey") +
  annotate("text", x = x2end, y = y2end - 2, label = "5%", col = "darkgrey") +
  labs(tag = "B", size = "", caption = "PCA-AGHQ")

fig4 <- fig4 +
  labs(tag = "A", size = "", caption = "AGHQ (spectral)")

fig4 + fig5

ggsave("figures/naomi-aghq/pca-demo.png", h = 4, w = 6.25)

# Doing the quadrature
source("figures/naomi-aghq/functions.R")

run_quad <- function(x) {
  logSumExpWeights(lp = apply(mvQuad::getNodes(x), 1, function(y) - obj$fn(y)), w = mvQuad::getWeights(x))
}

quad <- run_quad(gg)
quad2 <- run_quad(gg2)
quad3 <- run_quad(gg3)
quad4 <- run_quad(gg4)

c(quad, quad2, quad3, quad4)
