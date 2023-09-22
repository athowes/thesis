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

fig0 + labs(tag = "A") + fig0 + labs(tag = "B")

ggsave("figures/naomi-aghq/inla-grid-demo.png", h = 3.5, w = 6.25)
