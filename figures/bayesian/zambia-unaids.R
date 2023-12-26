library(ggplot2)
library(ggpubr)
library(patchwork)

set.seed(2)

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

sf <- sf::st_read("figures/bayesian/zmb_areas_admin3.geojson")

sf_constituency_full <- dplyr::filter(sf, area_level == 3)
sf_constituency <- sf::st_simplify(sf_constituency_full, dTolerance = 1000)

simulate_icar <- function(W, sd = 1) {
  n <- ncol(W)
  num <- rowSums(W)
  Q <- -W
  diag(Q) <- num
  Q_aux <- eigen(Q)$vectors[, order(eigen(Q)$values)]
  D_aux <- sort(eigen(Q)$values)
  rnd <- rnorm(n - 1, 0, sqrt(sd * (1/D_aux[-1])))
  rnd <- Q_aux %*% c(0, rnd)
  return(as.vector(rnd))
}

nb <- spdep::poly2nb(sf_constituency_full)

nb_sf <- spdep::nb2lines(nb, coords = sp::coordinates(as(sf_constituency, "Spatial"))) %>%
  as("sf") %>%
  sf::st_set_crs(sf::st_crs(sf_constituency))

W <- spdep::nb2mat(neighbours = nb, style = "B", zero.policy = TRUE)

beta <- -2
u <- simulate_icar(W, sd = 0.5)

eta <- beta + u
rho <- plogis(eta)
sf_constituency$rho <- rho

m1 <- 5
y1 <- rbinom(n = nrow(sf_constituency), size = m1, prob = rho)
sf_constituency$direct1 <- y1 / m1

m2 <- 25
y2 <- rbinom(n = nrow(sf_constituency), size = m2, prob = rho)
sf_constituency$direct2 <- y2 / m2

m3 <- 125
y3 <- rbinom(n = nrow(sf_constituency), size = m3, prob = rho)
sf_constituency$direct3 <- y3 / m3

fit_model <- function(dat) {
  tau_prior <- list(prec = list(prior = "logtnormal", param = c(0, 1/2.5^2), initial = 0, fixed = FALSE))
  beta_prior <- list(mean.intercept = -2, prec.intercept = 1)
  spdep::nb2INLA("sf.adj", nb)
  g <- INLA::inla.read.graph(filename = "sf.adj")
  formula <- y ~ 1 + f(id, model = "besag", graph = g, scale.model = TRUE, constr = TRUE, hyper = tau_prior)
  
  fit <- INLA::inla(
    formula, family = "binomial", control.family = list(control.link = list(model = "logit")),
    control.fixed = beta_prior, data = dat, Ntrials = m, control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
  )
  
  return(fit)
}

dat1 <- list(id = 1:nrow(sf_constituency), y = y1, m = m1)
fit1 <- fit_model(dat1)

dat2 <- list(id = 1:nrow(sf_constituency), y = y2, m = m2)
fit2 <- fit_model(dat2)

dat3 <- list(id = 1:nrow(sf_constituency), y = y3, m = m3)
fit3 <- fit_model(dat3)

sf_constituency$modelled1 <- fit1$summary.fitted.values$mean
sf_constituency$modelled2 <- fit2$summary.fitted.values$mean
sf_constituency$modelled3 <- fit3$summary.fitted.values$mean

maps <- sf_constituency %>%
  tidyr::pivot_longer(cols = c(direct1, direct2, direct3, modelled1, modelled2, modelled3), names_to = "id", values_to = "est") %>%
  tidyr::separate(id, c("type", "sample_size"), sep = "(?<=[a-zA-Z_])\\s*(?=[0-9])") %>%
  mutate(
    sample_size = 5^as.numeric(sample_size),
    type = stringr::str_to_title(type)
  ) %>%
  ggplot(aes(fill = est)) +
  geom_sf(size = 0.1, color = "grey30") +
  facet_grid(type ~ sample_size) +
  scale_fill_viridis_c(
    option = "C", direction = 1, limits = c(0, 0.6),
    labels = scales::label_percent(1), na.value = viridis::viridis(20, option = "C")[1]
  ) +
  labs(fill = "") +
  theme_void() +
  theme(
    legend.position = "left",
    legend.key.width = unit(0.75, 'cm'),
    strip.text.y = element_text(angle = 270)
  )

maps

ggsave("figures/bayesian/zmb-maps.png", h = 3.5, w = 6.25)

scatter <- sf_constituency %>%
  sf::st_drop_geometry() %>%
  tidyr::pivot_longer(cols = c(direct1, direct2, direct3, modelled1, modelled2, modelled3), names_to = "id", values_to = "est") %>%
  tidyr::separate(id, c("type", "sample_size"), sep = "(?<=[a-zA-Z_])\\s*(?=[0-9])") %>%
  mutate(
    sample_size = 5^as.numeric(sample_size),
    type = stringr::str_to_title(type)
  ) %>%
  dplyr::select(rho, sample_size, type, est) %>%
  ggplot(aes(x = rho, y = est)) +
  geom_point(alpha = 0.5, shape = 1) +
  facet_grid(type ~ sample_size) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  stat_cor(aes(label = after_stat(r.label)), method = "pearson", label.x = 0.6, label.y = 0.1, p.digits = 3) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  coord_fixed(ratio = 1) +
  labs(x = "Underlying truth", y = "Estimate") +
  theme_minimal()

scatter

ggsave("figures/bayesian/zmb-scatter.png", h = 4.25, w = 6.25, bg = "white")
