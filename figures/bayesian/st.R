library(tidyverse)
library(sf)
library(patchwork)

sf <- sf::st_read("figures/bayesian/zaf_areas.geojson") %>%
  filter(area_level == 2)

ct <- data.frame(name = "Cape Town", x = 18.375, y = -33.95) %>%
  sf::st_as_sf(coords = c("x", "y"), crs = "WGS84")

sf_lightgrey <- "#E6E6E6"
darkgrey <- "#4A4A4A"

i <- 9
indicator <- c(rep(0, i - 1), 1, rep(0, nrow(sf) - i))

space <- ggplot() +
  geom_sf(data = sf, aes(fill = as.factor(indicator))) +
  annotate("text", x = 15, y = -33.95, label = "Cape Town", size = 3) +
  geom_sf(data = ct, aes(geometry = geometry), size = 2, col = darkgrey) +
  annotate("text", x = 20, y = -23.5, label = "ZF Mgcawu\nDistrict Municipality", size = 3) +
  labs(x = "", y = "", fill = "") +
  lims(x = c(13, NA)) +
  scale_fill_manual(values = c(sf_lightgrey, darkgrey)) +
  theme_minimal() +
  theme(legend.position = "none")

time <- data.frame(
    x = c(5, 12),
    y = c(2, 1),
    alpha = c(0, 1)
  ) %>%
  ggplot(aes(x = x, y = as.factor(y), alpha = alpha)) +
  geom_point(col = darkgrey, size = 2) +
  geom_segment(aes(x = 4, y = 2, xend = 7, yend = 2), col = darkgrey, size = 1, arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  scale_x_continuous(breaks = 1:12, labels = month.abb, limits = c(1, 12)) +
  scale_y_discrete(labels = c("World AIDS Day", "Q2")) +
  coord_flip() +
  scale_alpha(range = c(0, 1)) +
  guides(alpha = "none") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

space + time + plot_layout(widths = c(1.25, 1))

ggsave("figures/bayesian/st.png", h = 2.5, w = 6.25)
