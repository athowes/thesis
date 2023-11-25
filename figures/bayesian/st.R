library(tidyverse)
library(sf)
library(patchwork)

sf <- sf::st_read("figures/bayesian/zaf_areas.geojson") %>%
  filter(area_level == 2)

ct <- data.frame(name = "Cape Town", x = 18.375, y = -33.95) %>%
  sf::st_as_sf(coords = c("x", "y"), crs = "WGS84")

sf_lightgrey <- "#E6E6E6"

i <- 9
indicator <- c(rep(0, i - 1), 1, rep(0, nrow(sf) - i))

space <- ggplot() +
  geom_sf(data = sf, aes(fill = as.factor(indicator))) +
  annotate("text", x = 14.5, y = -33.95, label = "Cape Town\n(point)", size = 3, col = "grey30") +
  geom_sf(data = ct, aes(geometry = geometry), size = 2, col = "#D55E00") +
  annotate("text", x = 20, y = -23, label = "ZF Mgcawu\nDM (area)", size = 3, col = "grey30") +
  labs(x = "", y = "", fill = "", tag = "A") +
  lims(x = c(12.5, NA)) +
  scale_fill_manual(values = c(sf_lightgrey, "#0072B2")) +
  theme_void() +
  theme(legend.position = "none")

time <- data.frame(
    x = c(5, 12),
    y = c(2, 1),
    alpha = c(0, 1)
  ) %>%
  ggplot(aes(x = x, y = as.factor(y), alpha = alpha)) +
  geom_point(col = "#D55E00", size = 2) +
  geom_segment(aes(x = 4, y = 2, xend = 7, yend = 2), col = "#0072B2", linewidth = 0.5, arrow = arrow(ends = "both", length = unit(0.1, "inches"))) +
  scale_x_continuous(breaks = 1:12, labels = month.abb, limits = c(1, 12)) +
  scale_y_discrete(labels = c("World AIDS\nDay (point)", "Q2 (period)")) +
  scale_alpha(range = c(0, 1)) +
  guides(alpha = "none") +
  labs(x = "", y = "", tag = "B") +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

space + time + plot_layout(widths = c(1.35, 1))

ggsave("figures/bayesian/st.png", h = 2.5, w = 6.25)
