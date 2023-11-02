library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(2)

sf <- sf::st_read("figures/beyond-borders/zwe_areas.geojson")
sf <- filter(sf, area_level == 2)

nb <- arealutils::sf_to_nb(sf)

nb_sf <- spdep::nb2lines(nb, coords = sp::coordinates(as(sf, "Spatial"))) %>%
  as("sf") %>%
  sf::st_set_crs(sf::st_crs(sf))

figA <- ggplot(sf) +
  geom_sf() +
  theme_minimal() +
  labs(tag = "A") +
  theme_void()

figB <- ggplot(sf) +
  geom_sf(data = nb_sf) +
  theme_minimal() +
  labs(tag = "B") +
  theme_void()

figA + figB

ggsave("figures/beyond-borders/geometry-graph.png", h = 3, w = 6.25)
