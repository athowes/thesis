library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(2)

sf <- sf::st_read("figures/beyond-borders/mwi_areas.geojson")
sf <- filter(sf, area_level == 3)

nrow(sf)

plot_samples <- function(samples, tag){
  ggplot(sf) +
    geom_sf(fill = "#E6E6E6") +
    geom_sf(data = samples, size = 0.75, col = "#009E73") +
    labs(x = "", y = "") +
    theme_minimal() +
    labs(tag = tag, fill = "") +
    theme_void()
}

L <- 10
n <- nrow(sf)

centroid <- sf::st_centroid(sf)
random <- sf::st_sample(sf, size = rep(L, n))
hexagonal <- sf::st_sample(sf::st_transform(sf, 3857), size = rep(L, n), type = "hexagonal", exact = TRUE)
regular <- sf::st_sample(sf, size = rep(L, n), type = "regular", exact = TRUE)

plot_samples(centroid, tag = "A") + plot_samples(random, tag = "B") + plot_samples(hexagonal, tag = "C") + plot_samples(regular, tag = "D") +  plot_layout(ncol = 4)

ggsave("figures/beyond-borders/integration-strategy.png", h = 3.5, w = 6.25)
