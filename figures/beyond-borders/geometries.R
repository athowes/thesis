library(ggplot2)
library(patchwork)
library(tidyverse)
library(sf)

area1 <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
area2 <- area1
area2[, 1] <- area2[, 1] + 1
area3 <- area1
area3[, 1] <- area3[, 1] + 2
geometry_1 <- lapply(list(area1, area2, area3), function(area) list(area) %>% sf::st_polygon()) %>%
  sf::st_sfc() %>%
  sf::st_sf()

figA <- ggplot() +
  geom_sf(data = geometry_1, fill = "#E6E6E6") +
  labs(tag = "A") +
  theme_void()

p <- st_point(c(0, 0))
radius1 <- st_buffer(p, dist = 1)
radius2 <- st_buffer(p, dist = 2)
radius3 <- st_buffer(p, dist = 3)
area1 <- radius1
area2 <- st_difference(radius2, radius1)
area3 <- st_difference(radius3, radius2)
geometry_2 <- lapply(list(area1, area2, area3), function(area) list(area) %>% sf::st_multipolygon()) %>%
  sf::st_sfc() %>%
  sf::st_sf()

figB <- ggplot() +
  geom_sf(data = geometry_2, fill = "#E6E6E6") +
  labs(tag = "B") +
  theme_void()

scale <- 0.7
area1 <- matrix(c(0, -0.5, 0, 0.5, 1, 1, 1, -1, 0, -0.5), ncol = 2, byrow = TRUE) * scale
area2 <- matrix(c(1, -1, 1, 1, 2, 1.5, 2, -1.5, 1, -1), ncol = 2, byrow = TRUE) * scale
area3 <- matrix(c(2, -1.5, 2, 1.5, 3, 2, 3, -2, 2, -1.5), ncol = 2, byrow = TRUE) * scale
geometry_3 <- lapply(list(area1, area2, area3), function(area) list(area) %>% sf::st_polygon()) %>%
  sf::st_sfc() %>%
  sf::st_sf()

figC <- ggplot() +
  geom_sf(data = geometry_3, fill = "#E6E6E6") +
  labs(tag = "C") +
  theme_void()

area1 <- matrix(c(0, 0, 0.5, 0, 0.5, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE)
area2 <- matrix(c(0.5, 0, 1.5, 0, 1.5, 1, 0.5, 1, 0.5, 0), ncol = 2, byrow = TRUE)
area3 <- matrix(c(1.5, 0, 3, 0, 3, 1, 1.5, 1, 1.5, 0), ncol = 2, byrow = TRUE)
geometry_4 <- lapply(list(area1, area2, area3), function(area) list(area) %>% sf::st_polygon()) %>%
  sf::st_sfc() %>%
  sf::st_sf()

figD <- ggplot() +
  geom_sf(data = geometry_4, fill = "#E6E6E6") +
  labs(tag = "D") +
  theme_void()

st_area(geometry_1)
st_area(geometry_2)
st_area(geometry_3)
st_area(geometry_4)

create_sf_grid <- function(height, width){
  sfc <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(width, 0), c(width, height), c(0, 0)))))
  grid <- sf::st_make_grid(sfc, cellsize = 1, square = TRUE)
  return(grid)
}

geometry_5 <- sf::st_sf(create_sf_grid(height = 6, width = 6))

figE <- ggplot() +
  geom_sf(data = geometry_5, fill = "#E6E6E6") +
  labs(tag = "E") +
  theme_void()

geometry_6 <- sf::st_geometry(sf::st_read("figures/beyond-borders/texas/U_S__House_District.shp"))
geometry_6 <- sf::st_sf(geometry_6)
sf::st_crs(geometry_6) <- NA

figF <- ggplot() +
  geom_sf(data = geometry_6, fill = "#E6E6E6") +
  labs(tag = "F") +
  theme_void()

geometry_7 <- readRDS("figures/beyond-borders/gadm36_CIV_2_sf.rds")
geometry_7 <- sf::st_sf(geometry_7)
sf::st_crs(geometry_7) <- NA

figG <- ggplot() +
  geom_sf(data = geometry_7, fill = "#E6E6E6") +
  labs(tag = "G") +
  theme_void()

layout <- "
AABBCCDD
AABBCCDD
#EEFFGG#
#EEFFGG#
"

figA + figB + figC + figD + figE + figF + figG + plot_layout(design = layout, widths = 1)

ggsave("figures/beyond-borders/geometries.png", h = 3.5, w = 6.25)
