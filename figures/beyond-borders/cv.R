library(dplyr)
library(ggplot2)
library(patchwork)

set.seed(2)

sf <- sf::st_read("figures/beyond-borders/zwe_areas.geojson")
sf <- filter(sf, area_level == 2) %>% mutate(y = 1)
  
training_sets_sloo <- arealutils::create_folds(sf, type = "SLOO")
training_sets_loo <- arealutils::create_folds(sf, type = "LOO")

i <- 32

df <- rbind(
  training_sets_loo[[i]]$data %>%
    mutate(
      type = "Leave-one-out (LOO)",
      left_out = as.numeric(is.na(y))
    ),
  training_sets_sloo[[i]]$data %>%
    mutate(
      type = "Spatial-leave-one-out (SLOO)",
      left_out = as.numeric(is.na(y))
    )
)

df$left_out[c(training_sets_loo[[i]]$predict_on, nrow(sf) + training_sets_sloo[[i]]$predict_on)] <- 2

ggplot(df, aes(fill = as.factor(left_out))) +
  geom_sf(aes(geometry = geometry)) +
  coord_sf() +
  facet_grid(~ type) +
  scale_fill_manual(
    values = c("#E6E6E6", "#56B4E9", "#009E73"),
    name = "",
    labels = c("Training", "Left out", "Left out\nand predicted on")
  ) +
  guides(fill = guide_legend(override.aes = list(colour = NA))) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )

ggsave("figures/beyond-borders/cv.png", h = 4, w = 6.25)
