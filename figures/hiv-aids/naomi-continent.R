library(tidyverse)
library(rnaturalearth)
library(sf)
library(rmapshaper)
library(scales)
library(patchwork)

sf::sf_use_s2(FALSE)

#' Obtained from HIV Inference Group Sharepoint
files <- list.files("resources/hiv-aids/2023 naomi final/", "zip$", full.names = TRUE)
out <- lapply(files, naomi::read_output_package)

get_iso3 <- function(x) {
  x$meta_area %>%
    filter(area_level == 0) %>%
    pull(area_id)
}

names(out) <- sapply(out, get_iso3)

redactions <- c("HTI", "NGA")
out[redactions] <- NULL

areas <- out %>%
  lapply("[[", "meta_area") %>%
  Map(mutate, ., iso3 = names(.)) %>%
  bind_rows() %>%
  select(iso3, everything())

areas_district <- areas %>%
  group_by(iso3) %>%
  filter(area_level == max(area_level))

data <- out %>%
  lapply("[[", "indicators") %>%
  Map(mutate, ., iso3 = names(.)) %>%
  bind_rows()

#' Estimates calendar_quarter only
data <- data %>%
  group_by(iso3) %>%
  filter(
    calendar_quarter == sort(unique(calendar_quarter))[2]
  ) %>%
  mutate(
    is_district = area_level == max(area_level),
    is_national = area_level == 0
  )

areas_district <- areas %>%
  semi_join(filter(data, is_district), by = "area_id")

areas_district_simple <- areas_district %>%
  ms_simplify(0.05, keep_shapes = TRUE)

nat_boundaries <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf") %>%
  mutate(
    longitude = st_coordinates(st_centroid(geometry))[,1],
    latitude = st_coordinates(st_centroid(geometry))[,2]
  )

#' Remove small island countries
nat_boundaries <- nat_boundaries %>%
  filter(!iso_a3 %in% c("CPV", "COM", "STP"))

#' Trim Prince Edward Island
bbox <- c(xmin = -17.5327797,
          ymin = -35, # -46.9697266,
          xmax = 51.4113159179688, 
          ymax = 37.3404121398926)

nat_boundaries <- st_crop(nat_boundaries, bbox)

iso3_sort <- nat_boundaries %>% arrange(-latitude) %>% pull(iso_a3)

lakes <- ne_download("medium", "lakes", category = "physical", returnclass = "sf")

ssa_lakes <- st_filter(lakes, nat_boundaries)

#' Retract Cabo Delgado Province
cabo_delago_areas <- areas %>%
  filter(parent_area_id == "MOZ_2_10") %>%
  pull(area_name)

data <- data %>%
  filter(!(area_name %in% cabo_delago_areas)) 

national_map_data <- data %>%
  filter(
    is_national,
    sex == "both",
    age_group == "Y015_049",
    indicator == "prevalence"
  ) %>%
  arrange(desc(mean)) %>%
  ungroup() %>%
  mutate(
    iso3_sort_order = row_number(),
    area_name = forcats::fct_recode(area_name,
                                    "Dem. Rep. Congo" = "Democratic Republic of The Congo",
                                    "Tanzania" = "United Republic of Tanzania",
                                    "Sao Tome and Pri." = "SAO TOME AND PRINCIPE",
                                    "Cen. Afr. Rep." = "Central African Republic",
                                    "Equatorial Guinea" = "Equatorial guinea")
  )

dotplot <- data %>%
  filter(
    is_district,
    sex == "both",
    age_group == "Y015_049",
    indicator == "prevalence"
  ) %>%
  left_join(
    select(national_map_data, iso3, area_name_national = area_name, iso3_sort_order),
    by = "iso3"
  ) %>%
  ggplot(aes(x = fct_rev(reorder(area_name_national, iso3_sort_order)), y = mean, col = mean)) +
    geom_jitter(width = 0, size = 1, alpha = 0.8) +
    viridis::scale_color_viridis(option = "A", direction = -1, begin = 0.3, end = 1.0) +
    geom_point(data = national_map_data,
      aes(x = fct_rev(reorder(area_name, iso3_sort_order)), y = mean),
          shape = 21, size = 2, fill = "white", col = "black", alpha = 0.9
    ) +
    scale_y_continuous(labels = function(x) paste0(100 * x, "%")) +
    coord_flip() +
    labs(y = "HIV prevalence (15-49)", x = "Country") +
    guides(col = FALSE) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6, hjust = 0), axis.title.x = element_text(vjust = -0.5))

map <- data %>%
  filter(
    is_district,
    sex == "both",
    age_group == "Y015_049",
    indicator == "prevalence"
  ) %>%
  left_join(areas_district_simple) %>%
  st_as_sf() %>%
  ggplot(aes(fill = mean)) +
  geom_sf(data = nat_boundaries, inherit.aes = FALSE, fill = "grey80", color = NA) +
  geom_sf(data = ssa_lakes, fill = "skyblue2", color = NA) +
  geom_sf(color = NA) +
  geom_sf(data = nat_boundaries, fill = NA, inherit.aes = FALSE, size = 0.3, color = "grey30") +
  expand_limits(fill = 0) +
  scale_fill_viridis_c(element_blank(), labels = label_percent(), option = "A", direction = -1, begin = 0.3, end = 1.0) +
  coord_sf(expand = FALSE) +
  theme_void(12) +
  labs(caption = "Source: UNAIDS Naomi model estimates, 2023") +
  theme(legend.position = c(0.1, 0.1),
        legend.just = c(0, 0),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(1, "lines"))

dotplot + map

ggsave("figures/hiv-aids/naomi-continent.png", h = 4.5, w = 6.25)

#' Numbers for text
data %>%
  filter(
    is_district,
    sex == "both",
    age_group == "Y015_049",
    indicator == "prevalence",
    substr(area_id, 1, 3) == "ZAF"
  ) %>%
  filter(
    mean == max(mean) | mean == min(mean)
  ) %>%
  mutate(mean_pct = round(100 * mean)) %>%
  pull(mean_pct, area_name)
