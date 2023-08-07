library(readxl)
library(stringr)
library(tidyverse)
library(patchwork)

ihme <- readr::read_csv("resources/background/IHME-GBD_2019_DATA-a295712c-1.csv")

unique(ihme$age_name)

ihme %>%
  filter(age_name != "28-364 days") %>%
  group_by(cause_name) %>%
  summarise(val = sum(val)) %>%
  arrange(desc(val)) %>%
  head(n = 8) %>%
  mutate(
    cause_name = fct_recode(
      cause_name,
      "Lower respiratory\ninfections" = "Lower respiratory infections",
      "Ischemic heart\ndisease" = "Ischemic heart disease",
      "Diarrheal\ndiseases" = "Diarrheal diseases"      
    ),
    is_hiv = ifelse(cause_name == "HIV/AIDS", TRUE, FALSE)
  ) %>%
  ggplot(aes(x = reorder(cause_name, val), y = val, fill = is_hiv)) +
    geom_col() +
    scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6), limits = c(0, NA)) +
    scale_fill_manual(values = c("grey75", "#CC79A7")) +
    coord_flip() +
    guides(fill = "none") +
    labs(x = "Cause", y = "Disability-adjusted life years (ages >1)", caption = "Source: Global Burden of Disease Study, 2019") +
    theme_minimal()

ggsave("figures/gbd.png", h = 3, w = 6.25)
