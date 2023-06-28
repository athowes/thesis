library(readr)
library(tidyverse)
library(patchwork)

incidence <- read_csv("resources/introduction/incidence-of-hiv-by-age.csv")

incidence_plot <- incidence %>%
  filter(Entity == "Sub-Saharan Africa (UN)") %>%
  select(
    year = Year,
    incidence = `3.3.1 - Number of new HIV infections per 1,000 uninfected population, by sex and age (per 1,000 uninfected population) - SH_HIV_INCD - All age ranges or no breaks by age - Both sexes`
  ) %>%
  filter(year >= 2000) %>%
  ggplot(aes(x = year, y = incidence)) +
    geom_point() +
    labs(x = "", y = "HIV incidence") +
    theme_minimal()

deaths <- read_csv("resources/introduction/hivaids-and-tuberculosis-deaths.csv")

deaths_plot <- deaths %>%
  filter(Entity == "Sub-Saharan Africa (WB)") %>%
  select(
    year = Year,
    deaths_direct = `Deaths - HIV/AIDS - Sex: Both - Age: All Ages (Number)`,
    deaths_indirect = `Deaths - HIV/AIDS resulting in other diseases - Sex: Both - Age: All Ages (Number)`
  ) %>%
  mutate(deaths = deaths_direct + deaths_indirect) %>%
  filter(year >= 2000) %>%
  ggplot(aes(x = year, y = deaths)) +
    geom_point() +
    lims(x = c(2000, 2020)) +
    labs(x = "Year", y = "AIDS related deaths") +
    theme_minimal()

ggsave("figures/overall-picture.png", incidence_plot / deaths_plot, h = 4, w = 6.25)
