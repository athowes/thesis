library(readxl)
library(stringr)
library(tidyverse)
library(patchwork)

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

as_numeric_spaces <- function(x) {
  as.numeric(str_replace_all(x, " ", ""))
}

#' Obtained from https://aidsinfo.unaids.org/
infections_esa <- read_xlsx("resources/hiv-aids/epidemic-transition-metrics-esa.xlsx", sheet = 1, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Eastern and southern Africa", Indicator = "New HIV infections")

deaths_esa <- read_xlsx("resources/hiv-aids/epidemic-transition-metrics-esa.xlsx", sheet = 2, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Eastern and southern Africa", Indicator = "AIDS-related deaths")

infections_wca <- read_xlsx("resources/hiv-aids/epidemic-transition-metrics-wca.xlsx", sheet = 1, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Western and central Africa", Indicator = "New HIV infections")

deaths_wca <- read_xlsx("resources/hiv-aids/epidemic-transition-metrics-wca.xlsx", sheet = 2, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Western and central Africa", Indicator = "AIDS-related deaths")

infections_global <- read_xlsx("resources/hiv-aids/epidemic-transition-metrics-global.xlsx", sheet = 1, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Global", Indicator = "New HIV infections")

deaths_global <- read_xlsx("resources/hiv-aids/epidemic-transition-metrics-global.xlsx", sheet = 2, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Global", Indicator = "AIDS-related deaths")

fct_reorg <- function(fac, ...) {
  fct_recode(fct_relevel(fac, ...), ...)
}

infections_deaths <- bind_rows(
  infections_esa,
  deaths_esa,
  infections_wca,
  deaths_wca,
  infections_global,
  deaths_global
) %>%
  mutate(
    Region = fct_reorg(Region, "Global" = "Global", "Eastern and\nsouthern Africa" = "Eastern and southern Africa", "Western and\ncentral Africa" = "Western and central Africa"),
    Indicator = fct_relevel(Indicator, "New HIV infections", "AIDS-related deaths")
  )

ggplot(infections_deaths, aes(x = Year, y = `All ages estimate`, ymax = `Upper Estimate`, ymin = `Lower Estimate`, fill = Region, col = Region)) +
  geom_ribbon(alpha = 0.25, colour = NA) +
  geom_line() +
  facet_grid(Indicator ~ .) +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6), limits = c(0, NA)) + 
  scale_fill_manual(values = cbpalette) +
  scale_color_manual(values = cbpalette) +
  theme_minimal() +
  labs(y = "Estimate")

ggsave("figures/hiv-aids/overall-picture.png", h = 4.25, w = 6.25)

infections_global %>%
  filter(`All ages estimate` == max(`All ages estimate`))

round(100 * (max(infections_global$`All ages estimate`) - tail(infections_global$`All ages estimate`, 1)) / max(infections_global$`All ages estimate`))

deaths_global %>%
  filter(`All ages estimate` == max(`All ages estimate`))

round(100 * (max(deaths_global$`All ages estimate`) - tail(deaths_global$`All ages estimate`, 1)) / max(deaths_global$`All ages estimate`))


