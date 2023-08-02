library(readxl)
library(stringr)
library(tidyverse)
library(patchwork)

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

as_numeric_spaces <- function(x) {
  as.numeric(str_replace_all(x, " ", ""))
}

infections_esa <- read_xlsx("resources/introduction/epidemic-transition-metrics-esa.xlsx", sheet = 1, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Eastern and southern Africa", Indicator = "New HIV infections")

deaths_esa <- read_xlsx("resources/introduction/epidemic-transition-metrics-esa.xlsx", sheet = 2, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Eastern and southern Africa", Indicator = "AIDS-related deaths")

infections_wca <- read_xlsx("resources/introduction/epidemic-transition-metrics-wca.xlsx", sheet = 1, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Western and central Africa", Indicator = "New HIV infections")

deaths_wca <- read_xlsx("resources/introduction/epidemic-transition-metrics-wca.xlsx", sheet = 2, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Western and central Africa", Indicator = "AIDS-related deaths")

infections_global <- read_xlsx("resources/introduction/epidemic-transition-metrics-global.xlsx", sheet = 1, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Global", Indicator = "New HIV infections")

deaths_global <- read_xlsx("resources/introduction/epidemic-transition-metrics-global.xlsx", sheet = 2, skip = 1) %>%
  select(-5) %>%
  mutate(across(1:4, as_numeric_spaces)) %>%
  mutate(Region = "Global", Indicator = "AIDS-related deaths")

infections_deaths <- bind_rows(
  infections_esa,
  deaths_esa,
  infections_wca,
  deaths_wca,
  infections_global,
  deaths_global
) %>%
  mutate(Region = fct_relevel(Region, "Global", "Eastern and southern Africa", "Western and central Africa"))

infections_plot <- infections_deaths %>%
  filter(Indicator == "New HIV infections") %>%
  ggplot(aes(x = Year, y = `All ages estimate`, ymax = `Upper Estimate`, ymin = `Lower Estimate`, fill = Region, col = Region)) +
  geom_ribbon(alpha = 0.25, colour = NA) +
  geom_line() +
  labs(x = "", y = "New HIV infections") +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6), limits = c(0, NA)) + 
  scale_fill_manual(values = cbpalette) +
  scale_color_manual(values = cbpalette) +
  theme_minimal() +
  theme(legend.position = "none")

deaths_plot <- infections_deaths %>%
  filter(Indicator == "AIDS-related deaths") %>%
  ggplot(aes(x = Year, y = `All ages estimate`, ymax = `Upper Estimate`, ymin = `Lower Estimate`, fill = Region, col = Region)) +
  geom_ribbon(alpha = 0.25, colour = NA) +
  geom_line() +
  labs(x = "Year", y = "AIDS-related deaths") +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6), limits = c(0, NA)) + 
  scale_fill_manual(values = cbpalette) +
  scale_color_manual(values = cbpalette) +
  theme_minimal() +
  theme(legend.position = "bottom")


infections_global %>%
  filter(`All ages estimate` == max(`All ages estimate`))

round(100 * (max(infections_global$`All ages estimate`) - tail(infections_global$`All ages estimate`, 1)) / max(infections_global$`All ages estimate`))

deaths_global %>%
  filter(`All ages estimate` == max(`All ages estimate`))

round(100 * (max(deaths_global$`All ages estimate`) - tail(deaths_global$`All ages estimate`, 1)) / max(deaths_global$`All ages estimate`))

ggsave("figures/overall-picture.png", infections_plot / deaths_plot, h = 4, w = 6.25)
