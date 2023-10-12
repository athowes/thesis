devtools::install_github("metacran/cranlogs")
cran_downloads(when = "last-week", packages = c("ggplot2", "httr"))

cbpalette <- c("#56B4E9","#009E73", "#E69F00", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

packages <- cranlogs::cran_downloads(from = "2016-01-01", to = "2023-09-30", c("TMB", "glmmTMB", "rstan", "nimble", "brms"))

packages %>%
  mutate(
    my = lubridate::floor_date(date, unit = "month")
  ) %>%
  group_by(my, package) %>%
  summarise(
    count = sum(count)
  ) %>%
  ggplot(aes(x = my, y = count, col = package)) +
    geom_line() +
    scale_y_continuous(labels = scales::unit_format(unit = "K", scale = 1e-3)) +
    scale_colour_manual(values = cbpalette) +
    labs(x = "", y = "Monthly CRAN downloads", col = "R package") +
    theme_minimal()

ggsave("figures/naomi-aghq/package-downloads.png", h = 3, w = 6.25, bg = "white")
