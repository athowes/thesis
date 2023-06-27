x <- y <- c(1:4)

expand.grid(x, y) %>%
  mutate(risk = ifelse(Var2 > 1, Var1 * Var2, 0)) %>%
  ggplot(aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = risk)) +
  scale_fill_gradientn(
    colours = c("white", "#E69F00"),
    breaks = c(2, 14),
    labels = c("Less risk", "More risk") 
  ) +
  scale_x_continuous(breaks = c(1:4), labels = c("Low: <0.3%", "Moderate: 0.3-1%", "High: 1-3%", "Very high: >3%")) +
  scale_y_continuous(breaks = c(1:4), labels = c("None", "Low", "High", "Very\n high")) +
  labs(x = "Population-level HIV incidence", y = "Individual-level risk behaviour", title = bquote("Important we consider" ~ bold("both") ~ "population setting and individual behaviour"), fill = "") +
  theme_minimal() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "grey92", linewidth = 0.5, arrow = arrow(length = unit(5, "pt"), type = "open"))
  )