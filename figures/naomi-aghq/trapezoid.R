library(tidyverse)
library(patchwork)

trapezoid_rule <- function(x, spacing) {
  # Assumes nodes are evenly spaced
  w <- rep(spacing, length(x)) # Weights given by space between nodes
  w[1] <- w[1] / 2 # Apart from the first which is halved
  w[length(x)] <- w[length(x)] / 2 # And the last, also halved
  sum(w * x) # Compute the weighted sum
}

finegrid <- seq(0, pi, length.out = 1000)

f <- function(x) sin(x) * x

trapezoid_df <- function(N) {
  grid <- seq(0, pi, length.out = N)
  df <- data.frame(x = grid, y = f(grid), N = N - 2)
  int <- trapezoid_rule(x = df$y, spacing = df$x[2] - df$x[1])
  df$int <- int
  df$spacing <- df$x[2] - df$x[1]
  return(df)
}

df <- bind_rows(trapezoid_df(N = 7), trapezoid_df(N = 27), trapezoid_df(N = 127))

ggplot(df) +
  geom_col(aes(x = x, y = y, width = spacing * 0.95), fill = "#56B4E9", alpha = 0.8) +
  geom_function(fun = f, n = 500, col = "#009E73") +
  geom_text(
    data = df %>% select(N, int) %>% unique(),
    aes(label = paste0("Estimate: ", signif(int, 5)), x = 0.45, y = 1.6), size = 3
  ) +
  facet_grid(N ~ .) +
  labs(x = "", y = "") +
  theme_minimal() 

ggsave("figures/naomi-aghq/trapezoid.png", h = 7, w = 6.25, bg = "white")
