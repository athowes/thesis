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
  trapezoid <- c(0, 0, rep(1:(N - 1), each = 4), 0, 0)
  grid <- seq(0, pi, length.out = N)
  int <- trapezoid_rule(x = f(grid), spacing = grid[2] - grid[1])
  
  df <- data.frame(
    index = 1:(4 * N),
    trapezoid = trapezoid,
    x = rep(grid, each = 4)
  )

  y <- f(grid)
  df$y <- c(rbind(y, 0, 0, y))
  df$int <- int
  df$spacing <- grid[2] - grid[1]
  df$N <- (N - 1)
  return(df)
}

df <- bind_rows(trapezoid_df(N = 6), trapezoid_df(N = 11), trapezoid_df(N = 21))

ggplot(df) +
  geom_polygon(aes(x, y, group = trapezoid), fill = "grey90", col = "grey80", alpha = 0.8) +
  geom_function(fun = f, n = 500, col = "#009E73") +
  geom_text(
    data = df %>% select(N, int) %>% unique(),
    aes(label = paste0("Estimate: ", signif(int, 5)), x = 0.45, y = 1.6), size = 3
  ) +
  facet_grid(N ~ .) +
  labs(x = "", y = "") +
  theme_minimal() 

ggsave("figures/naomi-aghq/trapezoid.png", h = 7, w = 6.25, bg = "white")
