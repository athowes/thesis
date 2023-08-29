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

plot <- data.frame(x = finegrid, y = sin(finegrid)) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(col = "#009E73") +
  theme_minimal() +
  labs(x = "", y = "")

trapezoid_plot <- function(N) {
  grid <- seq(0, pi, length.out = N)
  int <- trapezoid_rule(x = sin(grid), spacing = grid[2] - grid[1])
  
  plot + 
    geom_bar(
      data = data.frame(x = grid, y = sin(grid)),
      aes(x = x, y = y), alpha = 0.7, stat = "identity",
      inherit.aes = FALSE, fill = "#56B4E9") +
    theme_minimal() +
    labs(
      subtitle = paste0("Number of nodes: ", N - 2, "\nTrapezoid rule estimate: ", round(int, 3)),
      x = "", y = ""
    )
}

fig1 <- trapezoid_plot(N = 7)
fig2 <- trapezoid_plot(N = 27)
fig3 <- trapezoid_plot(N = 127)

fig1 + fig2 + fig3 + plot_layout(ncol = 1)

ggsave("trapezoid.png", h = 7, w = 6.25)
