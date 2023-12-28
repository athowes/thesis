library(pdftools)
library(magick)
library(ggplot2)
library(patchwork)

pdf_file <- "figures/conclusions/global-update.pdf"

title_page <- pdftools::pdf_convert(pdf_file, pages = 1, filenames = "figures/conclusions/title.png", dpi = 450)
figure_page <- pdftools::pdf_convert(pdf_file, pages = 53, filenames = "figures/conclusions/figure.png", dpi = 450)

title_ggplot <- ggplot() +
  annotation_custom(grid::rasterGrob(image_read(title_page)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  labs(tag = "A")

figure_ggplot <- ggplot() +
  annotation_custom(grid::rasterGrob(image_read(figure_page)), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  labs(tag = "B")

title_ggplot + figure_ggplot + plot_layout(ncol = 2)

ggsave("figures/conclusions/global-update-conclusion.png", height = 4.5, width = 6.25, dpi = 450)

file.remove("figures/conclusions/title.png")
file.remove("figures/conclusions/figure.png")
