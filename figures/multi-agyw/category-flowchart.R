setwd(here::here("figures/multi-agyw"))
tools::texi2dvi("category-flowchart.tex", pdf = TRUE, clean = TRUE)

convert_pdf_png <- function(name, dpi = 300) {
  command <- paste0(
    "convert -density ", dpi, " ", name, ".pdf -scene 1 -background white",
    " -alpha remove -alpha off -quality 90 ", name, ".png"
  )
  system(command)
}

convert_pdf_png(name = "category-flowchart", dpi = 300)

setwd(".")
