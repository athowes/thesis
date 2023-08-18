setwd(here::here("figures/introduction/"))
tools::texi2dvi("chapter-flowchart.tex", pdf = TRUE, clean = TRUE)

convert_pdf_png <- function(name, dpi = 300) {
  command <- paste0(
    "convert -density ", dpi, " ", name, ".pdf -scene 1 -background white",
    " -alpha remove -alpha off -quality 90 ", name, ".png"
  )
  system(command)
}

convert_pdf_png(name = "chapter-flowchart", dpi = 300)

setwd(".")