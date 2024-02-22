rmarkdown::render(input = "slides/slides.Rmd", output_dir = "docs", output_format = "beamer_presentation")
rmarkdown::render(input = "slides/slides-long.Rmd", output_dir = "docs", output_format = "beamer_presentation")

png_file <- paste0("docs/slide-", 1, ".png")
pdftools::pdf_convert("docs/slides.pdf", pages = 1, format = "png", filenames = png_file, dpi = 300)
