setwd(here::here("figures"))
tools::texi2dvi("hierarchical-structure.tex", pdf = TRUE, clean = TRUE)
setwd(".")
