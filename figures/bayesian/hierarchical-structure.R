setwd(here::here("figures/bayesian"))
tools::texi2dvi("hierarchical-structure.tex", pdf = TRUE, clean = TRUE)
setwd(".")
