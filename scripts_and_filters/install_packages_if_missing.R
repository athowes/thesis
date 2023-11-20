for (package in c("rmarkdown", "bookdown", "knitr", "kableExtra", "tidyverse", "here")) {
  print(paste0("checking for install of ", package))
  if (!requireNamespace(package)) install.packages(package, repos = "http://cran.rstudio.com")
  library(package, character.only = TRUE)
}
