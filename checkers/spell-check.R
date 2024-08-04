library(hunspell)

spell_check <- function(filename) {
  content <- tolower(readLines(filename, warn = FALSE))
  
  words <- unlist(strsplit(content, "\\W+"))
  words <- words[words != ""]

  misspelled_words <- hunspell(words, dict = 'en_GB')
  misspelled_words <- unlist(misspelled_words)
  correct_words <- c(
    "hiv", "unaids", "howes", "eaton", "flaxman", "athowes", "aghq", "agyw", "naomi", "plhiv", "bayesian", "laplace", "inla",
    "knitr", "bookdown", "mcmc", "github", "mathbf", "mathcal", "poisson", "bayes", "tmbstan", "stan", "phd", "texttt", "ldots",
    "tex", "biblatex", "documentclass", "adjustmtc", "mgcv", "waterloo", "stringer", "elgm", "hermite", "tmb", "malawi", "africa",
    "unnormalised", "besag", "daly", "dalys", "markboth", "png", "hmc", "antiretroviral", "pepfar", "saharan", "boldsymbol", "glm",
    "glmm", "waic", "bic", "dic", "cdot", "propto", "hamiltonian", "mathbb", "scipen", "zambia", "african", "argmin", "zimbabwe",
    "risher", "mozambique", "lesotho", "kronecker", "icar", "bym", "fck", "fik", "ck", "ik", "crps", "mse"
  )
  misspelled_words <- setdiff(misspelled_words, correct_words)

  return(sort(misspelled_words))
}

spell_check("front-and-back-matter/_welcome-ebook.Rmd")
spell_check("front-and-back-matter/_originality.Rmd")
spell_check("front-and-back-matter/_copyright.Rmd")
spell_check("front-and-back-matter/_acknowledgements.Rmd")
spell_check("front-and-back-matter/_abstract.Rmd")

spell_check("01-introduction.Rmd")
spell_check("02-hiv-aids.Rmd")
spell_check("03-bayesian.Rmd")
spell_check("04-beyond-borders.Rmd")
spell_check("05-multi-agyw.Rmd")
spell_check("06-naomi-aghq.Rmd")
spell_check("07-conclusions.Rmd")
spell_check("90-appendixA.Rmd")
spell_check("91-appendixB.Rmd")
spell_check("92-appendixC.Rmd")
spell_check("references.bib")

spell_check("README.Rmarkdown")

spell_check("corrections/corrections.Rmd")
spell_check("corrections/cover.Rmd")
