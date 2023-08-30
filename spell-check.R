library(hunspell)

spell_check <- function(filename) {
  content <- tolower(readLines(filename, warn = FALSE))
  
  words <- unlist(strsplit(content, "\\W+"))
  words <- words[words != ""]

  misspelled_words <- hunspell(words, dict = 'en_GB')
  misspelled_words <- unlist(misspelled_words)
  common_words <- c("hiv", "unaids", "howes", "eaton", "flaxman", "athowes", "aghq", "agyw", "naomi", "plhiv", "bayesian", "laplace", "inla")
  misspelled_words <- setdiff(misspelled_words, common_words)

  return(misspelled_words)
}

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