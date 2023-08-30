find_adjacent_duplicates <- function(filename) {
  content <- tolower(readLines(filename, warn = FALSE))
  
  words <- unlist(strsplit(content, "\\W+"))
  words <- words[words != ""]
  
  adjacent_duplicates <- character(0)
  
  for (i in 2:length(words)) {
    if (words[i] == words[i - 1]) {
      adjacent_duplicates <- c(adjacent_duplicates, words[i])
    }
  }
  
  return(adjacent_duplicates)
}

find_adjacent_duplicates("01-introduction.Rmd")
find_adjacent_duplicates("02-hiv-aids.Rmd")
find_adjacent_duplicates("03-bayesian.Rmd")
find_adjacent_duplicates("04-beyond-borders.Rmd")
find_adjacent_duplicates("05-multi-agyw.Rmd")
find_adjacent_duplicates("06-naomi-aghq.Rmd")
find_adjacent_duplicates("07-conclusions.Rmd")
find_adjacent_duplicates("90-appendixA.Rmd")
find_adjacent_duplicates("91-appendixB.Rmd")
find_adjacent_duplicates("92-appendixC.Rmd")