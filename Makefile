pdf:
	Rscript -e 'options(bookdown.render.file_scope = FALSE); bookdown::render_book("index.Rmd", output_format = "bookdown::pdf_book", output_dir = "docs")'
	Rscript -e 'file.remove(list.files(pattern = "*.(log|mtc|maf|aux|bcf|lof|lot|out|toc)"), here::here("front-and-back-matter", "abbreviations.aux"))'
	Rscript -e 'file.rename("docs/_main.pdf", "docs/main.pdf")'
	Rscript -e 'file.rename("docs/_main.tex", "docs/main.tex")'
	Rscript -e 'browseURL(here::here("docs", "main.pdf"))'

bs4book:
	Rscript -e 'options(bookdown.render.file_scope = FALSE); bookdown::render_book("index.Rmd", output_format = "bookdown::bs4_book", output_dir = "docs")'
	Rscript -e 'file.create(here::here("docs", ".nojekyll"))'
	Rscript -e 'browseURL(here::here("docs", "index.html"))'

gitbook:
	Rscript -e 'options(bookdown.render.file_scope = FALSE); bookdown::render_book("index.Rmd", output_format = "bookdown::gitbook", output_dir = "docs")'
	Rscript -e 'file.create(here::here("docs", ".nojekyll"))'
	Rscript -e 'browseURL(here::here("docs", "index.html"))'

word:
	Rscript -e 'options(bookdown.render.file_scope = FALSE); bookdown::render_book("index.Rmd", output_format = "bookdown::word_document2", output_dir = "docs")'
	Rscript -e 'browseURL(here::here("docs", "_main.docx"))'

clean:
	Rscript -e 'file.remove(list.files(pattern = "*.(log|mtc|maf|aux|bbl|blg|xml)"))'
	
clean-knits:
	Rscript -e 'file.remove(list.files(pattern = "*.(docx|html|pdf|log|maf|mtc|tex|toc|out|lof|lot|bcf|aux)"))'
	Rscript -e 'unlink(list.files(pattern = "*_(files|cache)"), recursive = TRUE)'
