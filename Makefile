default: html-docs

html-docs: index.html

%.html: %.Rmd
	Rscript --vanilla -e "rmarkdown::render(\"$*.Rmd\",output_format=\"html_document\")"

%.html: %.md
	pandoc $*.md > $*.html

%.pdf: %.tex
	pdflatex $*.tex

%.tex: %.Rnw
	Rscript -e "library(knitr); knit(\"$*.Rnw\")"

%.R: %.Rmd
	Rscript --vanilla -e "library(knitr); purl(\"$*.Rmd\",output=\"$*.R\")"

clean:
	$(RM) *.o *.so *.log *.aux *.out *.nav *.snm *.toc *.bak
	$(RM) Rplots.ps Rplots.pdf

fresh: clean
	$(RM) -r cache figure

