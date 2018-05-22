all:
	make main.pdf

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

## R_HOME=/opt/Rpatched/

%.tex: %.Rnw
	"$(R_HOME)/bin/Rscript" -e 'require("knitr"); knit("$^")'

%.md: %.tex
	pandoc $^ -o $@

%.R: %.Rnw
	"$(R_HOME)/bin/Rscript" -e 'require("knitr"); purl("$^")'

main.pdf: bioc-workflow-resubmit.tex
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex

main-annotated.pdf: bioc-workflow-annotated.tex
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex
	cp main.pdf main-annotated.pdf

clean:
	rm -f *~
	rm -rf .Rcache

replies-to-reviewers.pdf: replies-to-reviewers.Rmd
	"$(R_HOME)/bin/Rscript" -e 'library(rmarkdown); render("$^", output_format = pdf_document())'

.PHONY: clean all bibtex
