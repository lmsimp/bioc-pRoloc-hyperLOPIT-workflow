all:
	make -B bioc-workflow-resubmit.tex

R_HOME=/opt/Rpatched/

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

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

clean:
	rm -f *~
	rm -rf .Rcache

.PHONY: clean all bibtex
