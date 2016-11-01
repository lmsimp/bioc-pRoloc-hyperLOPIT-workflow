setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif


all:
	make bioc-workflow.tex
	make bioc-workflow.md
	make bioc-workflow.R

%.tex: %.Rnw
	"$(R_HOME)/bin/Rscript" -e 'require("knitr"); knit("$^")'

%.md: %.tex
	pandoc $^ -o $@ 

%.R: %.Rnw
	"$(R_HOME)/bin/Rscript" -e 'require("knitr"); purl("$^")'

clean:
	rm -f *~
	rm -rf .Rcache

.PHONY: clean all
