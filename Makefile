all:
	make bioc-workflow.tex
	make bioc-workflow.md
	make bioc-workflow.R

%.tex: %.Rnw
	"/Library/Frameworks/R.framework/Resources/Rscript" -e 'require("knitr"); knit("$^")'

%.md: %.tex
	pandoc $^ -o $@ 

%.R: %.Rnw
	"/Library/Frameworks/R.framework/Resources/Rscript" -e 'require("knitr"); purl("$^")'

clean:
	rm -f *~
	rm -rf .Rcache

.PHONY: clean all
