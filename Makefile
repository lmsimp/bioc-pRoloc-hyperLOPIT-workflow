all:
	make bioc-workflow.tex
	make bioc-workflow.md
	make bioc-workflow.R

%.tex: %.Rnw
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); knit("$^")'

%.md: %.tex
	pandoc $^ -o $@ 

%.R: %.Rnw
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); purl("$^")'

clean:
	rm -f *~
	rm -rf .Rcache

.PHONY: clean all
