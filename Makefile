.PHONY: clean all

all:
	make bioc-workflow.md

%.md: %.Rmd
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); knit("$^")'

%.tex: %.md
	pandoc $^ -o $@ -N --bibliography=bibliography.bib

clean:
	rm -f *~
	rm -rf .Rcache

