.PHONY: clean all

all:
	make bioc-workflow.md

%.md: %.Rmd
	Rscript -e 'require("knitr"); knit("$^")'

clean:
	rm -f *~
	rm -rf .Rcache

