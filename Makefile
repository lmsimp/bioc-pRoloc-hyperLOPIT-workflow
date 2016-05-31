.PHONY: clean

%.md: %.Rmd
	Rscript -e 'require("knitr"); knit("$^")'

clean:
	rm -f *~
	rm -rf .Rcache
