%.tex: %.Rnw
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); knit("$^")'

%.md: %.tex
	pandoc $^ -o $@ 

clean:
	rm -f *~
	rm -rf .Rcache

.PHONY: clean all
