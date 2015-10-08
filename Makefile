all: slides

clean:
	rm *.pdf

SRC := $(wildcard *.rmd)
REPORTPDFS = $(SRC:.rmd=.pdf)
PANDOC_OPTIONS = -t beamer --highlight-style=tango -V fontsize=8pt --toc 

slides: $(REPORTPDFS)

%.md: %.rmd
	Rscript -e "require(knitr); require(markdown); knit('$<', '$@');"

.SECONDARY:

%.pdf: %.md
	pandoc $(PANDOC_OPTIONS) $< -o $@ 

