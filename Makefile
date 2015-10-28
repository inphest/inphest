BIBFILES =
FIGS = $(wildcard src/figs/*)

.PHONY: all

all: index.html

index.html: src/inphest-development-journal.tex $(BIBFILES) $(FIGS)
	pandoc --mathjax='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' \
		   --standalone \
		   -s \
		   -o $@  \
		   src/inphest-development-journal.tex
