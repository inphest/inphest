BIBFILES =
FIGS = $(wildcard src/figs/*)

.PHONY: all

all: index.html

index.html: src/inphest-development-journal.tex $(BIBFILES) $(FIGS)
	pandoc -s -o $@ src/inphest-development-journal.tex
