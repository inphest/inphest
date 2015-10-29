BIBFILES =
FIGS = $(wildcard src/figs/*)

SOURCE_DIR = src
BUILD_DIR = build
PRODUCT_DIR = product

.PHONY: all

all: inphest-development-journal.pdf index.html

inphest-development-journal.pdf: src/inphest-development-journal.tex $(BIBFILES) $(FIGS)
	cd $(SOURCE_DIR) \
		&& latexmk -gg -halt-on-error -interaction=nonstopmode -file-line-error -pdf -use-make -output-directory=../$(BUILD_DIR) $(notdir $<) \
		&& cd .. \
		&& cp $(BUILD_DIR)/$(notdir $@) $(notdir $@)

index.html: inphest-development-journal.pdf
	# pandoc --mathjax='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' \
	# 	   --standalone \
	# 	   -s \
	# 	   -o $@  \
	# 	   src/inphest-development-journal.tex
	pdf2htmlEX inphest-development-journal.pdf index.html

