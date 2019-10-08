PANDOC=pandoc
BASEDIR=$(CURDIR)
INPUTDIR=$(BASEDIR)
OUTPUTDIR=$(BASEDIR)/out
BIBFILE=$(INPUTDIR)/references.bib

pdf:
	pandoc "$(INPUTDIR)"/*.md \
	-o "$(OUTPUTDIR)/Report, ex. 1.7, Tom Westerhout.pdf" \
	-H "$(INPUTDIR)/preamble.tex" \
	--template="$(INPUTDIR)/latex.template" \
	--bibliography="$(BIBFILE)" 2>pandoc.log \
	--csl="$(INPUTDIR)/ref_format.csl" \
	--highlight-style pygments \
	--standalone \
	-N \
	--pdf-engine=xelatex \
	--verbose

tex:
	pandoc "$(INPUTDIR)"/*.md \
	-o "$(OUTPUTDIR)/Report, ex. 1.7, Tom Westerhout.pdf" \
	--template="$(INPUTDIR)/latex.template" \
	--bibliography="$(BIBFILE)" 2>pandoc.log \
	--csl="$(INPUTDIR)/ref_format.csl" \
	--highlight-style pygments \
	--standalone \
	-N \
	--pdf-engine=xelatex \
	--verbose

.PHONY: pdf tex
