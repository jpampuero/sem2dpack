DOC = users_guide_sem2dpack

.SUFFIXES:

all: $(DOC).pdf

%.pdf: %.ps
	ps2pdf $<

%.ps: %.dvi
	dvips $* -o 

$(DOC).dvi: *.tex *.bib
	latex $(DOC)

selfdoc.tex: ../../SRC/*.f90 
	cd .. ; csh ./selfdoc.csh ; cd USERS_GUIDE

*.bib:
	bibtex $(DOC) ; latex $(DOC)
