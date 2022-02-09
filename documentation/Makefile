REV=documentation

# make or make dvi makes dvi file, including an update of the refs.
# make pdf creates pdf file
# make tar creates tgz-archive including all source-files, pdf-file and resources
# make install creates new dir and installs tgz file there
# make clean deletes all aux. files except tex-files, pdf file and tgz file

TEXFILE = $(REV).tex
AUXFILE = $(REV).aux

TEXFILES = definitions.tex
BIBFILES = bib_levin.bib

ARCHIVFILES = $(TEXFILES) $(BIBFILES) plots
ARCHIVFILES += Makefile aa.bst

dvi: $(REV).dvi
pdf: $(REV).pdf

$(REV).dvi: $(TEXFILE) $(TEXFILES) $(BIBFILES)
	    latex $(TEXFILE)
	    bibtex $(AUXFILE)
	    latex $(TEXFILE)
	    latex $(TEXFILE)
	    latex $(TEXFILE)

$(REV).pdf: $(REV).dvi
	    dvips -Ppdf $(REV).dvi
	    ps2pdf $(REV).ps
	    rm $(REV).ps


.PHONY: clean

clean:      
	    @rm -f $(REV).[^t][^e][^xf] *~ *.aux
