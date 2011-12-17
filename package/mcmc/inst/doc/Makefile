
all : bfst.pdf demo.pdf metrop.pdf temper.pdf debug.pdf clean

demo.tex : demo.Rnw
	$(R_HOME)/bin/R CMD Sweave demo.Rnw

demo.pdf : demo.tex
	pdflatex demo.tex
	pdflatex demo.tex

metrop.pdf : metrop.tex
	pdflatex metrop.tex
	pdflatex metrop.tex

temper.pdf : temper.tex
	pdflatex temper.tex
	pdflatex temper.tex

debug.tex : debug.Rnw
	$(R_HOME)/bin/R CMD Sweave debug.Rnw

debug.pdf : debug.tex
	pdflatex debug.tex
	pdflatex debug.tex

bfst.tex : bfst.Rnw
	$(R_HOME)/bin/R CMD Sweave bfst.Rnw

bfst.pdf : bfst.tex
	pdflatex bfst.tex
	pdflatex bfst.tex

clean :
	rm -f *.dvi *.aux *.log demo-fig* Rplots.*
