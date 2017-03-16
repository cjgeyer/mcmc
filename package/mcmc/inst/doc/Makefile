
all : metrop.pdf temper.pdf clean

metrop.pdf : metrop.tex
	pdflatex metrop.tex
	pdflatex metrop.tex

temper.pdf : temper.tex
	pdflatex temper.tex
	pdflatex temper.tex

clean :
	rm -f *.dvi *.aux *.log

