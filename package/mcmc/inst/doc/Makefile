
all : demo.pdf metrop.pdf temper.pdf clean

demo.pdf : demo.tex
	texi2dvi --pdf demo.tex

metrop.pdf : metrop.tex
	texi2dvi --pdf metrop.tex

temper.pdf : temper.tex
	texi2dvi --pdf temper.tex

clean :
	rm -f *.dvi *.aux *.log demo-fig* Rplots.*
