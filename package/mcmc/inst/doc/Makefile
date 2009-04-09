
all : demo.pdf metrop.pdf clean

demo.pdf : demo.tex
	texi2dvi --pdf demo.tex

metrop.pdf : metrop.tex
	texi2dvi --pdf metrop.tex

clean :
	rm -f *.dvi *.aux *.log demo-fig* Rplots.*
