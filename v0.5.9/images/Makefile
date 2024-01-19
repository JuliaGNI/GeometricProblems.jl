pdf:
	xelatex -shell-escape double-pendulum

png:
	pdftocairo  -png -r 150 -transp -singlefile  double-pendulum.pdf               double-pendulum              
clean:
	rm -f *.aux
	rm -f *.fdb_latexmk
	rm -f *.fls
	rm -f *.log
	rm -f *.pdf
	rm -f *.xdv

empty: clean
	rm -f *.png
	rm -f *.svg

all: pdf png
