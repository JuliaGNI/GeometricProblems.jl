pdf:
	xelatex -shell-escape double-pendulum
	xelatex -shell-escape double-pendulum_dark

png:
	pdftocairo  -png -r 150 -transp -singlefile  double-pendulum.pdf               double-pendulum              
	pdftocairo  -png -r 150 -transp -singlefile  double-pendulum_dark.pdf          double-pendulum_dark       

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
