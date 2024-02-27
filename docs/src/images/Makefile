pdf:
	xelatex -shell-escape double-pendulum
	xelatex -shell-escape double-pendulum_dark
	xelatex -shell-escape coupled_harmonic_oscillator
	xelatex -shell-escape coupled_harmonic_oscillator_dark

png:
	pdftocairo  -png -r 150 -transp -singlefile  double-pendulum.pdf               		double-pendulum              
	pdftocairo  -png -r 150 -transp -singlefile  double-pendulum_dark.pdf          		double-pendulum_dark       
	pdftocairo  -png -r 250 -transp -singlefile  coupled_harmonic_oscillator.pdf        coupled_harmonic_oscillator              
	pdftocairo  -png -r 250 -transp -singlefile  coupled_harmonic_oscillator_dark.pdf   coupled_harmonic_oscillator_dark       

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
