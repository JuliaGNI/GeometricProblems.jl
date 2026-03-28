pdf:
	xelatex -shell-escape double-pendulum
	xelatex -shell-escape double-pendulum_dark
	xelatex -shell-escape coupled_harmonic_oscillator
	xelatex -shell-escape coupled_harmonic_oscillator_dark
	xelatex -shell-escape third_degree_spline
	xelatex -shell_escape third_degree_spline_dark

png:
	pdftocairo  -png -r 150 -transp -singlefile  double-pendulum.pdf               		double-pendulum              
	pdftocairo  -png -r 150 -transp -singlefile  double-pendulum_dark.pdf          		double-pendulum_dark       
	pdftocairo  -png -r 250 -transp -singlefile  coupled_harmonic_oscillator.pdf        coupled_harmonic_oscillator              
	pdftocairo  -png -r 250 -transp -singlefile  coupled_harmonic_oscillator_dark.pdf   coupled_harmonic_oscillator_dark  
	pdftocairo  -png -r 150 -transp -singlefile  third_degree_spline.pdf               	third_degree_spline              
	pdftocairo  -png -r 150 -transp -singlefile  third_degree_spline_dark.pdf          	third_degree_spline_dark         

MYDIR := .
convert_with_sips: $(MYDIR)/*.pdf
	for file in $^ ; do \
		sips --setProperty format png --resampleHeightWidthMax 600 $${file} --out $${file%.*}.png ; \
	done

# this is converting pdfs to pngs using pdftocairo (linux version)
convert_with_pdftocairo: $(MYDIR)/*.pdf
	for file in $^ ; do \
		pdftocairo -png -r 250 -transp -singlefile $${file} $${file%.*} ; \
	done

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
