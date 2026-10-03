
.PHONY: documenter

all: documenter

doctests: 
	cd ..; julia --project=docs  -e ';									\
		using Pkg;														\
		Pkg.develop(PackageSpec(path=pwd()));							\
		Pkg.instantiate();												\
		Pkg.build();													\
		Pkg.precompile();												\
		using Documenter: doctest;										\
		using GeometricProblems;									\
		doctest(GeometricProblems)'					

documenter:
	julia --color=yes --project=. --threads=8 make.jl

clean:
	rm -Rf build
