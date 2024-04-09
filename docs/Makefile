
.PHONY: documenter images

all: images documenter

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
	julia --color=yes --project=. make.jl

images:
	$(MAKE) all   -C src/images
	$(MAKE) clean -C src/images
	
clean:
	$(MAKE) empty -C src/images
	rm -Rf build
