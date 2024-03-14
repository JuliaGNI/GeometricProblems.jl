
.PHONY: documenter images

all: images install_gml precompile documenter remove_gml

precompile:
	cd .. ; julia --project=docs -e '; \
		using Pkg; \
		Pkg.develop(PackageSpec(path=pwd())); \
		Pkg.instantiate(); \
		Pkg.build(); \
		Pkg.precompile()' ;
	cd docs

test_docs: 
	cd ..; julia --project=docs  -e ';									\
		using Pkg;														\
		Pkg.develop(PackageSpec(path=pwd()));							\
		Pkg.instantiate();												\
		Pkg.build();													\
		Pkg.precompile();												\
		using Documenter: doctest;										\
		using GeometricProblems;									\
		doctest(GeometricProblems)'					

install_gml:
	cd .. ;
	julia --project=. -e '; \
		using Pkg; \
		Pkg.add(url="https://github.com/JuliaGNI/GeometricMachineLearning.jl", rev="volume_preserving_transformer")'
	cd docs

remove_gml: 
	cd .. ; julia --project=. -e '; \
		using Pkg; \
		Pkg.rm("GeometricMachineLearning")'

documenter:
	julia --color=yes --project=. make.jl

images:
	$(MAKE) all   -C src/images
	$(MAKE) clean -C src/images

clean:
	$(MAKE) empty -C src/images
	rm -Rf build
