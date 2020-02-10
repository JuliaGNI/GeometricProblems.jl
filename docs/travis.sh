#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --project --color=yes -e 'using Pkg; Pkg.instantiate(); import GeometricProblems; include(joinpath(dirname(pathof(GeometricProblems)), "..", "docs", "make.jl"))';
fi
