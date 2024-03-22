# GeometricProblems.jl

*GeometricProblems.jl is a collection of ODEs and DAEs with interesting geometric structures
together with useful diagnostics and plotting tools.*

[![PkgEval Status](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricProblems.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/G/GeometricProblems.html)
![CI](https://github.com/JuliaGNI/GeometricProblems.jl/workflows/CI/badge.svg)
[![Build Status](https://travis-ci.org/JuliaGNI/GeometricProblems.jl.svg?branch=main)](https://travis-ci.org/JuliaGNI/GeometricProblems.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaGNI/GeometricProblems.jl/badge.svg)](https://coveralls.io/github/JuliaGNI/GeometricProblems.jl)
[![codecov](https://codecov.io/gh/JuliaGNI/GeometricProblems.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/GeometricProblems.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3740036.svg)](https://doi.org/10.5281/zenodo.3740036)

Typical structures are
* Variational structure, i.e., the equations can defined in terms of a Lagrangian function and be obtained from an action principle;
* Hamiltonian structure, i.e., the equations can be defined in terms of a Hamiltonian function together with a symplectic or Poisson matrix;
* Invariants, i.e., the equations have symmetries and associated conservation laws;
* Volume preservation, i.e., the flow of the equations is divergence-free.


## Contents

```@contents
Pages = [
"diagnostics.md",
"abc_flow.md",
"coupled_harmonic_oscillator.md",
"henon_heiles.md",
"kepler_problem.md",
"lorenz_attractor.md",
"lotka_volterra_2d.md",
"lotka_volterra_3d.md",
"lotka_volterra_4d.md",
"massless_charged_particle.md",
"harmonic_oscillator.md",
"nonlinear_oscillators.md",
"pendulum.md",
"double_pendulum.md",
"point_vortices.md",
"inner_solar_system.md",
"outer_solar_system.md",
"rigid_body.md",
"toda_lattice.md"
]
Depth = 1
```


## References

If you use the figures or implementations provided here, please consider citing GeometricIntegrators.jl as

```
@misc{Kraus:2020:GeometricIntegratorsRepo,
  title={GeometricIntegrators.jl: Geometric Numerical Integration in Julia},
  author={Kraus, Michael},
  year={2020},
  howpublished={\url{https://github.com/JuliaGNI/GeometricIntegrators.jl}},
  doi={10.5281/zenodo.3648325}
}
```

as well as this repository as

```
@misc{Kraus:2020:GeometricProblemsRepo,
  title={GeometricProblems.jl: Collection of Differential Equations with Geometric Structure.},
  author={Kraus, Michael},
  year={2020},
  howpublished={\url{https://github.com/JuliaGNI/GeometricProblems.jl}},
  doi={10.5281/zenodo.4285904}
}
```


## Figure License

> Copyright (c) Michael Kraus <michael.kraus@ipp.mpg.de>
>
> All figures are licensed under the Creative Commons [CC BY-NC-SA 4.0 License](https://creativecommons.org/licenses/by-nc-sa/4.0/).


## Software License

> Copyright (c) Michael Kraus <michael.kraus@ipp.mpg.de>
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.
