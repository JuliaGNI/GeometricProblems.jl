# Lotka-Volterra 3d

```@docs
GeometricProblems.LotkaVolterra3d
```

Integrating the default problem traces a trajectory on the Casimir surface:

```@example
using GeometricProblems.LotkaVolterra3d
using GeometricIntegrators: integrate, Gauss
using CairoMakie

sol = integrate(odeproblem(), Gauss(2))
plot_phase_portrait(sol)
```

```@autodocs
Modules = [GeometricProblems.LotkaVolterra3d]
Order   = [:constant, :type, :macro, :function]
```
