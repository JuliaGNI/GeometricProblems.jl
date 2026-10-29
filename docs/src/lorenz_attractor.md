# Lorenz Attractor in 3d

Integrating the default problem traces the well-known strange attractor:

```@example
using GeometricProblems.LorenzAttractor
using GeometricIntegrators: integrate, ImplicitMidpoint
using CairoMakie

sol = integrate(odeproblem(), ImplicitMidpoint())
plot_phase_portrait(sol)
```

```@autodocs
Modules = [GeometricProblems.LorenzAttractor]
Order   = [:constant, :type, :macro, :function]
```
