# Planar Point Vortices with Linear One-Form

```@docs
GeometricProblems.PointVorticesLinear
```

Integrating the default two-vortex problem traces the (co-rotating) paths of the two vortices:

```@example
using GeometricProblems.PointVorticesLinear
using GeometricIntegrators: integrate, ImplicitMidpoint
using CairoMakie

sol = integrate(odeproblem(), ImplicitMidpoint())

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Point vortices (linear one-form)")
lines!(ax, sol.q[:, 1], sol.q[:, 2]; label = "vortex 1")
lines!(ax, sol.q[:, 3], sol.q[:, 4]; label = "vortex 2")
axislegend(ax)
fig
```

## Library functions

```@autodocs
Modules = [GeometricProblems.PointVorticesLinear]
Order   = [:constant, :type, :macro, :function]
```
