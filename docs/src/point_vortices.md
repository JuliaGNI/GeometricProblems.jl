# Planar Point Vortices

Integrating the default two-vortex problem traces the paths of the two vortices in the plane:

```@example
using GeometricProblems.PointVortices
using GeometricIntegrators: integrate, ImplicitMidpoint
using CairoMakie

sol = integrate(odeproblem(), ImplicitMidpoint())

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y", aspect = DataAspect(), title = "Point vortices")
lines!(ax, sol.q[:, 1], sol.q[:, 2]; label = "vortex 1")
lines!(ax, sol.q[:, 3], sol.q[:, 4]; label = "vortex 2")
axislegend(ax)
fig
```

```@autodocs
Modules = [GeometricProblems.PointVortices]
Order   = [:constant, :type, :macro, :function]
```
