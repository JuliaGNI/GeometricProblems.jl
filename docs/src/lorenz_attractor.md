# Lorenz Attractor in 3d

Integrating the default problem traces the well-known strange attractor:

```@example
using GeometricProblems.LorenzAttractor
using GeometricIntegrators: integrate, ImplicitMidpoint
using CairoMakie

sol = integrate(odeproblem(), ImplicitMidpoint())

fig = Figure(size = (600, 600))
ax = Axis3(fig[1, 1]; xlabel = "x", ylabel = "y", zlabel = "z", title = "Lorenz attractor")
lines!(ax, sol.q[:, 1], sol.q[:, 2], sol.q[:, 3]; linewidth = 0.6)
fig
```

```@autodocs
Modules = [GeometricProblems.LorenzAttractor]
Order   = [:constant, :type, :macro, :function]
```
