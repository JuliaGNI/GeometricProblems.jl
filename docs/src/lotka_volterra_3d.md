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

fig = Figure(size = (600, 600))
ax = Axis3(fig[1, 1]; xlabel = "q₁", ylabel = "q₂", zlabel = "q₃", title = "Lotka–Volterra 3d")
lines!(ax, sol.q[:, 1], sol.q[:, 2], sol.q[:, 3])
fig
```

```@autodocs
Modules = [GeometricProblems.LotkaVolterra3d]
Order   = [:constant, :type, :macro, :function]
```
