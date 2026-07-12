# Lotka-Volterra 4d

Integrating the default problem gives a quasi-periodic trajectory; here is its projection onto the
``(q_1, q_2)`` plane:

```@example
using GeometricProblems.LotkaVolterra4d
using GeometricIntegrators: integrate, Gauss
using CairoMakie

sol = integrate(odeproblem(), Gauss(2))

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "q₁", ylabel = "q₂", title = "Lotka–Volterra 4d")
lines!(ax, sol.q[:, 1], sol.q[:, 2])
fig
```

```@autodocs
Modules = [GeometricProblems.LotkaVolterra4d]
Order   = [:constant, :type, :macro, :function]
```

