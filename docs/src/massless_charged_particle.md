# Massless Charged Particle

```@docs
GeometricProblems.MasslessChargedParticle
```

```@eval
using GeometricIntegrators
using GeometricProblems.MasslessChargedParticle
using CairoMakie

ode = odeproblem()
sol = integrate(ode, Gauss(1))

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "x₁", ylabel = "x₂", title = "Massless charged particle trajectory")
lines!(ax, sol.q[:, 1], sol.q[:, 2])
save("massless_charged_particle.svg", fig)

nothing
```

![](massless_charged_particle.svg)



```@autodocs
Modules = [GeometricProblems.MasslessChargedParticle]
Order   = [:constant, :type, :macro, :function]
```
