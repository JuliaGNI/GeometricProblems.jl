# Massless Charged Particle

```@docs
GeometricProblems.MasslessChargedParticle
```

```@eval
using Plots
using GeometricIntegrators
using GeometricProblems.MasslessChargedParticle
using GeometricProblems.MasslessChargedParticlePlots

ode = odeproblem()
sol = integrate(ode, Gauss(1))

plot_massless_charged_particle(sol, ode)
savefig("massless_charged_particle.svg")

nothing
```

![](massless_charged_particle.svg)



```@autodocs
Modules = [GeometricProblems.MasslessChargedParticle]
Order   = [:constant, :type, :macro, :function]
```
