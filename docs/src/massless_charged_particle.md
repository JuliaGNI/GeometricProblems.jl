
# Massless Charged Particle

```@docs
GeometricProblems.MasslessChargedParticle
```

```@eval
using Plots
using GeometricIntegrators
using GeometricProblems.MasslessChargedParticle
using GeometricProblems.MasslessChargedParticle: Δt, nt

ode = massless_charged_particle_ode()
sol = integrate(ode, getTableauGLRK(1), Δt, nt)

plot_massless_charged_particle(sol, ode.parameters)
savefig("massless_charged_particle.svg")

nothing
```

![](massless_charged_particle.svg)



```@autodocs
Modules = [GeometricProblems.MasslessChargedParticle]
Order   = [:constant, :type, :macro, :function]
```
