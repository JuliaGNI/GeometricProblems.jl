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

fig = plot_phase_portrait(sol)
save("massless_charged_particle.svg", fig)

nothing
```

![](massless_charged_particle.svg)



```@autodocs
Modules = [GeometricProblems.MasslessChargedParticle]
Order   = [:constant, :type, :macro, :function]
```


## Singular vector potential

```@docs
GeometricProblems.MasslessChargedParticleSingular
```

```@eval
using GeometricIntegrators
using GeometricProblems.MasslessChargedParticleSingular
using CairoMakie

ode = odeproblem()
sol = integrate(ode, Gauss(1))

fig = plot_phase_portrait(sol)
save("massless_charged_particle_singular.svg", fig)

nothing
```

![](massless_charged_particle_singular.svg)



```@autodocs
Modules = [GeometricProblems.MasslessChargedParticleSingular]
Order   = [:constant, :type, :macro, :function]
```
