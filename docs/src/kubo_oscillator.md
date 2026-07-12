# Kubo Oscillator

```@docs
GeometricProblems.KuboOscillator
```

The deterministic drift alone (no noise) is a plain harmonic oscillator, tracing a circle in phase
space:

```@example
using GeometricProblems.KuboOscillator
using GeometricIntegrators: integrate, ImplicitMidpoint
using CairoMakie

sol = integrate(kubo_oscillator_ode(; timespan = (0.0, 2π)), ImplicitMidpoint())

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "q₁", ylabel = "q₂", aspect = DataAspect(), title = "Kubo oscillator (deterministic drift)")
lines!(ax, sol.q[:, 1], sol.q[:, 2])
fig
```

## Library functions

```@autodocs
Modules = [GeometricProblems.KuboOscillator]
Order   = [:constant, :type, :macro, :function]
```
