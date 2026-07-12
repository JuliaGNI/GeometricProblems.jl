# Coupled Harmonic Oscillator 

This system describes two harmonic oscillators that are coupled nonlinearly. 

```@example
HTML("""<object type="image/svg+xml" class="display-light-only" data=$(joinpath(Main.buildpath, "images/coupled_harmonic_oscillator.png"))></object>""") # hide
```

```@example
HTML("""<object type="image/svg+xml" class="display-dark-only" data=$(joinpath(Main.buildpath, "images/coupled_harmonic_oscillator_dark.png"))></object>""") # hide
```

The following shows the ``q_1`` component of the system for different values of ``k``: 

```@eval
using GeometricIntegrators: integrate, ImplicitMidpoint 
using GeometricProblems.CoupledHarmonicOscillator: hodeensemble, default_parameters
using CairoMakie
 
const m₁ = default_parameters.m₁  
const m₂ = default_parameters.m₂ 
const k₁ = default_parameters.k₁ 
const k₂ = default_parameters.k₂ 
const k = [0.0, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0] 
 
params_collection = [(m₁ = m₁, m₂ = m₂, k₁ = k₁, k₂ = k₂, k = k_val) for k_val in k] 
# ensemble problem
ep = hodeensemble(; parameters = params_collection)
ensemble_solution = integrate(ep, ImplicitMidpoint())
 
t = ensemble_solution.t

q₁ = zeros(1, length(t), length(k))

for index in axes(k, 1)
    q₁[1, :, index] =  ensemble_solution.s[index].q[:, 1]
end

q₁ = q₁[1, :, :]

tvals = range(0.0, 100.0, length = size(q₁, 1))

fig = Figure(size = (900, 600))
for index in axes(k, 1)
    ax = Axis(fig[index, 1]; ylabel = "q₁", title = "k = $(k[index])", titlesize = 10)
    lines!(ax, tvals, q₁[:, index])
end

save("q_component.png", fig)

nothing
```

![](q_component.png)


## Library functions

```@docs
GeometricProblems.CoupledHarmonicOscillator
```

```@autodocs
Modules = [GeometricProblems.CoupledHarmonicOscillator]
Order   = [:constant, :type, :macro, :function]
```
