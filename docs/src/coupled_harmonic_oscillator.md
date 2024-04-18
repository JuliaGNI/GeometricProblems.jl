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
using Plots 
 
const m₁ = default_parameters.m₁  
const m₂ = default_parameters.m₂ 
const k₁ = default_parameters.k₁ 
const k₂ = default_parameters.k₂ 
const k = [0.0, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0] 
 
params_collection = [(m₁ = m₁, m₂ = m₂, k₁ = k₁, k₂ = k₂, k = k_val) for k_val in k] 
# ensemble problem
ep = hodeensemble(; params = params_collection)
ensemble_solution = [integrate(problem, ImplicitMidpoint()) for problem in problems]
 
t = ensemble_solution[1].t

q₁ = zeros(1, length(t), length(k))

for index in axes(k, 1)
    q₁[1, :, index] =  ensemble_solution[index].q[:, 1]
end

n_param_sets = length(params_collection) #hide 
labels = reshape(["k = "*string(params.k) for params in params_collection], 1, n_param_sets) 
 
q₁ = q₁[1, :, :]
const one_plot = false 
const psize = (900, 600) 
plot_q₁ = one_plot ? plot(0.0:0.4:100.0, q₁, size=psize) : plot(0.0:0.4:100.0, q₁, layout=(n_param_sets, 1), size=psize, label=labels, legend=:topright)

png(plot_q₁, "q_component")

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
