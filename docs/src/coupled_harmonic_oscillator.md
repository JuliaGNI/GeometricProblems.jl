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
using GeometricIntegrators: integrate, ImplicitMidpoint # hide
using GeometricProblems.CoupledHarmonicOscillator: hodeproblem, default_parameters, tspan, tstep, q₀, p₀
using GeometricEquations: EnsembleProblem # hide
using GeometricMachineLearning: DataLoader # hide
using Plots # hide

const m₁ = default_parameters.m₁
const m₂ = default_parameters.m₂
const k₁ = default_parameters.k₁
const k₂ = default_parameters.k₂
const k = [0.0, 0.5, 0.75, 1.0, 2.0, 3.0, 4.0]

params_collection = [(m₁ = m₁, m₂ = m₂, k₁ = k₁, k₂ = k₂, k = k_val) for k_val in k]
ensemble_problem = EnsembleProblem(hodeproblem().equation, tspan, tstep, (q = q₀, p = p₀), params_collection)
ensemble_solution = integrate(ensemble_problem, ImplicitMidpoint())

dl = DataLoader(ensemble_solution)
q₁ = dl.input.q[1, :, :]
h = tspan[2] / (size(q₁, 1) - 1)
t = 0.0:h:tspan[2]
n_param_sets = length(params_collection)
labels = reshape(["k = "*string(params.k) for params in params_collection], 1, n_param_sets)

const one_plot = false
const psize = (900,600)
plot_q₁ = one_plot ? plot(t, q₁, size=psize) : plot(t, q₁, layout=(n_param_sets, 1), size=psize, label=labels, legend=:topright)

savefig("plot.svg")

nothing
```

![](plot.svg)

## Library functions

```@docs
GeometricProblems.CoupledHarmonicOscillator
```

```@autodocs
Modules = [GeometricProblems.CoupledHarmonicOscillator]
Order   = [:constant, :type, :macro, :function]
```
