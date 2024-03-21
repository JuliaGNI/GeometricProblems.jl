# ABC Flow 

The ABC flow (see [hairer2006geometric](@cite)) is described by a divergence-free differential equation whose flow strongly depends on the initial condition. 

```@example
using GeometricIntegrators: integrate, ImplicitMidpoint
using GeometricProblems.ABCFlow
using GeometricEquations: EnsembleProblem
using GeometricSolutions: GeometricSolution
using Plots

ensemble_solution = integrate(odeensemble(), ImplicitMidpoint())

function plot_geometric_solution!(p::Plots.Plot, solution::GeometricSolution; kwargs...)
    plot!(p, solution.q[:, 1], solution.q[:, 2], solution.q[:, 3]; kwargs...)
end

p = plot()
for solution in ensemble_solution
    plot_geometric_solution!(p, solution)
end
p
```

## Library functions

```@docs 
GeometricProblems.ABCFlow
```

```@autodocs
Modules = [GeometricProblems.ABCFlow]
Order = [:constant, :type, :macro, :function]
```

```@bibliography
Pages = []

hairer2006geometric
```