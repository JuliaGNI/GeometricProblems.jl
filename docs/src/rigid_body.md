# The rigid body 

```@example 
using GeometricProblems.RigidBody: odeensemble
using GeometricIntegrators: integrate, ImplicitMidpoint
using GeometricEquations: EnsembleProblem
using GeometricSolutions: GeometricSolution
using GLMakie

ics = [
        [sin(1.1), 0., cos(1.1)],
        [sin(2.1), 0., cos(2.1)],
        [sin(2.2), 0., cos(2.2)],
        [0., sin(1.1), cos(1.1)],
        [0., sin(1.5), cos(1.5)], 
        [0., sin(1.6), cos(1.6)]
        ]

ensemble_problem = odeensemble(ics)
ensemble_solution = integrate(ensemble_problem, ImplicitMidpoint())

function plot_geometric_solution!(p, solution::GeometricSolution; kwargs...)
    lines!(p, solution.q[:, 1].parent, solution.q[:, 2].parent, solution.q[:, 3].parent; kwargs...)
end

function sphere(r, C)   # r: radius; C: center [cx,cy,cz]
           n = 100
           u = range(-π, π; length = n)
           v = range(0, π; length = n)
           x = C[1] .+ r * cos.(u) * sin.(v)'
           y = C[2] .+ r * sin.(u) * sin.(v)'
           z = C[3] .+ r * ones(n) * cos.(v)'
           return x, y, z
       end

fig, ax, plt = surface(sphere(1., [0., 0., 0.])..., alpha = .6)

for (i, solution) in zip(1:length(ensemble_solution), ensemble_solution.s)
    plot_geometric_solution!(ax, solution; label = "trajectory "*string(i), linewidth=2)
end

ax
```


## Library functions

```@docs 
GeometricProblems.RigidBody
```

```@autodocs
Modules = [GeometricProblems.RigidBody]
Order = [:constant, :type, :macro, :function]
```

```@bibliography
Pages = []

bajars2023locally
```