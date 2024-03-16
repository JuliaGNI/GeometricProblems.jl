# The rigid body 

```@example 
using GeometricProblems.RigidBody: odeproblem, tspan, tstep, default_parameters
using GeometricIntegrators: integrate, ImplicitMidpoint
using GeometricEquations: EnsembleProblem
using GeometricSolutions: GeometricSolution
using Plots; pyplot()

ics = [
        (q = [sin(1.1), 0., cos(1.1)], ),
        (q = [sin(2.1), 0., cos(2.1)], ),
        (q = [sin(2.2), 0., cos(2.2)], ),
        (q = [0., sin(1.1), cos(1.1)], ),
        (q = [0., sin(1.5), cos(1.5)], ), 
        (q = [0., sin(1.6), cos(1.6)], )
]

ensemble_problem = EnsembleProblem(odeproblem().equation, tspan, tstep, ics, default_parameters)
ensemble_solution = integrate(ensemble_problem, ImplicitMidpoint())

function plot_geometric_solution!(p::Plots.Plot, solution::GeometricSolution; kwargs...)
    plot!(p, solution.q[:, 1], solution.q[:, 2], solution.q[:, 3]; kwargs...)
end

function sphere(r, C)   # r: radius; C: center [cx,cy,cz]
           n = 100
           u = range(-π, π; length = n)
           v = range(0, π; length = n)
           x = C[1] .+ r*cos.(u) * sin.(v)'
           y = C[2] .+ r*sin.(u) * sin.(v)'
           z = C[3] .+ r*ones(n) * cos.(v)'
           return x, y, z
       end

p = surface(sphere(1., [0., 0., 0.]), alpha = .2, colorbar = false)

for (i, solution) in zip(1:length(ensemble_solution), ensemble_solution)
    plot_geometric_solution!(p, solution; color = i, label = "trajectory "*string(i))
end

p
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