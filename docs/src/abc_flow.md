# ABC Flow 

The ABC flow (see [hairer2006geometric](@cite)) is described by a divergence-free differential equation whose flow strongly depends on the initial condition. 

```@example
using GeometricIntegrators: integrate, ImplicitMidpoint
using GeometricProblems.ABCFlow
using CairoMakie

ensemble_solution = integrate(odeensemble(), ImplicitMidpoint())

fig = Figure()
ax = Axis3(fig[1, 1]; xlabel = "x", ylabel = "y", zlabel = "z")
for solution in ensemble_solution
    lines!(ax, solution.q[:, 1], solution.q[:, 2], solution.q[:, 3])
end
fig
```

## Library functions

```@docs 
GeometricProblems.ABCFlow
```

```@autodocs
Modules = [GeometricProblems.ABCFlow]
Order = [:constant, :type, :macro, :function]
```

## References

```@bibliography
Pages = []

hairer2006geometric
```
