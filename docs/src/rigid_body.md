# The rigid body 

```@example 
using GeometricProblems.RigidBody: odeproblem
using GeometricIntegrators: integrate, ImplicitMidpoint
using Plots; pyplot()

sol = integrate(odeproblem(), ImplicitMidpoint())

function sphere(r, C)   # r: radius; C: center [cx,cy,cz]
           n = 100
           u = range(-π, π; length = n)
           v = range(0, π; length = n)
           x = C[1] .+ r*cos.(u) * sin.(v)'
           y = C[2] .+ r*sin.(u) * sin.(v)'
           z = C[3] .+ r*ones(n) * cos.(v)'
           return x, y, z
       end

surface(sphere(1., [0., 0., 0.]), alpha = .2, legend = false)

plot!(sol.q[:, 1], sol.q[:, 2], sol.q[:, 3], color = 1, label = "trajectory")
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