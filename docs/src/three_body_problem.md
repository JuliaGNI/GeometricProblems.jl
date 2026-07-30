# Three Body Problem

The default initial condition is the figure-eight choreography of [chenciner2000remarkable](@cite),
in which the three equal masses chase one another along a single closed curve — all three orbits
below lie on top of each other. The default window is one of its periods, so the plot closes.

The 4096-member initial-condition grid of [jin2020sympnets](@cite) is available as
`initial_conditions`, but every one of its members ends in a collision, so those are usable only on a
window that stops short of it — see
[`GeometricProblems.ThreeBody.sympnets_initial_condition`](@ref).

```@example
using GeometricProblems.ThreeBody: hodeproblem
using GeometricIntegrators: ImplicitMidpoint, integrate
using CairoMakie
using LaTeXStrings

morange = RGBf(255 / 256, 127 / 256, 14 / 256)
mred = RGBf(214 / 256, 39 / 256, 40 / 256) 
mpurple = RGBf(148 / 256, 103 / 256, 189 / 256)
mblue = RGBf(31 / 256, 119 / 256, 180 / 256)
mgreen = RGBf(44 / 256, 160 / 256, 44 / 256)

pr = hodeproblem()
sol = integrate(pr, ImplicitMidpoint())
first_body = zeros(2, length(sol.q))
second_body = zeros(2, length(sol.q))
third_body = zeros(2, length(sol.q))

for index in axes(sol.q, 1)
  data_for_present_index = sol.q[index]
  first_body[:, index + 1] = data_for_present_index[1:2]
  second_body[:, index + 1] = data_for_present_index[3:4] 
  third_body[:, index + 1] = data_for_present_index[5:6] 
end

fig = Figure()
ax = Axis(fig[1, 1])
scatter!(ax, first_body, color = mred)
lines!(ax, first_body, color = mred, linestyle = :dash)
scatter!(ax, second_body, color = mblue)
lines!(ax, second_body, color = mblue, linestyle = :dash)
scatter!(ax, third_body, color = mgreen)
lines!(ax, third_body, color = mgreen, linestyle = :dash)

scatter!(ax, first_body[:, 1]', color = :black, label = L"t = 0")
scatter!(ax, second_body[:, 1]', color = :black)
scatter!(ax, third_body[:, 1]', color = :black)
axislegend(position = :rb)

fig
```


## Library functions

```@docs
GeometricProblems.ThreeBody
```

```@autodocs
Modules = [GeometricProblems.ThreeBody]
Order   = [:constant, :type, :macro, :function]
```

## References

```@bibliography
Pages = []

chenciner2000remarkable
jin2020sympnets
```
