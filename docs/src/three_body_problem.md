# Three Body Problem

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

pr = hodeproblem(; tstep = .2)
sol = integrate(pr, ImplicitMidpoint())
first_body = zeros(2, length(sol.s.q))
second_body = zeros(2, length(sol.s.q))
third_body = zeros(2, length(sol.s.q))

for index in axes(sol.s.q, 1)
  data_for_present_index = sol.s.q[index]
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

jin2020sympnets
```