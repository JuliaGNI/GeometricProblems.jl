# Lotka-Volterra 2d

Lotka–Volterra models are used in mathematical biology for modelling population dynamics of animal species, as well as many other fields where predator-prey and similar models appear. The dynamics of the growth of two interacting species can be modelled by the following noncanonical Hamiltonian system
```math
\dot{q} = \begin{pmatrix}
\hphantom{-} 0 & + q_1 q_2 \\
- q_1 q_2 & \hphantom{+} 0 \\
\end{pmatrix}
\nabla H (q) ,
\quad
H (q) = a_1 \, q_1 + a_2 \, q_2 + b_1 \, \log q_1 + b_2 \, \log q_2 .
```

```@eval
using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using CairoMakie

ode = odeproblem()
sol = integrate(ode, Gauss(1))

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "q₁", ylabel = "q₂", title = "Lotka–Volterra 2d")
lines!(ax, sol.q[:, 1], sol.q[:, 2])
save("lotka_volterra_2d.svg", fig)

nothing
```

![](lotka_volterra_2d.svg)



## Sub-models

The Euler-Lagrange equations of the Lotka-Volterra model can be obtained from different Lagrangians, which are connected by gauge transformations.
Although they all lead to the same equations of motion, they lead to different variational integrators. Therefore different models based on different Lagrangians are implemented.

```@docs
GeometricProblems.LotkaVolterra2d
```

```@docs
GeometricProblems.LotkaVolterra2dSymmetric
```

```@docs
GeometricProblems.LotkaVolterra2dSingular
```

```@docs
GeometricProblems.LotkaVolterra2dGauge
```


## User Functions

The plotting functions (`plot_solution`, `plot_phase_portrait`, `plot_traces`)
are provided by the `LotkaVolterra2dPlots` extension and become available once a
Makie backend such as `CairoMakie` is loaded.

```@autodocs
Modules = [GeometricProblems.LotkaVolterra2d]
Order   = [:constant, :type, :macro, :function]
```
