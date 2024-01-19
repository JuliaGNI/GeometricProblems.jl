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
using Plots
using GeometricIntegrators
using GeometricProblems.LotkaVolterra2d
using GeometricProblems.LotkaVolterra2dPlots

ode = lotka_volterra_2d_ode()
sol = integrate(ode, Gauss(1))

plot_lotka_volterra_2d(sol, ode)
savefig("lotka_volterra_2d.svg")

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

```@autodocs
Modules = [GeometricProblems.LotkaVolterra2d]
Order   = [:constant, :type, :macro, :function]
```

## Plotting Functions

```@autodocs
Modules = [GeometricProblems.LotkaVolterra2dPlots]
Order   = [:constant, :type, :macro, :function]
```
