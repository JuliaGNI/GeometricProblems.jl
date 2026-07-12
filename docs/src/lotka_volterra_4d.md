# Lotka-Volterra 4d

Integrating the default problem gives a quasi-periodic trajectory; here is its projection onto the
``(q_1, q_2)`` plane:

```@example
using GeometricProblems.LotkaVolterra4d
using GeometricIntegrators: integrate, Gauss
using CairoMakie

sol = integrate(odeproblem(), Gauss(2))
plot_phase_portrait(sol)
```

```@autodocs
Modules = [GeometricProblems.LotkaVolterra4d]
Order   = [:constant, :type, :macro, :function]
```

