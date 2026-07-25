# Planar Point Vortices

Integrating the default two-vortex problem traces the paths of the two vortices in the plane:

```@example
using GeometricProblems.PointVortices
using GeometricIntegrators: integrate, ImplicitMidpoint
using CairoMakie

sol = integrate(odeproblem(), ImplicitMidpoint())

plot_phase_portrait(sol)
```

Both point-vortex modules conserve the energy `:h` and the angular momentum `:p`, so the
generic [diagnostic plots](diagnostics.md) work on their solutions for either invariant.

```@autodocs
Modules = [GeometricProblems.PointVortices]
Order   = [:constant, :type, :macro, :function]
```
