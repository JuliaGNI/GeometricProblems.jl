# The Discretized Linear Wave 

The linear wave equation in one dimension has the following Hamiltonian (see e.g. [buchfink2023symplectic](@cite)): 

```math
    \mathcal{H}_\mathrm{cont}(q, p; \mu) := \frac{1}{2}\int_\Omega \mu^2(\partial_\xi q(t, \xi; \mu))^2 + p(t, \xi; \mu)^2 d\xi,
```
where the domain is ``\Omega = (-1/2, 1/2)``. We then divide the interior of the domain into ``\tilde{N}`` equidistantly spaced points ``\xi_i = i\Delta_\xi - 1/2`` for ``i = 1, \ldots, \tilde{N}`` and ``\Delta_\xi := 1/(\tilde{N} + 1)``. The homogeneous Dirichlet boundary conditions fix ``q_0 = q_{\tilde{N}+1} = 0``, so the two boundary points are not part of the state, which has ``\tilde{N}`` positions and ``\tilde{N}`` momenta. The resulting (semi-discrete) Hamiltonian is the one of [buchfink2023symplectic](@cite):

```math
    \mathcal{H}_h(q, p) = \frac{1}{2} \sum_{i = 1}^{\tilde{N}} p_i^2 + \frac{\mu^2}{2\Delta_\xi^2} \sum_{i = 0}^{\tilde{N}} (q_{i+1} - q_i)^2
    = \frac{1}{2} p^Tp - \frac{\mu^2}{2} q^T D_{\xi\xi} q,
```

where ``D_{\xi\xi} = \mathrm{tridiag}(1, -2, 1)/\Delta_\xi^2`` is the central-difference matrix, and the equations of motion ``\dot{q} = p``, ``\dot{p} = \mu^2 D_{\xi\xi} q`` are Hamilton's equations with the canonical Poisson matrix.

The discretized linear wave equation is an example of a *completely integrable system*, i.e. a Hamiltonian system evolving in ``\mathbb{R}^{2n}`` that has ``n`` Poisson-commuting invariants of motion (see [arnold1978mathematical](@cite)). 

For evaluating the system we specify the following initial[^2] and boundary conditions: 

```math
\begin{aligned}
	q_0(\omega;\mu) := & q(0, \omega; \mu) \\ 
	p(0, \omega; \mu) = \partial_tq(0,\xi;\mu) = & -\mu\partial_\omega{}q_0(\xi;\mu) \\
	q(t,\omega;\mu) = & 0, \text{ for } \omega\in\partial\Omega.
\end{aligned}
```

[^2]: The precise shape of ``q_0(\cdot;\cdot)`` is described in [the chapter on initial conditions](initial_condition.md).

## Solution

By default `GeometricProblems` uses the following parameters:

```@example linear_wave
using GeometricIntegrators
using CairoMakie
import GeometricProblems.LinearWave as lw

lw.default_parameters()
```

And if we integrate we get:

```@example linear_wave
# A smaller lattice than the default Ñ = 256. The figure is indistinguishable, but an implicit step
# solves a dense 2n-dimensional nonlinear system, so the cost of the integration grows as n³.
N = 128

sol = integrate(lw.hodeproblem(N), ImplicitMidpoint())

Ω  = lw.interior_domain(N)        # the N interior points; q = 0 at ξ = ±1/2
nt = length(sol.t) - 1             # the time axis of a DataSeries runs 0:nt

fig = Figure(size = (800, 700))

ax1 = Axis(fig[1, 1]; ylabel = "q", title = "Snapshots")
for step in round.(Int, range(0, nt, length = 6))   # `0:nt÷5:nt` would stop short of the last step
    lines!(ax1, Ω, sol.q[step, :]; label = "t = $(round(sol.t[step]; digits = 2))")
end
axislegend(ax1; position = :lt, labelsize = 10)

# q over the whole (ξ, t) plane: a single band with one slope, which is the one-directional travel
# that the snapshots above only hint at.
Q = [sol.q[n][i] for i in eachindex(Ω), n in 0:nt]

ax2 = Axis(fig[2, 1]; xlabel = "ξ", ylabel = "t", title = "Space-time")
hm = heatmap!(ax2, Ω, collect(sol.t), Q)
Colorbar(fig[2, 2], hm; label = "q")

linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false)

fig
```

As we can see the pulse travels in one direction, at the speed ``\mu`` set by the initial momentum
``p(0, \omega; \mu) = -\mu \partial_\omega q_0``, and keeps its shape: the discretization is
non-dispersive to the resolution of the figure.

## Implementation

The equations of motion are written by hand. The force ``-\partial V/\partial q = \mu^2 D_{\xi\xi} q``
is the three-point stencil ``\mu^2 (q_{j-1} - 2q_j + q_{j+1})/\Delta_\xi^2``, in which the Dirichlet
values ``q_0 = q_{\tilde{N}+1} = 0`` drop one neighbor from the first and the last row; this is how
`∇V!` is written.

The Lagrangian ``L = \tfrac{1}{2} \sum_k \dot q_k^2 - V(q)`` is *regular*: ``\vartheta = \partial
L/\partial\dot q = \dot q``, so the mass matrix ``M = \partial\vartheta/\partial\dot q`` is the
identity. It is therefore a second-order system of ``n = \tilde{N}`` equations, equivalently first
order in ``2n``, and its Lagrange two-form is the ``2n \times 2n`` canonical
``\omega = [\,0\; -I;\; I\; 0\,]``.

Passing `symbolic = true` to any of the four constructors generates the equations of motion with
EulerLagrange instead. The two agree to round-off and the tests check that they do, but the symbolic
route does not scale: EulerLagrange builds ``\omega`` by differentiating a dense ``2n \times 2n``
matrix, so at the default ``\tilde{N} = 256`` `lodeproblem(; symbolic = true)` takes 155 s and emits
14 MB of code for a two-form that no integrator evaluates. `benchmark/linear_wave.jl` has the
measurements.

## Library functions

```@docs
GeometricProblems.LinearWave
```

```@autodocs
Modules = [GeometricProblems.LinearWave]
Order   = [:constant, :type, :macro, :function]
```

```@bibliography
Pages = []
 
buchfink2023symplectic
```