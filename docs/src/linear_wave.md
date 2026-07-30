# The Discretized Linear Wave 

The linear wave equation in one dimension has the following Hamiltonian (see e.g. [buchfink2023symplectic](@cite)): 

```math
    \mathcal{H}_\mathrm{cont}(q, p; \mu) := \frac{1}{2}\int_\Omega \mu^2(\partial_\xi q(t, \xi; \mu))^2 + p(t, \xi; \mu)^2 d\xi,
```
where the domain is ``\Omega = (-1/2, 1/2)``. We then divide the domain into ``\tilde{N}`` equidistantly spaced points[^1] ``\xi_i = i\Delta_\xi - 1/2`` for ``i = 1, \ldots, \tilde{N}`` and ``\Delta_\xi := 1/(\tilde{N} + 1)``.

[^1]: In total the system is therefore described by ``N = \tilde{N} + 2`` coordinates, since we also have to consider the two boundary points. The resulting (semi-discrete) Hamiltonian, matching the implementation, then is:

```math
    \mathcal{H}_h(z) = \frac{1}{2} \sum_{i = 1}^{\tilde{N} + 2} p_i^2 + \frac{\mu^2}{4\Delta_\xi^2} \sum_{i = 2}^{\tilde{N} + 1} \left[ (q_i - q_{i - 1})^2 + (q_{i+1} - q_i)^2 \right].
```

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

## Implementation

The equations of motion are written by hand. Because the sum above runs over the *interior* points
``i = 2, \ldots, \tilde{N}+1``, it counts every interior difference twice, so with
``d_k := q_{k+1} - q_k`` the potential is

```math
    V(q) = \frac{c}{2} \sum_{k=1}^{\tilde{N}+1} w_k d_k^2, \qquad
    c := \frac{\mu^2}{2\Delta_\xi^2}, \qquad
    w_k := 2 - [k = 1] - [k = \tilde{N}+1],
```

i.e. the two boundary weights are ``1`` rather than ``2``, and
``\partial V/\partial q_j = c \, (w_{j-1} d_{j-1} - w_j d_j)`` with ``d_0 = d_{\tilde{N}+2} = 0``.
Substituting the weights splits that into a *uniform* weight-two stencil plus four scalar corrections
at the boundaries, which is how `∇V!` is written: the loop is then branch-free, each output is
written exactly once, and the expression is correct for every ``\tilde{N} \ge 1`` without any case
analysis at the edges.

The Lagrangian ``L = \tfrac{1}{2} \sum_k \dot q_k^2 - V(q)`` is *regular* — ``\vartheta = \partial
L/\partial\dot q = \dot q``, so the mass matrix ``M = \partial\vartheta/\partial\dot q`` is the
identity. It is therefore a second-order system of ``n = \tilde{N} + 2`` equations, equivalently first
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