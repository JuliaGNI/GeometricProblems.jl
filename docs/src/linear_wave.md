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



## Library functions

```@docs
GeometricProblems.LinearWave
```

```@bibliography
Pages = []
 
buchfink2023symplectic
```