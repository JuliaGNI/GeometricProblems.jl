# The Discretized Linear Wave 

The linear wave equation in one dimension has the following Hamiltonian (see e.g. [buchfink2023symplectic](@cite)): 

```math
    \mathcal{H}_\mathrm{cont}(q, p; \mu) := \frac{1}{2}\int_\Omega \mu^2(\partial_\xi q(t, \xi; \mu))^2 + p(t, \xi; \mu)^2 d\xi,
```
where the domain is ``\Omega = (-1/2, 1/2)``. We then divide the domain into ``\tilde{N}`` equidistantly spaces points[^1] ``\xi_i = i\Delta_\xi - 1/2`` for ``i = 1, \ldots, \tilde{N}`` and ``\Delta_xi := 1/(\tilde{N} + 1)``.

[^1]: In total the system is therefore described by ``N = \tilde{N} + 2`` coordinates, since we also have to consider the boundary. The resulting Hamiltonian then is:

```math
    \mathcal{H}_h(z) = \sum_{i = 1}^{\tilde{N}}\frac{\Delta{}x}{2}\left[ p_i^2 + \mu^2 \frac{(q_i - q_{i - 1})^2 + (q_{i+1} - q_i)^2}{2\Delta{}x} \right].
```

The discretized linear wave equation example of an *completely-integrable system*, i.e. a Hamiltonian system evolving in ``\mathbb{R}^{2n}`` that has ``n`` Poisson-commuting invariants of motion (see [arnold1978mathematical](@cite)). 

For evaluating the system we specify the following initial[^2] and boundary conditions: 

```math
\begin{aligned}
	q_0(\omega;\mu) := & q(0, \omega; \mu) \\ 
	p(0, \omega; \mu) = \partial_tq(0,\xi;\mu) = & -\mu\partial_\omega{}q_0(\xi;\mu) \\
	q(t,\omega;\mu) = & 0, \text{ for } \omega\in\partial\Omega.
\end{aligned}
```

[^2]: The precise shape of ``q_0(\cdot;\cdot)`` is described in [the chapter on initial conditions](initial_condition.md).


By default `GeometricProblems` uses the following parameters: 
```@example linear_wave
using GeometricIntegrators, Plots # hide
import GeometricProblems.LinearWave as lw

lw.default_parameters
```

And if we integrate we get: 
```@example linear_wave
problem = lw.hodeproblem() 
sol = integrate(problem, ImplicitMidpoint())

# plot 6 time steps
time_steps = 0 : (length(sol.t) - 1)  ÷ 5 : (length(sol.t) - 1)
p = plot()
for time_step in time_steps
    plot!(p, lw.get_domain(lw.Ñ + 2), sol.q[time_step, :], label = "t = "*string(round(sol.t[time_step]; digits = 2)))
end

p
```

As we can see the thin pulse travels in one direction. 

## Library functions

```@docs
GeometricProblems.LinearWave
```

```@bibliography
Pages = []
 
buchfink2023symplectic
```