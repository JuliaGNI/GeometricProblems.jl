# The Discretized Linear Wave 

The discretized linear wave equation example of an *completely-integrable system*, i.e. a Hamiltonian system evolving in ``\mathbb{R}^{2n}`` that has ``n`` Poisson-commuting invariants of motion (see [arnold1978mathematical](@cite)). 

The Hamiltonian takes the following form: 

```math
    H(q, p) = \sum_{n\in\mathbb{Z}}\left(  \frac{p_n^2}{2} + \alpha (q_n - q_{n+1}) ^ 2 \right).
```

In practice we impose periodic boundary conditions: 
```math
\begin{aligned}
    q_{n+N} &  \equiv q_n \\ 
    p_{n+N} &   \equiv p_n.
\end{aligned}
```

Hence we have: 

```math 
    H(q, p) = \sum_{n=1}^{N-1} \left(  \frac{p_n^2}{2} + \alpha (q_n - q_{n+1}) ^ 2 \right) + \frac{p_N^2}{2} + \alpha (q_N - q_1)^2.
```

We can model the evolution of a thin pulse in this system:

```@example
using GeometricProblems, GeometricIntegrators, Plots # hide

problem = GeometricProblems.LinearWave.hodeproblem() 
sol = integrate(problem, ImplicitMidpoint())

time_steps = (0, 50, 150, 190)
p = plot()
for time_step in time_steps
    plot!(p, sol.q[time_step, :], label = "t = $(sol.t[time_step])")
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