# Toda Lattice 

The Toda lattice is a prime example of an *completely-integrable system*, i.e. a Hamiltonian system evolving in ``\mathbb{R}^{2n}`` that has ``n`` Poisson-commuting invariants of motion (see [arnold1978mathematical](@cite)). It is named after Morikazu Toda who used it to model a one-dimensional crystal [toda1967vibration](@cite).

The Hamiltonian of the Toda lattice takes the following form: 

```math
    H(q, p) = \sum_{n\in\mathbb{Z}}\left(  \frac{p_n^2}{2} + \alpha e^{q_n - q_{n+1}} \right).
```

In practice we work with a finite number of particles ``N`` and impose periodic boundary conditions: 
```math
\begin{aligned}
    q_{n+N} &  \equiv q_n \\ 
    p_{n+N} &   \equiv p_n.
\end{aligned}
```

Hence we have: 

```math 
    H(q, p) = \sum_{n=1}^{N-1} \left(  \frac{p_n^2}{2} + \alpha e^{q_n - q_{n+1}} \right) + \frac{p_N^2}{2} + \alpha e^{q_N - q_1}.
```

We can model the evolution of a thin pulse in this system:

```julia
using GeometricProblems, GeometricIntegrators, GLMakie # hide

problem = GeometricProblems.TodaLattice.hodeproblem(; tspan = (0.0, 2000.)) 
sol = integrate(problem, ImplicitMidpoint())

time_steps = 0:10:length(sol.q)

fig = Figure()
ax = Axis(fig[1, 1])
mblue = RGBf(31 / 256, 119 / 256, 180 / 256)
lines!(ax, sol.q[0, :], label = "t = $(sol.t[0])", color = mblue)
framerate = 30
mblue = 
record(fig, "toda_animation.mp4", time_steps;
    framerate = framerate) do time_step
    empty!(ax)
    lines!(ax, sol.q[time_step, :], label = "t = $(sol.t[time_step])", color = mblue)
    ylims!(ax, 0., 1.)
    axislegend(ax; position = (1.01, 1.5), labelsize = 8)
end
Docs.HTML("""<video mute autoplay loop controls src="toda_animation.mp4" />""")
```

As we can see the thin pulse separates into two smaller pulses an they start traveling in opposite directions until they meet again at time ``t\approx120``. But it is important to note that the right peak at time ``120`` is below the one at time ``0``. This is not a numerical artifact but a feature of the Toda lattice! 

## Library functions

```@docs
GeometricProblems.TodaLattice
```

```@bibliography
Pages = []

arnold1978mathematical 
toda1967vibration
```