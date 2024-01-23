# Toda Lattice 

The Toda lattice is a prime example of an *completely-integrable system*, i.e. a Hamiltonian system evolving in ``\mathbb{R}^{2n}`` that has ``n`` Poisson-commuting invariants of motion (see [arnol1978mathematical](@cite)). It is named after Morikazu Toda who used it to model a one-dimensional crystal [toda1967vibration](@cite).

The Hamiltonian of the Toda lattice takes the following form: 

```math
    H(q, p) = \sum_{n\in\mathbb{Z}}\left(  \frac{p_n^2}{2} + \alpha e^{q_n - q_{n+1}} \right).
```

In practice we work with a finite number of particles ``N`` and impose periodic boundary conditions: 
```math
\begin{aligned}
    q_{n+N} & = \equiv q_n \\ 
    p_{n+N} & = \equiv p_n.
\end{aligned}
```

Hence we have: 

```math 
    H(q, p) = \sum_{n=1}^{N-1} \left(  \frac{p_n^2}{2} + \alpha e^{q_n - q_{n+1}} \right) + \frac{p_N^2}{2} + \alpha e^{q_N - q_1}.
```

## Library functions

```@docs
GeometricProblems.TodaLattice
```