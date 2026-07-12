# Inner Solar System

An ``N``-body gravitational model of the inner solar system (Sun and the terrestrial planets),
governed by the Hamiltonian
```math
H(q, p) = \sum_i \frac{\| p_i \|^2}{2 m_i}
          - G \sum_{i < j} \frac{m_i m_j}{\| q_i - q_j \|} .
```

!!! note "Not yet implemented"
    This system is documented for completeness but is **not yet provided** as a problem type in
    `GeometricProblems.jl`. For a small gravitational ``N``-body example, see the
    [Three Body Problem](three_body_problem.md).
