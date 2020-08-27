@doc raw"""
# Lotka-Volterra model in 2D with "singular" Lagrangian

```math
\begin{aligned}
L (q, \dot{q}) &= \frac{\log q_2}{q_1} \, \dot{q_1} - H(q) , \\
H(q) &= a_1 \, q_1 + a_2 \, q_2 + b_1 \, \log q_1 + b_2 \, \log q_2
\end{aligned}
```

This Lagrangian is equivalent to the Lagrangian of the symmetric Lotka-Volterra model. It
differs only by a gauge transformation with the term ``- 1/2 \, d(\log(q_1) \log(q_2))/dt``.
It leads to the same Euler-Lagrange equations but to a different variational integrator.
"""
module LotkaVolterra2dSingular

    ϑ₁(t, q) = + log(q[2]) / q[1]
    ϑ₂(t, q) = zero(eltype(q))

    dϑ₁dx₁(t, q) = - log(q[2]) / q[1]^2
    dϑ₁dx₂(t, q) = + 1 / (q[1] * q[2])

    dϑ₂dx₁(t, q) = zero(eltype(q))
    dϑ₂dx₂(t, q) = zero(eltype(q))


    # ϑ₁(t, q) = zero(eltype(q))
    # ϑ₂(t, q) = - log(q[1]) / q[2]

    # dϑ₁dx₁(t, q) = zero(eltype(q))
    # dϑ₁dx₂(t, q) = zero(eltype(q))

    # dϑ₂dx₁(t, q) = - 1 / (q[2] * q[1])
    # dϑ₂dx₂(t, q) = + log(q[1]) / q[2]^2


    include("lotka_volterra_2d_common.jl")
    include("lotka_volterra_2d_equations.jl")

end
