@doc raw"""
# Lotka-Volterra model in 2D with symmetric Lagrangian with gauge term

```math
\begin{aligned}
L (q, \dot{q}) &= \bigg( q_2 + \frac{1}{2} \frac{\log q_2}{q_1} \bigg) \, \dot{q_1} + \bigg( q_1 - \frac{1}{2} \frac{\log q_1}{q_2} \bigg) \, \dot{q_2} - H(q) , \\
H(q) &= a_1 \, q_1 + a_2 \, q_2 + b_1 \, \log q_1 + b_2 \, \log q_2
\end{aligned}
```

This Lagrangian is equivalent to the Lagrangian of the symmetric Lotka-Volterra model. It
differs only by a gauge transformation with the term ``d(q_1 q_2)/dt``. It leads to the same
Euler-Lagrange equations but to a different variational integrator.

"""
module LotkaVolterra2dGauge

    ϑ₁(t, q) = + log(q[2]) / q[1] / 2
    ϑ₂(t, q) = - log(q[1]) / q[2] / 2

    dϑ₁dx₁(t, q) = - log(q[2]) / q[1]^2 / 2
    dϑ₁dx₂(t, q) = + 1 / (q[1] * q[2]) / 2

    dϑ₂dx₁(t, q) = - 1 / (q[2] * q[1]) / 2
    dϑ₂dx₂(t, q) = + log(q[1]) / q[2]^2 / 2


    include("lotka_volterra_2d_common.jl")
    include("lotka_volterra_2d_equations.jl")

end
