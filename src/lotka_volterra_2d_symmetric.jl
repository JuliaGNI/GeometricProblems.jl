@doc raw"""
# Lotka-Volterra model in 2D with symmetric Lagrangian

```math
\begin{aligned}
L (q, \dot{q}) &= \frac{1}{2} \frac{\log q_2}{q_1} \, \dot{q_1} - \frac{1}{2} \frac{\log q_1}{q_2} \, \dot{q_2} - H(q) , \\
H(q) &= a_1 \, q_1 + a_2 \, q_2 + b_1 \, \log q_1 + b_2 \, \log q_2
\end{aligned}
```

This Lagrangian is a slight generalization of Equation (5) in José Fernández-Núñez,
Lagrangian Structure of the Two-Dimensional Lotka-Volterra System, International
Journal of Theoretical Physics, Vol. 37, No. 9, pp. 2457-2462, 1998.

"""
module LotkaVolterra2dSymmetric

    ϑ₁(t, q) = + log(q[2]) / q[1] / 2
    ϑ₂(t, q) = - log(q[1]) / q[2] / 2

    dϑ₁dx₁(t, q) = - log(q[2]) / q[1]^2 / 2
    dϑ₁dx₂(t, q) = + 1 / (q[1] * q[2]) / 2

    dϑ₂dx₁(t, q) = - 1 / (q[2] * q[1]) / 2
    dϑ₂dx₂(t, q) = + log(q[1]) / q[2]^2 / 2


    include("lotka_volterra_2d_common.jl")
    include("lotka_volterra_2d_equations.jl")

end
