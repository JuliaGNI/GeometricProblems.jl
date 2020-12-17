@doc raw"""
# Lotka-Volterra model in 2D

```math
\begin{aligned}
L (q, \dot{q}) &= \bigg( q_2 + \frac{\log q_2}{q_1} \bigg) \, \dot{q_1} + q_1 \, \dot{q_2} - H(q) , \\
H(q) &= a_1 \, q_1 + a_2 \, q_2 + b_1 \, \log q_1 + b_2 \, \log q_2
\end{aligned}
```

"""
module LotkaVolterra2d

    export lotka_volterra_2d_dg_gauge

    ϑ₁(t, q) = q[2] + log(q[2]) / q[1]
    ϑ₂(t, q) = q[1]

    dϑ₁dx₁(t, q) = - log(q[2]) / q[1]^2
    dϑ₁dx₂(t, q) = 1 + 1 / (q[1] * q[2])

    dϑ₂dx₁(t, q) = one(eltype(q))
    dϑ₂dx₂(t, q) = zero(eltype(q))

    
    include("lotka_volterra_2d_common.jl")
    include("lotka_volterra_2d_equations.jl")


    function d²ϑ₁d₁d₁(t, q)
        + 2 * log(q[2]) / q[1]^3
    end

    function d²ϑ₁d₁d₂(t, q)
        - 1 / (q[1]^2 * q[2])
    end

    function d²ϑ₁d₂d₁(t, q)
        - 1 / (q[1]^2 * q[2])
    end

    function d²ϑ₁d₂d₂(t, q)
        - 1 / (q[1] * q[2]^2)
    end

    function d²ϑ₂d₁d₁(t, q)
        zero(eltype(q))
    end

    function d²ϑ₂d₁d₂(t, q)
        zero(eltype(q))
    end

    function d²ϑ₂d₂d₁(t, q)
        zero(eltype(q))
    end

    function d²ϑ₂d₂d₂(t, q)
        zero(eltype(q))
    end


    function g̅₁(t, q, v)
        d²ϑ₁d₁d₁(t,q) * q[1] * v[1] + d²ϑ₁d₂d₁(t,q) * q[1] * v[2] + d²ϑ₂d₁d₁(t,q) * q[2] * v[1] + d²ϑ₂d₂d₁(t,q) * q[2] * v[2]
    end

    function g̅₂(t, q, v)
        d²ϑ₁d₁d₂(t,q) * q[1] * v[1] + d²ϑ₁d₂d₂(t,q) * q[1] * v[2] + d²ϑ₂d₁d₂(t,q) * q[2] * v[1] + d²ϑ₂d₂d₂(t,q) * q[2] * v[2]
    end


    function lotka_volterra_2d_ϑ_κ(κ, t, q, v, Θ, params)
        Θ[1] = (1-κ) * ϑ₁(t,q) - κ * f₁(t,q,q)
        Θ[2] = (1-κ) * ϑ₂(t,q) - κ * f₂(t,q,q)
        nothing
    end

    function lotka_volterra_2d_f_κ(κ::Real, t::Real, q::Vector, v::Vector, f::Vector, params)
        f[1] = (1-κ) * f₁(t,q,v) - κ * (g₁(t,q,v) + g̅₁(t,q,v)) - dHd₁(t, q, params)
        f[2] = (1-κ) * f₂(t,q,v) - κ * (g₂(t,q,v) + g̅₂(t,q,v)) - dHd₂(t, q, params)
        nothing
    end

    function lotka_volterra_2d_g_κ(κ::Real, t::Real, q::Vector, v::Vector, g::Vector, params)
        g[1] = (1-κ) * f₁(t,q,v) - κ * (g₁(t,q,v) + g̅₁(t,q,v))
        g[2] = (1-κ) * f₂(t,q,v) - κ * (g₂(t,q,v) + g̅₂(t,q,v))
        nothing
    end

    # function lotka_volterra_2d_g(κ::Real, t::Real, q::Vector, v::Vector, g::Vector)
    #     g[1] = (1-κ) * g₁(t,q,v) - κ * g̅₁(t,q,v) - κ * f₁(t,q,v)
    #     g[2] = (1-κ) * g₂(t,q,v) - κ * g̅₂(t,q,v) - κ * f₂(t,q,v)
    #     nothing
    # end


    function lotka_volterra_2d_dg_gauge(q₀=q₀, p₀=ϑ(0, q₀), κ=0; params=parameters)
        lotka_volterra_2d_ϑ = (t, q, v, p, params) -> lotka_volterra_2d_ϑ_κ(κ, t, q, v, p, params)
        lotka_volterra_2d_f = (t, q, v, f, params) -> lotka_volterra_2d_f_κ(κ, t, q, v, f, params)
        lotka_volterra_2d_g = (t, q, λ, g, params) -> lotka_volterra_2d_g_κ(κ, t, q, λ, g, params)

        IODE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
             lotka_volterra_2d_g, q₀, p₀;
             parameters=params, v̄=lotka_volterra_2d_v)
    end

end
