@doc raw"""

"""
module LotkaVolterra4d

    using GeometricEquations
    using Parameters

    export hamiltonian, ϑ, ϑ₁, ϑ₂, ω

    export lotka_volterra_4d_ode, 
           lotka_volterra_4d_pode, lotka_volterra_4d_pdae,
           lotka_volterra_4d_iode, lotka_volterra_4d_idae,
           lotka_volterra_4d_lode, lotka_volterra_4d_ldae,
           lotka_volterra_4d_dg


    const Δt = 0.01
    const nt = 1000
    const tspan = (0.0, Δt*nt)

    const q₀ = [2.0, 1.0, 1.0, 1.0]

    const default_parameters = (a₁=1.0, a₂=1.0, a₃=1.0, a₄=1.0, b₁=-1.0, b₂=-2.0, b₃=-1.0, b₄=-1.0)
    const reference_solution = [0.5988695239096916, 2.068567531039674, 0.2804351458645534, 1.258449091830993]

    # const q₀ = [2.0, 1.0, 2.0, 1.0]
    # const p  = (a₁=1.0, a₂=1.0, a₃=1.0, a₄=1.0, b₁=-1.0, b₂=-2.0, b₃=-1.0, b₄=-2.0)

    # const q₀ = [2.0, 1.0, 2.0, 1.0]
    # const p  = (a₁=1.0, a₂=1.0, a₃=1.0, a₄=1.0, b₁=-1.0, b₂=-4.0, b₃=-2.0, b₄=-3.0)


    ϑ₁(t, q) = 0.5 * ( + log(q[2]) - log(q[3]) + log(q[4]) ) / q[1]
    ϑ₂(t, q) = 0.5 * ( - log(q[1]) + log(q[3]) - log(q[4]) ) / q[2]
    ϑ₃(t, q) = 0.5 * ( + log(q[1]) - log(q[2]) + log(q[4]) ) / q[3]
    ϑ₄(t, q) = 0.5 * ( - log(q[1]) + log(q[2]) - log(q[3]) ) / q[4]

    # ϑ₁(t, q) = ( log(q[2]) + log(q[4]) ) / q[1]
    # ϑ₂(t, q) = ( log(q[3]) ) / q[2]
    # ϑ₃(t, q) = ( log(q[1]) + log(q[4]) ) / q[3]
    # ϑ₄(t, q) = ( log(q[2]) ) / q[4]

    # ϑ₁(t, q) = ( + log(q[2]) - log(q[3]) + log(q[4]) ) / q[1] / 2 + q[2] + q[3] + q[4]
    # ϑ₂(t, q) = ( - log(q[1]) + log(q[3]) - log(q[4]) ) / q[2] / 2 + q[1] + q[3] + q[4]
    # ϑ₃(t, q) = ( + log(q[1]) - log(q[2]) + log(q[4]) ) / q[3] / 2 + q[1] + q[2] + q[4]
    # ϑ₄(t, q) = ( - log(q[1]) + log(q[2]) - log(q[3]) ) / q[4] / 2 + q[1] + q[2] + q[3]

    # ϑ₁(t, q) = ( + log(q[2]) - log(q[3]) + log(q[4]) ) / q[1]
    # ϑ₂(t, q) = ( + log(q[3]) - log(q[4]) ) / q[2]
    # ϑ₃(t, q) = log(q[4]) / q[3]
    # ϑ₄(t, q) = zero(eltype(q))


    function v₁(t, q, params)
        @unpack a₁, a₂, a₃, a₄, b₁, b₂, b₃, b₄ = params
        q[1] * (+ a₂*q[2] + a₃*q[3] + a₄*q[4] + b₂ + b₃ + b₄)
    end

    function v₂(t, q, params)
        @unpack a₁, a₂, a₃, a₄, b₁, b₂, b₃, b₄ = params
        q[2] * (- a₁*q[1] + a₃*q[3] + a₄*q[4] - b₁ + b₃ + b₄)
    end

    function v₃(t, q, params)
        @unpack a₁, a₂, a₃, a₄, b₁, b₂, b₃, b₄ = params
        q[3] * (- a₁*q[1] - a₂*q[2] + a₄*q[4] - b₁ - b₂ + b₄)
    end

    function v₄(t, q, params)
        @unpack a₁, a₂, a₃, a₄, b₁, b₂, b₃, b₄ = params
        q[4] * (- a₁*q[1] - a₂*q[2] - a₃*q[3] - b₁ - b₂ - b₃)
    end


    dϑ₁dx₁(t, q) = ( - log(q[2]) + log(q[3]) - log(q[4]) ) / q[1]^2 / 2
    dϑ₁dx₂(t, q) = + 1 / (q[1] * q[2]) / 2
    dϑ₁dx₃(t, q) = - 1 / (q[1] * q[3]) / 2
    dϑ₁dx₄(t, q) = + 1 / (q[1] * q[4]) / 2

    dϑ₂dx₁(t, q) = - 1 / (q[2] * q[1]) / 2
    dϑ₂dx₂(t, q) = ( + log(q[1]) - log(q[3]) + log(q[4]) ) / q[2]^2 / 2
    dϑ₂dx₃(t, q) = + 1 / (q[2] * q[3]) / 2
    dϑ₂dx₄(t, q) = - 1 / (q[2] * q[4]) / 2

    dϑ₃dx₁(t, q) = + 1 / (q[3] * q[1]) / 2
    dϑ₃dx₂(t, q) = - 1 / (q[3] * q[2]) / 2
    dϑ₃dx₃(t, q) = ( - log(q[1]) + log(q[2]) - log(q[4]) ) / q[3]^2 / 2
    dϑ₃dx₄(t, q) = + 1 / (q[3] * q[4]) / 2

    dϑ₄dx₁(t, q) = - 1 / (q[4] * q[1]) / 2
    dϑ₄dx₂(t, q) = + 1 / (q[4] * q[2]) / 2
    dϑ₄dx₃(t, q) = - 1 / (q[4] * q[3]) / 2
    dϑ₄dx₄(t, q) = ( + log(q[1]) - log(q[2]) + log(q[3]) ) / q[4]^2 / 2


    # dϑ₁dx₁(t, q) = - ( log(q[2]) + log(q[4]) ) / q[1]^2
    # dϑ₁dx₂(t, q) = 1 / (q[1] * q[2])
    # dϑ₁dx₃(t, q) = zero(eltype(q))
    # dϑ₁dx₄(t, q) = 1 / (q[1] * q[4])

    # dϑ₂dx₁(t, q) = zero(eltype(q))
    # dϑ₂dx₂(t, q) = - ( log(q[3]) ) / q[2]^2
    # dϑ₂dx₃(t, q) = 1 / (q[2] * q[3])
    # dϑ₂dx₄(t, q) = zero(eltype(q))

    # dϑ₃dx₁(t, q) = + 1 / (q[3] * q[1])
    # dϑ₃dx₂(t, q) = zero(eltype(q))
    # dϑ₃dx₃(t, q) = - ( log(q[1]) + log(q[4]) ) / q[3]^2
    # dϑ₃dx₄(t, q) = + 1 / (q[3] * q[4])

    # dϑ₄dx₁(t, q) = zero(eltype(q))
    # dϑ₄dx₂(t, q) = 1 / (q[4] * q[2])
    # dϑ₄dx₃(t, q) = zero(eltype(q))
    # dϑ₄dx₄(t, q) = - ( log(q[2]) ) / q[4]^2


    # dϑ₁dx₁(t, q) = ( - log(q[2]) + log(q[3]) - log(q[4]) ) / q[1]^2 / 2
    # dϑ₁dx₂(t, q) = 1 + 1 / (q[1] * q[2]) / 2
    # dϑ₁dx₃(t, q) = 1 - 1 / (q[1] * q[3]) / 2
    # dϑ₁dx₄(t, q) = 1 + 1 / (q[1] * q[4]) / 2

    # dϑ₂dx₁(t, q) = 1 - 1 / (q[2] * q[1]) / 2
    # dϑ₂dx₂(t, q) = ( + log(q[1]) - log(q[3]) + log(q[4]) ) / q[2]^2 / 2
    # dϑ₂dx₃(t, q) = 1 + 1 / (q[2] * q[3]) / 2
    # dϑ₂dx₄(t, q) = 1 - 1 / (q[2] * q[4]) / 2

    # dϑ₃dx₁(t, q) = 1 + 1 / (q[3] * q[1]) / 2
    # dϑ₃dx₂(t, q) = 1 - 1 / (q[3] * q[2]) / 2
    # dϑ₃dx₃(t, q) = ( - log(q[1]) + log(q[2]) - log(q[4]) ) / q[3]^2 / 2
    # dϑ₃dx₄(t, q) = 1 + 1 / (q[3] * q[4]) / 2

    # dϑ₄dx₁(t, q) = 1 - 1 / (q[4] * q[1]) / 2
    # dϑ₄dx₂(t, q) = 1 + 1 / (q[4] * q[2]) / 2
    # dϑ₄dx₃(t, q) = 1 - 1 / (q[4] * q[3]) / 2
    # dϑ₄dx₄(t, q) = ( + log(q[1]) - log(q[2]) + log(q[3]) ) / q[4]^2 / 2


    # dϑ₁dx₁(t, q) = ( - log(q[2]) + log(q[3]) - log(q[4]) ) / q[1]^2
    # dϑ₁dx₂(t, q) = + 1 / (q[1] * q[2])
    # dϑ₁dx₃(t, q) = - 1 / (q[1] * q[3])
    # dϑ₁dx₄(t, q) = + 1 / (q[1] * q[4])

    # dϑ₂dx₁(t, q) = zero(eltype(q))
    # dϑ₂dx₂(t, q) = ( - log(q[3]) + log(q[4]) ) / q[2]^2
    # dϑ₂dx₃(t, q) = + 1 / (q[2] * q[3])
    # dϑ₂dx₄(t, q) = - 1 / (q[2] * q[4])

    # dϑ₃dx₁(t, q) = zero(eltype(q))
    # dϑ₃dx₂(t, q) = zero(eltype(q))
    # dϑ₃dx₃(t, q) = ( - log(q[4]) ) / q[3]^2
    # dϑ₃dx₄(t, q) = + 1 / (q[3] * q[4])

    # dϑ₄dx₁(t, q) = zero(eltype(q))
    # dϑ₄dx₂(t, q) = zero(eltype(q))
    # dϑ₄dx₃(t, q) = zero(eltype(q))
    # dϑ₄dx₄(t, q) = zero(eltype(q))


    function ϑ(Θ::AbstractVector, t::Number, q::AbstractVector)
        Θ[1] = ϑ₁(t,q)
        Θ[2] = ϑ₂(t,q)
        Θ[3] = ϑ₃(t,q)
        Θ[4] = ϑ₄(t,q)
        nothing
    end

    function ϑ(t::Number, q::AbstractVector)
        Θ = zero(q)
        ϑ(Θ, t, q)
        return Θ
    end

    function ϑ(t::Number, q::AbstractVector, k::Int)
        if k == 1
            ϑ₁(t, q)
        elseif k == 2
            ϑ₂(t, q)
        elseif k == 3
            ϑ₃(t, q)
        elseif k == 4
            ϑ₄(t, q)
        else
            throw(BoundsError(ϑ,k))
        end
    end


    function ω(Ω, t, q)
        Ω[1,1] = 0
        Ω[1,2] = dϑ₁dx₂(t,q) - dϑ₂dx₁(t,q)
        Ω[1,3] = dϑ₁dx₃(t,q) - dϑ₃dx₁(t,q)
        Ω[1,4] = dϑ₁dx₄(t,q) - dϑ₄dx₁(t,q)

        Ω[2,1] = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)
        Ω[2,2] = 0
        Ω[2,3] = dϑ₂dx₃(t,q) - dϑ₃dx₂(t,q)
        Ω[2,4] = dϑ₂dx₄(t,q) - dϑ₄dx₂(t,q)

        Ω[3,1] = dϑ₃dx₁(t,q) - dϑ₁dx₃(t,q)
        Ω[3,2] = dϑ₃dx₂(t,q) - dϑ₂dx₃(t,q)
        Ω[3,3] = 0
        Ω[3,4] = dϑ₃dx₄(t,q) - dϑ₄dx₃(t,q)

        Ω[4,1] = dϑ₄dx₁(t,q) - dϑ₁dx₄(t,q)
        Ω[4,2] = dϑ₄dx₂(t,q) - dϑ₂dx₄(t,q)
        Ω[4,3] = dϑ₄dx₃(t,q) - dϑ₃dx₄(t,q)
        Ω[4,4] = 0

        nothing
    end


    f₁(t, q, v) = dϑ₁dx₁(t,q) * v[1] + dϑ₂dx₁(t,q) * v[2] + dϑ₃dx₁(t,q) * v[3] + dϑ₄dx₁(t,q) * v[4]
    f₂(t, q, v) = dϑ₁dx₂(t,q) * v[1] + dϑ₂dx₂(t,q) * v[2] + dϑ₃dx₂(t,q) * v[3] + dϑ₄dx₂(t,q) * v[4]
    f₃(t, q, v) = dϑ₁dx₃(t,q) * v[1] + dϑ₂dx₃(t,q) * v[2] + dϑ₃dx₃(t,q) * v[3] + dϑ₄dx₃(t,q) * v[4]
    f₄(t, q, v) = dϑ₁dx₄(t,q) * v[1] + dϑ₂dx₄(t,q) * v[2] + dϑ₃dx₄(t,q) * v[3] + dϑ₄dx₄(t,q) * v[4]

    g₁(t, q, v) = dϑ₁dx₁(t,q) * v[1] + dϑ₁dx₂(t,q) * v[2] + dϑ₁dx₃(t,q) * v[3] + dϑ₁dx₄(t,q) * v[4]
    g₂(t, q, v) = dϑ₂dx₁(t,q) * v[1] + dϑ₂dx₂(t,q) * v[2] + dϑ₂dx₃(t,q) * v[3] + dϑ₂dx₄(t,q) * v[4]
    g₃(t, q, v) = dϑ₃dx₁(t,q) * v[1] + dϑ₃dx₂(t,q) * v[2] + dϑ₃dx₃(t,q) * v[3] + dϑ₃dx₄(t,q) * v[4]
    g₄(t, q, v) = dϑ₄dx₁(t,q) * v[1] + dϑ₄dx₂(t,q) * v[2] + dϑ₄dx₃(t,q) * v[3] + dϑ₄dx₄(t,q) * v[4]


    function hamiltonian(t, q, params)
        @unpack a₁, a₂, a₃, a₄, b₁, b₂, b₃, b₄ = params
        a = [a₁, a₂, a₃, a₄]
        b = [b₁, b₂, b₃, b₄]
        sum(a .* q) + sum(b .* log.(q))
    end

    hamiltonian_iode(t, q, v, params) = hamiltonian(t, q, params)

    hamiltonian_pode(t, q, p, params) = hamiltonian(t, q, params)


    function dHd₁(t, q, params)
        @unpack a₁, b₁ = params
        a₁ + b₁ / q[1]
    end

    function dHd₂(t, q, params)
        @unpack a₂, b₂ = params
        a₂ + b₂ / q[2]
    end

    function dHd₃(t, q, params)
        @unpack a₃, b₃ = params
        a₃ + b₃ / q[3]
    end

    function dHd₄(t, q, params)
        @unpack a₄, b₄ = params
        a₄ + b₄ / q[4]
    end

    function lotka_volterra_4d_dH(dH, t, q, params)
        dH[1] = dHd₁(t, q, params)
        dH[2] = dHd₂(t, q, params)
        dH[3] = dHd₃(t, q, params)
        dH[4] = dHd₄(t, q, params)
        nothing
    end


    lotka_volterra_4d_ϑ(Θ, t, q, params) = ϑ(Θ, t, q)
    lotka_volterra_4d_ϑ(Θ, t, q, v, params) = ϑ(Θ, t, q)
    lotka_volterra_4d_ω(Ω, t, q, params) = ω(Ω, t, q)


    function lotka_volterra_4d_v(v, t, q, params)
        v[1] = v₁(t, q, params)
        v[2] = v₂(t, q, params)
        v[3] = v₃(t, q, params)
        v[4] = v₄(t, q, params)
        nothing
    end

    function lotka_volterra_4d_v(v, t, q, p, params)
        lotka_volterra_4d_v(v, t, q, params)
    end

    function lotka_volterra_4d_v_ham(v, t, q, p, params)
        v .= 0
        nothing
    end

    function lotka_volterra_4d_f(f::AbstractVector, t::Real, q::AbstractVector, v::AbstractVector, params)
        f[1] = f₁(t,q,v) - dHd₁(t, q, params)
        f[2] = f₂(t,q,v) - dHd₂(t, q, params)
        f[3] = f₃(t,q,v) - dHd₃(t, q, params)
        f[4] = f₄(t,q,v) - dHd₄(t, q, params)
        nothing
    end

    function lotka_volterra_4d_f_ham(f::AbstractVector, t::Real, q::AbstractVector, v::AbstractVector, params)
        f[1] = - dHd₁(t, q, params)
        f[2] = - dHd₂(t, q, params)
        f[3] = - dHd₃(t, q, params)
        f[4] = - dHd₄(t, q, params)
    end

    function lotka_volterra_4d_g(g::AbstractVector, t::Real, q::AbstractVector, v::AbstractVector, params)
        g[1] = f₁(t,q,v)
        g[2] = f₂(t,q,v)
        g[3] = f₃(t,q,v)
        g[4] = f₄(t,q,v)
        nothing
    end

    lotka_volterra_4d_g(g, t, q, p, λ, params) = lotka_volterra_4d_g(g, t, q, λ, params)
    lotka_volterra_4d_g(g, t, q, v, p, λ, params) = lotka_volterra_4d_g(g, t, q, p, λ, params)

    function lotka_volterra_4d_g̅(g::AbstractVector, t::Real, q::AbstractVector, v::AbstractVector, params)
        g[1] = g₁(t,q,v)
        g[2] = g₂(t,q,v)
        g[3] = g₃(t,q,v)
        g[4] = g₄(t,q,v)
        nothing
    end

    function lotka_volterra_4d_g̅(g::AbstractVector, t::Real, q::AbstractVector, p::AbstractVector, v::AbstractVector, params)
        lotka_volterra_4d_g̅(t, q, v, g, params)
    end

    function lotka_volterra_4d_u(u, t, q, λ, params)
        u .= λ
        nothing
    end

    lotka_volterra_4d_u(u, t, q, p, λ, params) = lotka_volterra_4d_u(u, t, q, λ, params)
    lotka_volterra_4d_u(u, t, q, v, p, λ, params) = lotka_volterra_4d_u(u, t, q, p, λ, params)


    function lotka_volterra_4d_ϕ(ϕ, t, q, p, params)
        ϕ[1] = p[1] - ϑ₁(t,q)
        ϕ[2] = p[2] - ϑ₂(t,q)
        ϕ[3] = p[3] - ϑ₃(t,q)
        ϕ[4] = p[4] - ϑ₄(t,q)
        nothing
    end

    lotka_volterra_4d_ϕ(ϕ, t, q, v, p, params) = lotka_volterra_4d_ϕ(ϕ, t, q, p, params)

    function lotka_volterra_4d_ψ(ψ, t, q, p, q̇, ṗ, params)
        ψ[1] = f[1] - g₁(t,q,q̇)
        ψ[2] = f[2] - g₂(t,q,q̇)
        ψ[3] = f[3] - g₃(t,q,q̇)
        ψ[4] = f[4] - g₄(t,q,q̇)
        nothing
    end

    lotka_volterra_4d_ψ(ψ, t, q, v, p, q̇, ṗ, params) = lotka_volterra_4d_ψ(ψ, t, q, p, q̇, ṗ, params)


    function lotka_volterra_4d_ode(q₀=q₀; tspan=tspan, tstep=Δt, parameters=default_parameters)
        ODEProblem(lotka_volterra_4d_v, tspan, tstep, q₀; parameters=parameters, invariants=(h=hamiltonian,))
    end


    function lotka_volterra_4d_pode(q₀=q₀, p₀=nothing; tspan=tspan, tstep=Δt, parameters=default_parameters)
        p₀ === nothing ? p₀ = ϑ(tspan[begin], q₀) : nothing
        PODEProblem(lotka_volterra_4d_v, lotka_volterra_4d_f,
                    tspan, tstep, q₀, p₀; parameters=parameters, invariants=(h=hamiltonian_pode,))
    end

    function lotka_volterra_4d_iode(q₀=q₀, p₀=nothing; tspan=tspan, tstep=Δt, parameters=default_parameters)
        p₀ === nothing ? p₀ = ϑ(tspan[begin], q₀) : nothing
        IODEProblem(lotka_volterra_4d_ϑ, lotka_volterra_4d_f,
                    lotka_volterra_4d_g, tspan, tstep, q₀, p₀;
                    parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_4d_v)
    end

    function lotka_volterra_4d_lode(q₀=q₀, p₀=nothing; tspan=tspan, tstep=Δt, parameters=default_parameters)
        p₀ === nothing ? p₀ = ϑ(tspan[begin], q₀) : nothing
        LODEProblem(lotka_volterra_4d_ϑ, lotka_volterra_4d_f,
                    lotka_volterra_4d_g, tspan, tstep, q₀, p₀;
                    parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_4d_v,
                    Ω=lotka_volterra_4d_ω, ∇H=lotka_volterra_4d_dH)
    end

    function lotka_volterra_4d_idae(q₀=q₀, p₀=nothing, λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
        p₀ === nothing ? p₀ = ϑ(tspan[begin], q₀) : nothing
        IDAEProblem(lotka_volterra_4d_ϑ, lotka_volterra_4d_f,
                    lotka_volterra_4d_u, lotka_volterra_4d_g,
                    lotka_volterra_4d_ϕ, tspan, tstep, q₀, p₀, λ₀;
                    parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_4d_v)
    end

    function lotka_volterra_4d_pdae(q₀=q₀, p₀=nothing, λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
        p₀ === nothing ? p₀ = ϑ(tspan[begin], q₀) : nothing
        PDAEProblem(lotka_volterra_4d_v_ham, lotka_volterra_4d_f_ham,
                    lotka_volterra_4d_u, lotka_volterra_4d_g,
                    lotka_volterra_4d_ϕ, tspan, tstep, q₀, p₀, λ₀;
                    v̄=lotka_volterra_4d_v, f̄=lotka_volterra_4d_f,
                    parameters=parameters, invariants=(h=hamiltonian_pode,))
    end

    function lotka_volterra_4d_ldae(q₀=q₀, p₀=nothing, λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
        p₀ === nothing ? p₀ = ϑ(tspan[begin], q₀) : nothing
        LDAEProblem(lotka_volterra_4d_ϑ, lotka_volterra_4d_f_ham,
                    lotka_volterra_4d_g, lotka_volterra_4d_g̅,
                    lotka_volterra_4d_ϕ, lotka_volterra_4d_ψ,
                    tspan, tstep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                    v̄=lotka_volterra_4d_v, f̄=lotka_volterra_4d_f,)
    end

    function lotka_volterra_4d_dg(q₀=q₀, p₀=nothing; tspan=tspan, tstep=Δt, parameters=default_parameters)
        p₀ === nothing ? p₀ = ϑ(tspan[begin], q₀) : nothing
        IODEProblem(lotka_volterra_4d_ϑ, lotka_volterra_4d_f,
                    lotka_volterra_4d_g, tspan, tstep, q₀, p₀;
                    parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_4d_v)
    end

end
