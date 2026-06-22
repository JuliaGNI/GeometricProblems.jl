@doc raw"""
# Mathematical Pendulum

"""
module Pendulum

    using GeometricEquations

    export odeproblem, podeproblem, iodeproblem, idaeproblem

    export hamiltonian

    export plot_pendulum, plot_solution
    export labels_ode, labels_hamiltonian


    const t₀ = 0.0
    const Δt = 0.1
    const nt = 10

    timestep(::Type{T}=Float64, Δt=Δt) where {T} = T(Δt)
    timespan(::Type{T}=Float64, t₀=t₀, t₁=Δt * nt) where {T} = (T(t₀), T(t₁))


    const q₀ = [acos(0.4)]
    const p₀ = [0.0]
    const x₀ = [q₀[1], p₀[1]]
    const p₀_iode = [0.0, 0.0]


    hamiltonian(t, q::Number, p::Number) = p^2 / 2 + cos(q)
    hamiltonian(t, x::AbstractArray) = hamiltonian(t, x[1], x[2])
    hamiltonian(t, q::AbstractArray, p::AbstractArray) = hamiltonian(t, q[1], p[1])

    const labels_ode = (t = "t", q = "θ", p = "θ̇", h = "E")
    const labels_hamiltonian = (t = "t", q = "θ", p = "p", h = "H")

    function plot_pendulum end
    function plot_solution end


    function pendulum_ode_v(v, t, q, params)
        v[1] = q[2]
        v[2] = sin(q[1])
        nothing
    end

    function odeproblem(x₀::AbstractArray{DT}=x₀, ::Type{T}=DT; timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(x₀) == 2
        ODEProblem(pendulum_ode_v, timespan, timestep, T.(x₀))
    end

    odeproblem(::Type{T}, args...; kwargs...) where {T} = odeproblem(x₀, T, args...; kwargs...)


    function pendulum_pode_v(v, t, q, p, params)
        v[1] = p[1]
        nothing
    end

    function pendulum_pode_f(f, t, q, p, params)
        f[1] = sin(q[1])
        nothing
    end

    function podeproblem(q₀::AbstractArray{DT}=q₀, p₀::AbstractArray{DT}=p₀, ::Type{T}=DT; timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == 1
        PODEProblem(pendulum_pode_v, pendulum_pode_f, timespan, timestep, T.(q₀), T.(p₀))
    end

    podeproblem(::Type{T}, args...; kwargs...) where {T} = podeproblem(q₀, p₀, T, args...; kwargs...)


    function pendulum_iode_ϑ(p, t, q, v, params)
        p[1] = q[2]
        p[2] = 0
        nothing
    end

    function pendulum_iode_f(f, t, q, v, params)
        f[1] = sin(q[1])
        f[2] = v[1] - q[2]
        nothing
    end

    function pendulum_iode_g(g, t, q, v, λ, params)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function pendulum_iode_v(v, t, q, p, params)
        v[1] = q[2]
        v[2] = sin(q[1])
        nothing
    end

    function iodeproblem(q₀::AbstractArray{DT}=x₀, p₀::AbstractArray{DT}=p₀_iode, ::Type{T}=DT; timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == 2
        IODEProblem(pendulum_iode_ϑ, pendulum_iode_f,
            pendulum_iode_g, timespan, timestep, T.(q₀), T.(p₀);
            v̄=pendulum_iode_v)
    end

    iodeproblem(::Type{T}, args...; kwargs...) where {T} = iodeproblem(x₀, p₀_iode, T, args...; kwargs...)


    function pendulum_idae_u(u, t, q, v, p, λ, params)
        u[1] = λ[1]
        u[2] = λ[2]
        nothing
    end

    function pendulum_idae_g(g, t, q, v, p, λ, params)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function pendulum_idae_ϕ(ϕ, t, q, v, p, params)
        ϕ[1] = p[1] - q[2]
        ϕ[2] = p[2]
        nothing
    end

    function idaeproblem(q₀::AbstractArray{DT}=x₀, p₀::AbstractArray{DT}=p₀_iode, λ₀::AbstractArray{DT}=zero(q₀), ::Type{T}=DT; timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == length(λ₀) == 2
        IDAEProblem(pendulum_iode_ϑ, pendulum_iode_f,
            pendulum_idae_u, pendulum_idae_g, pendulum_idae_ϕ,
            timespan, timestep, T.(q₀), T.(p₀), T.(λ₀);
            v̄=pendulum_iode_v)
    end

    idaeproblem(::Type{T}, args...; kwargs...) where {T} = idaeproblem(x₀, p₀_iode, zero(x₀), T, args...; kwargs...)

end
