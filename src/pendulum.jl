@doc raw"""
# Mathematical Pendulum

The mathematical pendulum with Hamiltonian
```math
H(q, p) = \frac{p^2}{2 m l^2} + m g l \cos(q) ,
```
where `q` is the angle and `p` the conjugate momentum.

System parameters:
* `l`: length of the pendulum
* `m`: mass of the pendulum
* `g`: gravitational acceleration

The default parameters `l = m = g = 1` reduce this to the normalised
Hamiltonian ``H(q, p) = p^2 / 2 + \cos(q)``.
"""
module Pendulum

    using GeometricEquations
    using Parameters

    export odeproblem, podeproblem, hodeproblem, iodeproblem, idaeproblem

    export hodeensemble

    export hamiltonian

    export default_parameters

    export plot_pendulum, plot_solution
    export labels_ode, labels_hamiltonian


    const t₀ = 0.0
    const Δt = 0.1
    const nt = 10

    timestep(::Type{T}=Float64, Δt=Δt) where {T} = T(Δt)
    timespan(::Type{T}=Float64, t₀=t₀, t₁=Δt * nt) where {T} = (T(t₀), T(t₁))


    # Physical parameters of the mathematical pendulum.
    # The defaults l = m = g = 1 reproduce the normalised Hamiltonian H = p²/2 + cos(q).
    const l = 1.0   # length of the pendulum
    const m = 1.0   # mass of the pendulum
    const g = 1.0   # gravitational acceleration

    default_parameters(::Type{T}=Float64) where {T} = (l = T(l), m = T(m), g = T(g))


    const q₀ = [acos(0.4)]
    const p₀ = [0.0]
    const x₀ = [q₀[1], p₀[1]]
    const p₀_iode = [0.0, 0.0]

    const qmin = [0.0]
    const qmax = [2π]
    const pmin = [-2.0]
    const pmax = [+2.0]
    const qsamples = [10]
    const psamples = [10]


    function _pode_samples(qmin, qmax, pmin, pmax, qsamples, psamples)
        qs = [range(qmin[i], qmax[i]; length=qsamples[i]) for i in eachindex(qmin, qmax, qsamples)]
        ps = [range(pmin[i], pmax[i]; length=psamples[i]) for i in eachindex(pmin, pmax, psamples)]

        qsamples = vec(collect.(collect(Base.Iterators.product(qs...))))
        psamples = vec(collect.(collect(Base.Iterators.product(ps...))))
        zsamples = Base.Iterators.product(qsamples, psamples)

        return (
            q=vec([zs[1] for zs in zsamples]),
            p=vec([zs[2] for zs in zsamples]),
        )
    end


    function hamiltonian(t, q::Number, p::Number, params)
        @unpack l, m, g = params
        p^2 / (2 * m * l^2) + m * g * l * cos(q)
    end

    hamiltonian(t, x::AbstractArray, params) = hamiltonian(t, x[1], x[2], params)
    hamiltonian(t, q::AbstractArray, p::AbstractArray, params) = hamiltonian(t, q[1], p[1], params)

    # parameter-free convenience methods using the default parameters
    hamiltonian(t, q::Number, p::Number) = hamiltonian(t, q, p, default_parameters())
    hamiltonian(t, x::AbstractArray) = hamiltonian(t, x[1], x[2], default_parameters())
    hamiltonian(t, q::AbstractArray, p::AbstractArray) = hamiltonian(t, q[1], p[1], default_parameters())

    const labels_ode = (t = "t", q = "θ", p = "θ̇", h = "E")
    const labels_hamiltonian = (t = "t", q = "θ", p = "p", h = "H")

    function plot_pendulum end
    function plot_solution end


    function pendulum_ode_v(v, t, q, params)
        @unpack l, g = params
        v[1] = q[2]
        v[2] = (g / l) * sin(q[1])
        nothing
    end

    function odeproblem(x₀::AbstractArray{DT}=x₀, ::Type{T}=DT; parameters=default_parameters(T), timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(x₀) == 2
        ODEProblem(pendulum_ode_v, timespan, timestep, T.(x₀); parameters=parameters)
    end

    odeproblem(::Type{T}, args...; kwargs...) where {T} = odeproblem(x₀, T, args...; kwargs...)


    function pendulum_pode_v(v, t, q, p, params)
        @unpack l, m = params
        v[1] = p[1] / (m * l^2)
        nothing
    end

    function pendulum_pode_f(f, t, q, p, params)
        @unpack l, m, g = params
        f[1] = m * g * l * sin(q[1])
        nothing
    end

    function podeproblem(q₀::AbstractArray{DT}=q₀, p₀::AbstractArray{DT}=p₀, ::Type{T}=DT; parameters=default_parameters(T), timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == 1
        PODEProblem(pendulum_pode_v, pendulum_pode_f, timespan, timestep, T.(q₀), T.(p₀); parameters=parameters)
    end

    podeproblem(::Type{T}, args...; kwargs...) where {T} = podeproblem(q₀, p₀, T, args...; kwargs...)


    function hodeproblem(q₀::AbstractArray{DT}=q₀, p₀::AbstractArray{DT}=p₀, ::Type{T}=DT; parameters=default_parameters(T), timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == 1
        HODEProblem(pendulum_pode_v, pendulum_pode_f, hamiltonian, timespan, timestep, T.(q₀), T.(p₀); parameters=parameters)
    end

    hodeproblem(::Type{T}, args...; kwargs...) where {T} = hodeproblem(q₀, p₀, T, args...; kwargs...)

    # Ensemble over a grid of initial conditions (and optionally a vector of parameter sets).
    function hodeensemble(
        qmin::AbstractArray{DT}=qmin,
        qmax::AbstractArray{DT}=qmax,
        pmin::AbstractArray{DT}=pmin,
        pmax::AbstractArray{DT}=pmax,
        qsamples=qsamples,
        psamples=psamples,
        ::Type{T}=DT;
        parameters=default_parameters(T),
        timespan=timespan(T),
        timestep=timestep(T)) where {DT,T}
        samples = _pode_samples(T.(qmin), T.(qmax), T.(pmin), T.(pmax), qsamples, psamples)
        HODEEnsemble(pendulum_pode_v, pendulum_pode_f, hamiltonian, timespan, timestep, samples.q, samples.p; parameters=parameters)
    end

    # Ensemble over a vector of parameter sets sharing a single initial condition (q₀, p₀).
    function hodeensemble(q₀::AbstractArray{DT}, p₀::AbstractArray{DT}, parameters::AbstractVector{<:NamedTuple}, ::Type{T}=DT; timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == 1
        HODEEnsemble(pendulum_pode_v, pendulum_pode_f, hamiltonian, timespan, timestep, T.(q₀), T.(p₀); parameters=parameters)
    end


    function pendulum_iode_ϑ(p, t, q, v, params)
        @unpack l, m = params
        p[1] = m * l^2 * q[2]
        p[2] = 0
        nothing
    end

    function pendulum_iode_f(f, t, q, v, params)
        @unpack l, m, g = params
        f[1] = m * g * l * sin(q[1])
        f[2] = m * l^2 * (v[1] - q[2])
        nothing
    end

    function pendulum_iode_g(g, t, q, v, λ, params)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function pendulum_iode_v(v, t, q, p, params)
        @unpack l, g = params
        v[1] = q[2]
        v[2] = (g / l) * sin(q[1])
        nothing
    end

    function iodeproblem(q₀::AbstractArray{DT}=x₀, p₀::AbstractArray{DT}=p₀_iode, ::Type{T}=DT; parameters=default_parameters(T), timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == 2
        IODEProblem(pendulum_iode_ϑ, pendulum_iode_f,
            pendulum_iode_g, timespan, timestep, T.(q₀), T.(p₀);
            parameters=parameters, v̄=pendulum_iode_v)
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
        @unpack l, m = params
        ϕ[1] = p[1] - m * l^2 * q[2]
        ϕ[2] = p[2]
        nothing
    end

    function idaeproblem(q₀::AbstractArray{DT}=x₀, p₀::AbstractArray{DT}=p₀_iode, λ₀::AbstractArray{DT}=zero(q₀), ::Type{T}=DT; parameters=default_parameters(T), timespan=timespan(T), timestep=timestep(T)) where {DT,T}
        @assert length(q₀) == length(p₀) == length(λ₀) == 2
        IDAEProblem(pendulum_iode_ϑ, pendulum_iode_f,
            pendulum_idae_u, pendulum_idae_g, pendulum_idae_ϕ,
            timespan, timestep, T.(q₀), T.(p₀), T.(λ₀);
            parameters=parameters, v̄=pendulum_iode_v)
    end

    idaeproblem(::Type{T}, args...; kwargs...) where {T} = idaeproblem(x₀, p₀_iode, zero(x₀), T, args...; kwargs...)

end
