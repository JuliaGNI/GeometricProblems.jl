@doc raw"""
The Toda lattice is a model for a one-dimensional crystal named after its discoverer Morikazu Toda [toda1967vibration](@cite).

It is a prime example of a non-trivial completely integrable system.

The only system parameters are the *number of points* ``N`` in the periodic lattice and ``\alpha`` which adjusts the strength of the interactions in the lattice.
"""
module TodaLattice 

    using EulerLagrange
    using LinearAlgebra 
    using Parameters 
    using GeometricEquations: HODEEnsemble, HODEProblem, LODEEnsemble, LODEProblem

    export hodeproblem, hodeensemble, hamiltonian, hamiltonian_system
    export lodeproblem, lodeensemble, lagrangian, lagrangian_system

    include("bump_initial_condition.jl")

    const default_parameters = (
        α = .64, 
    )

    function potential(q, params, N)
        params.α * sum(exp(q[n] - q[n % N + 1]) for n in 1:N)
    end

    hamiltonian(t, q, p, params, N) = p ⋅ p / 2 + potential(q, params, N)
    lagrangian(t, q, q̇, params, N) = q̇ ⋅ q̇ / 2 - potential(q, params, N)

    const tstep = .1
    const tspan = (0.0, 120.0)

    # parameter for the default initial conditions
    const Ñ = 200
    const μ = .3

    const q₀ = compute_initial_q(μ, Ñ)
    const p₀ = zero(q₀)
    const Ω = compute_domain(Ñ, typeof(μ))


    function hamiltonian_system(N, parameters)
        t, q, p = hamiltonian_variables(N)
        sparams = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams, N), t, q, p, sparams; simplify = N ≤ 10)
    end

    function lagrangian_system(N, parameters)
        t, x, v = lagrangian_variables(N)
        sparams = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams, N), t, x, v, sparams; simplify = N ≤ 10)
    end

    _parameters(p::NamedTuple) = p
    _parameters(p::AbstractVector) = p[begin]

    """
    Hamiltonian problem for the Toda lattice.
    """
    function hodeproblem(N::Int = Ñ, q₀ = compute_initial_q(μ, N), p₀ = zero(q₀); tspan = tspan, tstep = tstep, parameters = default_parameters)
        HODEProblem(hamiltonian_system(N, parameters), tspan, tstep, q₀, p₀; parameters = parameters)
    end

    function hodeproblem(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        hodeproblem(length(q₀), q₀, p₀; kwargs...)
    end

    function hodeensemble(N::Int = Ñ, q₀ = compute_initial_q(μ, N), p₀ = zero(q₀); tspan = tspan, tstep = tstep, parameters = default_parameters)
        eqs = functions(hamiltonian_system(N, _parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    function hodeensemble(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        hodeensemble(length(q₀), q₀, p₀; kwargs...)
    end

    """
    Lagrangian problem for the Toda lattice.
    """
    function lodeproblem(N::Int = Ñ, q₀ = compute_initial_q(μ, N), p₀ = zero(q₀); tspan = tspan, tstep = tstep, parameters = default_parameters)
        LODEProblem(lagrangian_system(N, parameters), tspan, tstep, q₀, p₀; parameters = parameters)
    end

    function lodeproblem(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        lodeproblem(length(q₀), q₀, p₀; kwargs...)
    end

    function lodeensemble(N::Int = Ñ, q₀ = compute_initial_q(μ, N), p₀ = zero(q₀); tspan = tspan, tstep = tstep, parameters = default_parameters)
        eqs = functions(lagrangian_system(N, _parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    function lodeensemble(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        lodeensemble(length(q₀), q₀, p₀; kwargs...)
    end

end
