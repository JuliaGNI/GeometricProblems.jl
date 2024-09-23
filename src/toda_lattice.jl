@doc raw"""
The Toda lattice is a model for a one-dimensional crystal named after its discoverer Morikazu Toda [toda1967vibration](@cite).

It is a prime example of a non-trivial completely integrable system.

The only system parameters are the *number of points* ``N`` in the periodic lattice and ``\alpha`` which adjusts the strength of the interactions in the lattice.
"""
module TodaLattice 

    using EulerLagrange
    using LinearAlgebra 
    using Parameters 
    using GeometricEquations: HODEEnsemble

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem  
    export hodeensemble

    include("bump_initial_condition.jl")

    const α̃  = .64
    const Ñ = 200

    const default_parameters = (
        α = α̃, 
        N = Ñ
    )

    function hamiltonian(t, q, p, params)
        @unpack N, α = params
        
        sum(p[n] ^ 2 / 2 + α * exp(q[n] - q[n % Ñ + 1]) for n in 1:Ñ)
    end

    function lagrangian(t, q, q̇, params)
        @unpack N, α = params 

        sum(q̇[n] ^ 2 / 2 - α * exp(q[n] - q[n % Ñ + 1]) for n in 1:Ñ)
    end

    const tstep = .1 
    const tspan = (0.0, 120.0)

    # parameter for the initial conditions
    const μ = .3

    const q₀ = compute_initial_condition(μ, Ñ).q 
    const p₀ = compute_initial_condition(μ, Ñ).p 
    const Ω = compute_domain(N̄, Float64)

    """
    Hamiltonian problem for the Toda lattice.
    """
    function hodeproblem(q₀ = q₀, p₀ = p₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(Ñ)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    """
    Lagrangian problem for the Toda lattice.
    """
    function lodeproblem(q₀ = q₀, p₀ = p₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, x, v = lagrangian_variables(Ñ)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        lodeproblem(lag_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    function hodeensemble(q₀ = q₀, p₀ = p₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        eq = hodeproblem().equation

        HODEEnsemble(eq.v, eq.f, eq.hamiltonian, tspan, tstep, q₀, p₀; parameters = parameters)
    end

end
