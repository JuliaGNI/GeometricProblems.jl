@doc raw"""
The Toda lattice is a model for a one-dimensional crystal named after its discoverer Morikazu Toda [toda1967vibration](@cite).

It is a prime example of a non-trivial completely integrable system.

The only system parameters are the *number of points* ``N`` in the periodic lattice and `\alpha` which adjusts the strength of the interactions in the lattice.
"""
module TodaLattice 

    using EulerLagrange
    using LinearAlgebra 
    using Parameters 

    export hamiltonian
    export hodeproblem 

    include("initial_conditions.jl")

    const α̃ = .5 
    const Ñ = 100

    const default_parameters = (
        α = α̃, 
        N = Ñ
    )

    function hamiltonian(t, q, p, params)
        @unpack α = params
        
        sum(p[n] ^ 2 / 2 + α ^ 2 * exp(q[n] - q[n % Ñ + 1]) for n in 1:Ñ)
    end

    const tstep = .01 
    const tspan = (0.0, 10.0)

    # parameter for the initial conditions
    const μ = .3

    const q̃₀ = get_initial_condition(μ, Ñ).q 
    const p̃₀ = get_initial_condition(μ, Ñ).p 

    """
    Hamiltonian problem for the Toda lattice.
    """
    function hodeproblem(q₀ = q̃₀, p₀ = p̃₀; tspan = tspan, tstep = tstep, params = default_parameters)
        t, q, p = hamiltonian_variables(Ñ)
        sparams = symbolize(params)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = sparams)
    end

end