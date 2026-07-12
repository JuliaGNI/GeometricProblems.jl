@doc raw"""
The discretized version of the 1d linear wave equation.

It is a prime example of a non-trivial completely integrable system.

The only system parameters are the *number of points* ``N`` for which the system is discretized and ``\mu``.
"""
module LinearWave 

    using EulerLagrange
    using LinearAlgebra
    using Parameters
    using GeometricEquations: HODEEnsemble, LODEEnsemble

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hodeensemble, lodeensemble
    export hamiltonian_system, lagrangian_system

    include("bump_initial_condition.jl")

    const μ̃ = .6
    const Ñ = 256

    const default_parameters = (μ = μ̃, N = Ñ)

    function hamiltonian(t, q, p, parameters)
        @unpack N, μ = parameters
        
        Δx = one(μ) / (Ñ + 1)
        Δx² = Δx ^ 2
        μ² = μ ^ 2
        sum(p[n] ^ 2 for n in 1 : (Ñ + 2)) / 2 + μ² / 4Δx² * sum(((q[i] - q[i - 1]) ^ 2 + (q[i + 1] - q[i]) ^ 2) for i in 2 : (Ñ + 1))   
    end

    function lagrangian(t, q, q̇, parameters)
        @unpack N, μ = parameters 

        Δx = one(μ) / (Ñ + 1)
        Δx² = Δx ^ 2 
        μ² = μ ^ 2
        sum(q̇[n] ^ 2 for n in 1 : (Ñ + 2)) / 2 - μ² / 4Δx² * sum(((q[i] - q[i - 1]) ^ 2 + (q[i + 1] - q[i]) ^ 2) for i in 2 : (Ñ + 1))
    end

    _timestep(timespan::Tuple, n_time_steps::Integer) = (timespan[2] - timespan[1]) / (n_time_steps-1)


    const DEFAULT_TIMESPAN = (0, 1)
    const n_time_steps = 200
    const DEFAULT_TIMESTEP = _timestep(DEFAULT_TIMESPAN, n_time_steps)

    const q₀ = compute_initial_condition2(μ̃, Ñ + 2).q 
    const p₀ = compute_initial_condition2(μ̃, Ñ + 2).p 

    function hamiltonian_system(parameters::NamedTuple)
        t, q, p = hamiltonian_variables(Ñ + 2)
        sparams = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams; simplify = false)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(Ñ + 2)
        sparams = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams; simplify = false)
    end

    # Build the symbolic system from a single parameter set, while a vector of parameter
    # sets is passed on to the ensemble unchanged (see issue #64).
    _parameters(p::NamedTuple) = p
    _parameters(p::AbstractVector) = p[begin]

    """
    Hamiltonian problem for the linear wave equation.
    """
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    """
    Hamiltonian ensemble for the linear wave equation (varying initial conditions and/or parameters).
    """
    function hodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        eqs = functions(hamiltonian_system(_parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    """
    Lagrangian problem for the linear wave equation.
    """
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    """
    Lagrangian ensemble for the linear wave equation (varying initial conditions and/or parameters).
    """
    function lodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        eqs = functions(lagrangian_system(_parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀; parameters = parameters)
    end

end
