@doc raw"""
The discretized version of the 1d linear wave equation.

It is a prime example of a non-trivial completely integrable system.

The only system parameters are the *number of points* ``N`` for which the system is discretized and ``\mu``.
"""
module LinearWave 

    using EulerLagrange
    using LinearAlgebra 
    using Parameters 
    using GeometricIntegrators

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem  

    include("bump_initial_condition.jl")

    const μ̃ = 0.5
    const Ñ = 256

    const default_parameters = (μ = μ̃, N = Ñ)

    function hamiltonian(t, q, p, parameters)
        @unpack N, μ = parameters
        
        Δx = one(μ) / (Ñ + 1)
        Δx² = Δx ^ 2
        μ² = μ ^ 2
        # sum(p[n] ^ 2 / 2 for n in 1 : (Ñ + 2)) + sum(μ² / 4Δx² * ((q[i] - q[i - 1]) ^ 2 + (q[i + 1] - q[i]) ^ 2) for i in 2 : (Ñ + 1))   
        sum(p[n] ^ 2 / 2 for n in 1 : (Ñ + 2)) + sum(μ² / 2Δx² * ((q[i] - q[i - 1]) ^ 2) for i in 2 : (Ñ + 2)) +   μ² / 2Δx² * ((q[1] - q[Ñ + 2]) ^ 2)

    end

    function lagrangian(t, q, q̇, parameters)
        @unpack N, μ = parameters 

        Δx = one(μ) / (Ñ + 1)
        Δx² = Δx ^ 2 
        μ² = μ ^ 2
        # sum(q̇[n] ^ 2 / 2 for n in 1 : (Ñ + 2)) - sum(μ² / 4Δx² * ((q[i] - q[i - 1]) ^ 2 + (q[i + 1] - q[i]) ^ 2) for i in 2 : (Ñ + 1))
        sum(q̇[n] ^ 2 / 2 for n in 1 : (Ñ + 2)) - sum(μ² / 2Δx² * ((q[i] - q[i - 1]) ^ 2) for i in 2 : (Ñ + 2)) +   μ² / 2Δx² * ((q[1] - q[Ñ + 2]) ^ 2)

    end

    _tstep(tspan::Tuple, n_time_steps::Integer) = (tspan[2] - tspan[1]) / (n_time_steps-1)


    const tspan = (0, 1)
    const n_time_steps = 200
    const tstep = _tstep(tspan, n_time_steps)

    const q₀ = compute_initial_condition2(μ̃, Ñ + 2).q 
    const p₀ = compute_initial_condition2(μ̃, Ñ + 2).p 

    # const q₀ = initial_position(Ñ + 2)
    # const p₀ = initial_velocity(Ñ + 2)

    """
    Hamiltonian problem for the linear wave equation.
    """
    function hodeproblem(q₀ = q₀, p₀ = p₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(Ñ + 2)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    """
    Lagrangian problem for the linear wave equation.
    """
    function lodeproblem(q₀ = q₀, p₀ = p₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, x, v = lagrangian_variables(Ñ + 2)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        lodeproblem(lag_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

end