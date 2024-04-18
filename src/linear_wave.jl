@doc raw"""
The discretized version of the 1d linear wave equation.

It is a prime example of a non-trivial completely integrable system.

The only system parameters are the *number of points* ``N`` for which the system is discretized and `\alpha`.
"""
module LinearWave 

    using EulerLagrange
    using LinearAlgebra 
    using Parameters 

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem  

    include("bump_initial_condition.jl")
    include("assemble_matrix_for_linear_wave_equation.jl")

    const μ̃ = .6
    const Ñ = 128

    const default_parameters = (
        μ = μ̃, 
        N = Ñ
    )

    function hamiltonian(t, q, p, params)
        @unpack N, μ = params
        
        Δx = typeof(μ)(1) / (Ñ + 1)
        K = assemble_matrix(μ, Ñ)
        sum(Δx * p[n] ^ 2 / 2 for n in 1 : Ñ) + q ⋅ (K.parent * q)
    end

    function lagrangian(t, q, q̇, params)
        @unpack N, α = params 

        sum(q̇[n] ^ 2 / (2 * Δx) for n in 1 : Ñ) - q ⋅ (K.parent * q)
    end

    const tstep = .1 
    const tspan = (0.0, 120.0)

    const q̃₀ = get_initial_condition(μ̃, Ñ).q 
    const p̃₀ = get_initial_condition(μ̃, Ñ).p 

    """
    Hamiltonian problem for the linear wave equation.
    """
    function hodeproblem(q₀ = q̃₀, p₀ = p̃₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(Ñ + 2)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    """
    Lagrangian problem for the linear wave equation.
    """
    function lodeproblem(q₀ = q̃₀, p₀ = p̃₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, x, v = lagrangian_variables(Ñ + 2)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        lodeproblem(lag_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

end
