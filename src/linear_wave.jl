@doc raw"""
The discretized version of the 1d linear wave equation.

It is a prime example of a non-trivial completely integrable system.

The system is discretized on `N` interior points — the ``\tilde{N}`` of the documentation page —
so the state has `N + 2` components once the two boundary points are included.

Like the number of lattice sites of `TodaLattice`, `N` fixes the *size* of the system: it sets the
number of degrees of freedom and the summation bounds, so it cannot survive `symbolize` and hence
cannot be a system parameter. It is instead a leading positional argument of the problem
constructors, defaulting to `Ñ = 256`, and `hamiltonian`/`lagrangian` take it as a trailing
argument. The only system parameter proper is ``\mu``.
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
    const Ñ = 256

    default_parameters(::Type{T}=Float64) where {T} = (μ = T(μ̃),)

    # The discretized gradient energy on N interior points. Both the Hamiltonian and the Lagrangian
    # take N as a trailing argument (as in TodaLattice) rather than reading it from `parameters`:
    # it fixes the number of degrees of freedom and the summation bounds, so it has to be a plain
    # integer and cannot survive `symbolize`.
    function potential(q, parameters, N)
        @unpack μ = parameters

        Δx² = (one(μ) / (N + 1)) ^ 2
        μ ^ 2 / 4Δx² * sum(((q[i] - q[i - 1]) ^ 2 + (q[i + 1] - q[i]) ^ 2) for i in 2 : (N + 1))
    end

    function hamiltonian(t, q, p, parameters, N)
        sum(p[n] ^ 2 for n in 1 : (N + 2)) / 2 + potential(q, parameters, N)
    end

    function lagrangian(t, q, q̇, parameters, N)
        sum(q̇[n] ^ 2 for n in 1 : (N + 2)) / 2 - potential(q, parameters, N)
    end

    _timestep(timespan::Tuple, n_time_steps::Integer) = (timespan[2] - timespan[1]) / (n_time_steps-1)


    const DEFAULT_TIMESPAN = (0, 1)
    const n_time_steps = 200
    const DEFAULT_TIMESTEP = _timestep(DEFAULT_TIMESPAN, n_time_steps)

    const q₀ = compute_initial_condition2(μ̃, Ñ + 2).q 
    const p₀ = compute_initial_condition2(μ̃, Ñ + 2).p 

    function hamiltonian_system(N::Int, parameters::NamedTuple)
        t, q, p = hamiltonian_variables(N + 2)
        sparams = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams, N), t, q, p, sparams; simplify = false)
    end

    function lagrangian_system(N::Int, parameters::NamedTuple)
        t, x, v = lagrangian_variables(N + 2)
        sparams = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams, N), t, x, v, sparams; simplify = false)
    end

    # Build the symbolic system from a single parameter set, while a vector of parameter
    # sets is passed on to the ensemble unchanged (see issue #64).
    _parameters(p::NamedTuple) = p
    _parameters(p::AbstractVector) = p[begin]

    # The state has N + 2 components (N interior points plus the two boundary points), so an
    # initial condition of length n corresponds to N = n - 2.
    _nint(q::AbstractArray{<:Number}) = length(q) - 2
    _nint(q::AbstractVector{<:AbstractArray}) = length(q[begin]) - 2

    """
    Hamiltonian problem for the linear wave equation on `N` interior points.
    """
    function hodeproblem(N::Int = Ñ, q₀ = compute_initial_condition2(μ̃, N + 2).q, p₀ = compute_initial_condition2(μ̃, N + 2).p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        HODEProblem(hamiltonian_system(N, parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function hodeproblem(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        hodeproblem(_nint(q₀), q₀, p₀; kwargs...)
    end

    """
    Hamiltonian ensemble for the linear wave equation (varying initial conditions and/or parameters).
    """
    function hodeensemble(N::Int = Ñ, q₀ = compute_initial_condition2(μ̃, N + 2).q, p₀ = compute_initial_condition2(μ̃, N + 2).p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        eqs = functions(hamiltonian_system(N, _parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function hodeensemble(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        hodeensemble(_nint(q₀), q₀, p₀; kwargs...)
    end

    """
    Lagrangian problem for the linear wave equation on `N` interior points.
    """
    function lodeproblem(N::Int = Ñ, q₀ = compute_initial_condition2(μ̃, N + 2).q, p₀ = compute_initial_condition2(μ̃, N + 2).p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        LODEProblem(lagrangian_system(N, parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function lodeproblem(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        lodeproblem(_nint(q₀), q₀, p₀; kwargs...)
    end

    """
    Lagrangian ensemble for the linear wave equation (varying initial conditions and/or parameters).
    """
    function lodeensemble(N::Int = Ñ, q₀ = compute_initial_condition2(μ̃, N + 2).q, p₀ = compute_initial_condition2(μ̃, N + 2).p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        eqs = functions(lagrangian_system(N, _parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function lodeensemble(q₀, p₀; kwargs...)
        @assert length(q₀) == length(p₀)
        lodeensemble(_nint(q₀), q₀, p₀; kwargs...)
    end

end
