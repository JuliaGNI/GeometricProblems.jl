@doc raw"""
    ThreeBody

System parameters:
* `m₁`: mass of body 1
* `m₂`: mass of body 2
* `m₃`: mass of body 3
* `G`: gravitational constant
"""
module ThreeBody

    using EulerLagrange
    using LinearAlgebra
    using Parameters
    using GeometricEquations: HODEEnsemble, LODEEnsemble

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hodeensemble, lodeensemble
    export hamiltonian_system, lagrangian_system

    "Turn array into a vector."
    _reshape(arr::AbstractArray) = reshape(arr, length(arr))

    @doc raw"Here we get 4096 trajectories that should be similar to the data used in [jin2020sympnets](@cite)."
    const initial_conditions = [
                                (   q = [sin(π * i₁ / 5), cos(π * i₂ / 5), sin(π * j₁ / 5 + π / 2), sin(π * j₂ / 5), cos(π * k₁ / 5), sin(π * k₂ / 5 + π / 2)], 
                                    p = [sin(π * iᵖ₁ / 5), cos(π * iᵖ₂ / 5), sin(π * jᵖ₁ / 5 + π / 2), sin(π * jᵖ₂ / 5), cos(π * kᵖ₁ / 5), sin(π * kᵖ₂ / 5 + π / 2)]
                                ) 
                                for i₁ ∈ 1:2, i₂ ∈ 1:2, j₁ ∈ 1:2, j₂ ∈ 1:2, k₁ ∈ 1:2, k₂ ∈ 1:2, iᵖ₁ ∈ 1:2, iᵖ₂ ∈ 1:2, jᵖ₁ ∈ 1:2, jᵖ₂ ∈ 1:2, kᵖ₁ ∈ 1:2, kᵖ₂ ∈ 1:2] |> _reshape
    const initial_condition = initial_conditions[1]

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const m₁ = 1.

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const m₂ = 1.

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const m₃ = 1.

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const G = 1.

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const DEFAULT_TIMESTEP = .5

    @doc raw"Range is taken from [jin2020sympnets](@cite). In that reference the integration is only done for ten time steps."
    const DEFAULT_TIMESPAN = (0.0, 10 * DEFAULT_TIMESTEP)

    @doc raw"Default parameters taken from [jin2020sympnets](@cite)."
    const default_parameters = (
        m₁ = m₁,
        m₂ = m₂,
        m₃ = m₃,
        G = G
    )

    T(p::AbstractVector, params::NamedTuple) = (p[1] ^ 2 + p[2] ^ 2) / (2 * params.m₁) + (p[3] ^ 2 + p[4] ^ 2) / (2 * params.m₂) + (p[5] ^ 2 + p[6] ^ 2) / (2 * params.m₃)
    V(q::AbstractVector, params::NamedTuple) = -params.G * params.m₁ * params.m₂ / √((q[1] - q[3]) ^ 2 + (q[2] - q[4]) ^ 2) - params.G * params.m₂ * params.m₃ / √((q[3] - q[5]) ^ 2 + (q[4] - q[6]) ^ 2) - params.G * params.m₁ * params.m₃ / √((q[1] - q[5]) ^ 2 + (q[2] - q[6]) ^ 2)

    function hamiltonian(t, q, p, params)
        T(p, params) + V(q, params)
    end


    # kinetic energy in terms of the velocities, T = Σ mᵢ q̇ᵢ² / 2 (not the momentum form T(p))
    function lagrangian(t, q, q̇, params)
        ( params.m₁ * (q̇[1] ^ 2 + q̇[2] ^ 2)
        + params.m₂ * (q̇[3] ^ 2 + q̇[4] ^ 2)
        + params.m₃ * (q̇[5] ^ 2 + q̇[6] ^ 2) ) / 2 - V(q, params)
    end

    # initial guess for the velocity given the momentum, q̇ = M⁻¹ p
    function v̄(v, t, q, p, params)
        v[1] = p[1] / params.m₁
        v[2] = p[2] / params.m₁
        v[3] = p[3] / params.m₂
        v[4] = p[4] / params.m₂
        v[5] = p[5] / params.m₃
        v[6] = p[6] / params.m₃
        nothing
    end

    function hamiltonian_system(parameters::NamedTuple)
        t, q, p = hamiltonian_variables(6)
        sparams = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(6)
        sparams = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
    end

    # Build the symbolic system from a single parameter set, while a vector of parameter
    # sets is passed on to the ensemble unchanged (see issue #64).
    _parameters(p::NamedTuple) = p
    _parameters(p::AbstractVector) = p[begin]


    """
        hodeproblem(q₀, p₀; timespan, timestep, parameters)

    Hamiltonian version of the three-body problem

    Constructor with default arguments:
    ```
    hodeproblem(
        q₀ = $(initial_condition.q),
        p₀ = $(initial_condition.p);
        timespan = $(DEFAULT_TIMESPAN),
        timestep = $(DEFAULT_TIMESTEP),
        params = $(default_parameters)
    )
    ```
    """
    function hodeproblem(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Hamiltonian ensemble for the three-body problem (varying initial conditions and/or parameters)."
    function hodeensemble(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        eqs = functions(hamiltonian_system(_parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    """
        lodeproblem(q₀, p₀; timespan, timestep, parameters)

    Lagrangian version of the three-body problem

    Constructor with default arguments:
    ```
    lodeproblem(
        q₀ = $(initial_condition.q),
        p₀ = $(initial_condition.p);
        timespan = $(DEFAULT_TIMESPAN),
        timestep = $(DEFAULT_TIMESTEP),
        params = $(default_parameters)
    )
    ```
    """
    function lodeproblem(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

    "Lagrangian ensemble for the three-body problem (varying initial conditions and/or parameters)."
    function lodeensemble(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        eqs = functions(lagrangian_system(_parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
