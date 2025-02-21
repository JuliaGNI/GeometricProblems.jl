@doc raw"""
    ThreeBody

System parameters:
* `m₁`: mass of body 1
* `m₂`: mass of body 2
"""
module ThreeBody

    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

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
    const tstep = .5

    @doc raw"Range is taken from [jin2020sympnets](@cite). In that reference the integration is only done for ten time steps."
    const tspan = (0.0, 10 * tstep)

    @doc raw"Default parameters taken from [jin2020sympnets](@cite)."
    const default_parameters = (
        m₁ = m₁,
        m₂ = m₂,
        m₃ = m₃,
        G = G
    )

    T(p::AbstractVector, params::NamedTuple) = (p[1] ^ 2 + p[2] ^ 2) / (2 * params.m₁) + (p[3] ^ 2 + p[4] ^ 2) / (2 * params.m₂) + (p[5] ^ 2 + p[6] ^ 2) / (2 * params.m₃)
    V(q::AbstractVector, params::NamedTuple) = -params.G * params.m₁ * params.m₂ / √((q[1] - q[3]) ^ 2 + (q[2] - q[4]) ^ 2) - params.G * params.m₂ * params.m₃ / √((q[3] - q[5]) ^ 2 + (q[4] - q[6]) ^ 2)

    function hamiltonian(t, q, p, params)
        T(p, params) + V(q, params)
    end


    function lagrangian(t, q, q̇, params)
        T(q̇, params) - V(q, params)
    end


    """
        hodeproblem(q₀, p₀; tspan, tstep, parameters)

    Hamiltonian version of the three-body problem

    Constructor with default arguments:
    ```
    hodeproblem(
        q₀ = $(initial_condition.q),
        p₀ = $(initial_condition.p);
        tspan = $(tspan),
        tstep = $(tstep),
        params = $(default_parameters)
    )
    ```
    """
    function hodeproblem(q₀ = initial_condition.q, p₀ = initial_condition.p; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(6)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    """
        lodeproblem(q₀, p₀; tspan, tstep, parameters)

    Lagrangian version of the three-body problem

    Constructor with default arguments:
    ```
    lodeproblem(
        q₀ = $(initial_condition.q),
        p₀ = $(initial_condition.p);
        tspan = $(tspan),
        tstep = $(tstep),
        params = $(default_parameters)
    )
    ```
    """
    function lodeproblem(q₀ = θ₀, p₀ = p₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, x, v = lagrangian_variables(2)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, tspan, tstep, q₀, p₀; v̄ = θ̇, parameters = parameters)
    end

end
