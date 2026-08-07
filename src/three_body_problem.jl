@doc raw"""
    ThreeBody

System parameters:
* `m₁`: mass of body 1
* `m₂`: mass of body 2
* `m₃`: mass of body 3
* `G`: gravitational constant

The default initial condition is the figure-eight choreography (`figure_eight`, aliased as
`initial_condition`) and the default window is one of its periods. The 4096-member
`initial_conditions` grid of [jin2020sympnets](@cite) is also provided, but every one of its members
ends in a collision — see [`sympnets_initial_condition`](@ref) — so those are usable only on a window
that stops short of it.
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
    @doc raw"""
    The first member of `initial_conditions`. It is **not** the default initial condition, because it
    ends in a collision: the bodies collide at ``t \approx 0.08867``, where the solution ceases to
    exist. Successive RK4 refinements resolve a smallest mutual distance of ``1.95 \times 10^{-3}``,
    ``5.39 \times 10^{-4}``, ``4.78 \times 10^{-5}`` and ``2.93 \times 10^{-6}`` for
    ``\Delta{}t = 10^{-4} \ldots 10^{-7}`` — the distance keeps collapsing at a fixed time rather
    than bottoming out — and two RK4 references at ``\Delta{}t = 10^{-6}`` and ``5 \times 10^{-7}``
    disagree by 36 in position, so there is no computable trajectory to compare against either.

    Every member of `initial_conditions` behaves this way: all 4096 come within 0.04 of a collision
    within ``t \in [0, 5]``. Integrating any of them past its collision is not a matter of the step
    size or of the nonlinear solver (see [`DEFAULT_TIMESTEP`](@ref)); it needs regularization
    (Kustaanheimo–Stiefel, Levi-Civita) or an adaptive time transformation, neither of which this
    package provides. Use it only on a window that ends before the collision.
    """
    const sympnets_initial_condition = initial_conditions[1]

    @doc raw"""
    The default initial condition: the figure-eight choreography of [chenciner2000remarkable](@cite),
    in which three equal masses chase one another along a single figure-eight curve.

    It is the default because, unlike every member of `initial_conditions` (see
    [`sympnets_initial_condition`](@ref)), it has no close encounters — the smallest mutual distance
    over a period is 0.69 — so it can be integrated over many periods. Over one period
    `DEFAULT_TIMESTEP` conserves energy to ``4 \times 10^{-9}`` and closes the orbit to
    ``1.3 \times 10^{-3}`` with `ImplicitMidpoint`.

    The momenta equal the velocities because the masses are one: ``v_3 = (-0.93240737, -0.86473146)``
    and ``v_1 = v_2 = -v_3 / 2``.
    """
    const figure_eight = (
        q = [ 0.97000436, -0.24308753, -0.97000436,  0.24308753,  0.0,         0.0        ],
        p = [ 0.46620368,  0.43236573,  0.46620368,  0.43236573, -0.93240737, -0.86473146],
    )

    @doc raw"Alias for [`figure_eight`](@ref), the default initial condition of the problem constructors."
    const initial_condition = figure_eight

    @doc raw"Period of the `figure_eight` choreography for ``G = m_1 = m_2 = m_3 = 1``."
    const figure_eight_period = 6.32591398

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const m₁ = 1.

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const m₂ = 1.

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const m₃ = 1.

    @doc raw"Constant taken from [jin2020sympnets](@cite)."
    const G = 1.

    @doc raw"""
    Four hundred steps per period of the default `figure_eight` initial condition. Over one period
    `ImplicitMidpoint` then conserves energy to ``4 \times 10^{-9}`` and closes the orbit to
    ``1.3 \times 10^{-3}``; `Gauss(2)` reaches ``3 \times 10^{-15}`` and ``1.3 \times 10^{-7}``.

    Neither this nor any other step size makes `sympnets_initial_condition` integrable past its
    collision, and neither does exchanging the nonlinear solver. Over ``t \in [0, 1]`` at
    ``\Delta{}t = 0.5, 0.05, 0.01, 0.005, 0.001`` the default Newton method with a backtracking line
    search commits an energy error of 19, 15, 36, 399 and 249, and the trust-region `DogLeg` solver
    one of 9, 45, 12, 1615 and ``2.2 \times 10^6`` — with the two ending up 0.6, 3.8, 5.0, 50 and
    1209 apart in ``q``. Neither is quieter than the other about it — they emit 5, 9, 4, 4, 4 and
    6, 1, 1, 5, 6 solver messages respectively — though the count is not the thing to compare on:
    repeated reports are capped by `maxlog` and a stagnating solve stops after `max_stalls`, so it
    measures how a failure is *reported* rather than how badly it fails. Since neither solver can
    integrate the collision and neither is systematically better behaved, the default Newton solver
    is kept; on problems that *are* solvable the two agree bit for bit.
    """
    const DEFAULT_TIMESTEP = figure_eight_period / 400

    @doc raw"One period of the default `figure_eight` initial condition."
    const DEFAULT_TIMESPAN = (0.0, figure_eight_period)

    @doc raw"Default parameters taken from [jin2020sympnets](@cite)."
    default_parameters(::Type{DT}=Float64) where {DT} = (
        m₁ = DT(m₁),
        m₂ = DT(m₂),
        m₃ = DT(m₃),
        G = DT(G)
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
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams; nanmath = true)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(6)
        sparams = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams; nanmath = true)
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
        parameters = $(default_parameters())
    )
    ```
    """
    function hodeproblem(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Hamiltonian ensemble for the three-body problem (varying initial conditions and/or parameters)."
    function hodeensemble(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
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
        parameters = $(default_parameters())
    )
    ```
    """
    function lodeproblem(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

    "Lagrangian ensemble for the three-body problem (varying initial conditions and/or parameters)."
    function lodeensemble(q₀ = initial_condition.q, p₀ = initial_condition.p; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        eqs = functions(lagrangian_system(_parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
