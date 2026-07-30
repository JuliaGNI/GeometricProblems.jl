@doc raw"""
    OuterSolarSystem

Gravitational N-body model of the Sun, the four outer planets (Jupiter, Saturn,
Uranus, Neptune) and Pluto. The system has ``N = 6`` bodies in ``d = 3`` spatial
dimensions, giving 18 degrees of freedom.

`hodeproblem` and `lodeproblem` use hand-written vector fields by default, in which each of the
``N(N-1)/2`` pairwise distances is computed once per force evaluation. Passing `symbolic = true`
generates the equations of motion with EulerLagrange instead, which agrees to round-off and is
kept for cross-checking; see `benchmark/outer_solar_system.jl` for the comparison.

The Hamiltonian is
```math
H(q, p) = \frac{1}{2} \sum_{i=1}^{N} \frac{p_i \cdot p_i}{m_i}
         - G \sum_{1 \le j < i \le N} \frac{m_i \, m_j}{\| q_i - q_j \|}
```
where ``q_i, p_i \in \mathbb{R}^3`` are the position and momentum of body ``i``,
and ``G`` is the gravitational constant.

System parameters:
* `G`: gravitational constant in units A.U.³ / (M_☉ · day²)
* `m`: vector of body masses ``[m_1, \ldots, m_6]`` relative to the solar mass

Initial conditions (positions in A.U., velocities in A.U./day) and constants are
taken from E. Hairer, C. Lubich, G. Wanner, *Geometric Numerical Integration*, Chapter I.2.4.
"""
module OuterSolarSystem

    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hamiltonian_system, lagrangian_system
    export default_parameters

    const N_BODIES = 6
    const N_DIM    = 3

    const DEFAULT_TIMESTEP = 0.5
    const DEFAULT_TIMESPAN = (0.0, 200.0)

    # Sun
    const m₁ = 1.0
    const q₁ = [0.0,        0.0,        0.0       ]
    const q̇₁ = [0.0,        0.0,        0.0       ]

    # Jupiter
    const m₂ = 0.000954786104043
    const q₂ = [-3.5023653,  -3.8169847,  -1.5507963]
    const q̇₂ = [ 0.00565429,  -0.0041249,  -0.00190589]

    # Saturn
    const m₃ = 0.000285583733151
    const q₃ = [ 9.0755314,  -3.0458353,  -1.6483708]
    const q̇₃ = [ 0.00168318,  0.00483525,  0.00192462]

    # Uranus
    const m₄ = 0.0000437273164546
    const q₄ = [ 8.3101420, -16.2901086,  -7.2521278]
    const q̇₄ = [ 0.00354178,  0.00137102,  0.00055029]

    # Neptune
    const m₅ = 0.0000517759138449
    const q₅ = [11.4707666, -25.7294829, -10.8169456]
    const q̇₅ = [ 0.00288930,  0.00114527,  0.00039677]

    # Pluto
    const m₆ = 1.0 / 1.3e8
    const q₆ = [-15.5387357, -25.2225594,  -3.1902382]
    const q̇₆ = [  0.00276725,  -0.00170702,  -0.00136504]

    const q₀ = [q₁; q₂; q₃; q₄; q₅; q₆]
    const p₀ = [m₁ * q̇₁; m₂ * q̇₂; m₃ * q̇₃; m₄ * q̇₄; m₅ * q̇₅; m₆ * q̇₆]

    @doc raw"Default gravitational constant G in A.U.³ / (M_☉ · day²)."
    const G = 2.95912208286e-4

    default_parameters(::Type{T} = Float64) where {T} = (
        G = T(G),
        m = T[m₁, m₂, m₃, m₄, m₅, m₆],
    )

    # In-place velocity estimate from momentum: q̇ᵢ = pᵢ / mᵢ (diagonal mass matrix)
    function v̄(v, t, q, p, params)
        @unpack m = params
        for i in 1:N_BODIES
            v[N_DIM*i-2] = p[N_DIM*i-2] / m[i]
            v[N_DIM*i-1] = p[N_DIM*i-1] / m[i]
            v[N_DIM*i  ] = p[N_DIM*i  ] / m[i]
        end
        nothing
    end

    function T(p::AbstractVector, params::NamedTuple)
        @unpack m = params
        P = reshape(p, N_DIM, N_BODIES)
        s = zero(eltype(p))
        for i in 1:N_BODIES
            s += (P[1,i]^2 + P[2,i]^2 + P[3,i]^2) / m[i]
        end
        s / 2
    end

    function V(q::AbstractVector, params::NamedTuple)
        @unpack G, m = params
        Q = reshape(q, N_DIM, N_BODIES)
        s = zero(eltype(q))
        for i in 2:N_BODIES
            for j in 1:i-1
                Δ₁ = Q[1,i] - Q[1,j]
                Δ₂ = Q[2,i] - Q[2,j]
                Δ₃ = Q[3,i] - Q[3,j]
                s += G * m[i] * m[j] / sqrt(Δ₁^2 + Δ₂^2 + Δ₃^2)
            end
        end
        -s
    end

    @doc raw"""
        ∇V!(dV, q, params)

    Gradient ``∂V/∂q`` of the gravitational potential, written into `dV`.

    Each pairwise distance ``r_{ij}`` is computed once and its contribution accumulated into both
    bodies at once (Newton's third law), so the whole gradient costs ``N(N-1)/2`` square roots.
    This is the kernel shared by the Hamiltonian and Lagrangian force functions.
    """
    function ∇V!(dV, q, params)
        @unpack G, m = params
        fill!(dV, zero(eltype(dV)))
        for i in 2:N_BODIES
            for j in 1:i-1
                Δ₁ = q[N_DIM*i-2] - q[N_DIM*j-2]
                Δ₂ = q[N_DIM*i-1] - q[N_DIM*j-1]
                Δ₃ = q[N_DIM*i  ] - q[N_DIM*j  ]
                r² = Δ₁^2 + Δ₂^2 + Δ₃^2
                # ∂/∂qᵢ (-G mᵢ mⱼ / r) = G mᵢ mⱼ Δ / r³, and the opposite sign for body j
                c  = G * m[i] * m[j] / (r² * sqrt(r²))
                dV[N_DIM*i-2] += c * Δ₁
                dV[N_DIM*i-1] += c * Δ₂
                dV[N_DIM*i  ] += c * Δ₃
                dV[N_DIM*j-2] -= c * Δ₁
                dV[N_DIM*j-1] -= c * Δ₂
                dV[N_DIM*j  ] -= c * Δ₃
            end
        end
        nothing
    end

    function hamiltonian(t, q, p, params)
        T(p, params) + V(q, params)
    end

    # kinetic energy in velocity form: Σ mᵢ q̇ᵢ·q̇ᵢ / 2  (not the momentum form T(p))
    function lagrangian(t, q, q̇, params)
        @unpack m = params
        Q̇ = reshape(q̇, N_DIM, N_BODIES)
        s = zero(eltype(q̇))
        for i in 1:N_BODIES
            s += m[i] * (Q̇[1,i]^2 + Q̇[2,i]^2 + Q̇[3,i]^2)
        end
        s / 2 - V(q, params)
    end

    # ---------------------------------------------------------------------------------------------
    # Hand-written vector fields
    #
    # These are what `hodeproblem`/`lodeproblem` use by default. They follow the pattern of the
    # other hand-written problems in this package (see `lotka_volterra_2d_common.jl` and
    # `harmonic_oscillator.jl`), including the extra arity methods GeometricEquations calls.
    # ---------------------------------------------------------------------------------------------

    # ∂H/∂pᵢ = pᵢ / mᵢ — the same expression as the initial guess `v̄`.
    outer_solar_system_v(v, t, q, p, params) = v̄(v, t, q, p, params)

    # -∂H/∂q = -∂V/∂q, and ∂L/∂q = -∂V/∂q as well, so the Hamiltonian and Lagrangian force
    # functions coincide and differ only in whether the third slot holds p or q̇.
    function outer_solar_system_f(f, t, q, w, params)
        ∇V!(f, q, params)
        f .= .-f
        nothing
    end

    # ϑ = ∂L/∂q̇ = mᵢ q̇ᵢ
    function outer_solar_system_ϑ(Θ, t, q, w, params)
        @unpack m = params
        for i in 1:N_BODIES
            Θ[N_DIM*i-2] = m[i] * w[N_DIM*i-2]
            Θ[N_DIM*i-1] = m[i] * w[N_DIM*i-1]
            Θ[N_DIM*i  ] = m[i] * w[N_DIM*i  ]
        end
        nothing
    end

    # Projection direction for enforcing p = ϑ(q, q̇). The Lagrangian is regular, so the constraint
    # is non-degenerate and the projection acts in every direction; as in `harmonic_oscillator.jl`
    # the multiplier is taken to absorb the mass matrix, giving g = λ.
    function outer_solar_system_g(g, t, q, λ, params)
        g .= λ
        nothing
    end

    outer_solar_system_g(g, t, q, w, λ, params) = outer_solar_system_g(g, t, q, λ, params)

    @doc raw"""
        ω!(Ω, t, q, params)

    The Lagrange two-form ``ω = dθ`` of this regular Lagrangian, where ``θ = ϑ_i \, dq^i``.

    A regular Lagrangian is a second-order system of ``n`` equations, equivalent to a first-order
    system of ``2n``, so `Ω` is the ``2n × 2n`` form on ``(q, \dot q)``. In block form it is
    ``[∂ϑ_i/∂q_j - ∂ϑ_j/∂q_i \; -M^T; \; M \; 0]``; since ``ϑ = ∂L/∂\dot q = m_i \dot q_i`` depends
    on the velocities alone, the upper-left block vanishes and what remains is the canonical
    ``[0 \; -M; \; M \; 0]`` with ``M = \mathrm{diag}(m)`` repeated over the spatial components.

    This is the convention EulerLagrange's `LagrangianSystem` produces, so `lodeproblem` agrees
    whether or not `symbolic = true`.
    """
    function ω!(Ω, t, q, params)
        @unpack m = params
        fill!(Ω, zero(eltype(Ω)))
        D = N_BODIES * N_DIM
        for i in 1:N_BODIES
            for k in 0:N_DIM-1
                j = N_DIM * i - 2 + k
                Ω[j, D+j] = -m[i]
                Ω[D+j, j] = +m[i]
            end
        end
        nothing
    end

    ω!(Ω, t, q, w, params) = ω!(Ω, t, q, params)

    # ---------------------------------------------------------------------------------------------
    # Symbolic formulation (EulerLagrange)
    # ---------------------------------------------------------------------------------------------

    # EulerLagrange defaults to `simplify = false`, which matters a great deal for this problem:
    # `Symbolics.simplify` costs almost nothing on the Lagrangian itself, but it rewrites the sum of
    # 15 inverse-distance terms into a common-denominator form whose derivatives are vastly more
    # expensive to build, to compile and to evaluate. `LagrangianSystem` takes 0.31 s as it stands,
    # against 40 s with `simplify = true` (and 277 s with `simplify = true, cse = false`, which used
    # to be the default), and the generated force is 12× slower to call.
    function hamiltonian_system(parameters::NamedTuple)
        t, q, p = hamiltonian_variables(N_BODIES * N_DIM)
        sparams  = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams; nanmath = true)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(N_BODIES * N_DIM)
        sparams  = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams; nanmath = true)
    end

    """
        hodeproblem(q₀, p₀; timespan, timestep, parameters)

    Hamiltonian formulation of the outer solar system
    (Sun, Jupiter, Saturn, Uranus, Neptune, Pluto).

    Constructor with default arguments:
    ```
    hodeproblem(
        q₀ = q₀,
        p₀ = p₀;
        timespan  = $(DEFAULT_TIMESPAN),
        timestep  = $(DEFAULT_TIMESTEP),
        parameters = $(default_parameters()),
        symbolic  = false
    )
    ```

    With `symbolic = true` the equations of motion are generated with EulerLagrange via
    [`hamiltonian_system`](@ref) instead of using the hand-written vector fields. The two agree to
    round-off; the symbolic route is kept for cross-checking and costs a few seconds to build.
    """
    function hodeproblem(q₀ = q₀, p₀ = p₀;
                         timespan   = DEFAULT_TIMESPAN,
                         timestep   = DEFAULT_TIMESTEP,
                         parameters = default_parameters(),
                         symbolic   = false)
        if symbolic
            HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀;
                        parameters = parameters)
        else
            HODEProblem(outer_solar_system_v, outer_solar_system_f, hamiltonian,
                        timespan, timestep, q₀, p₀; parameters = parameters)
        end
    end

    """
        lodeproblem(q₀, p₀; timespan, timestep, parameters)

    Lagrangian formulation of the outer solar system
    (Sun, Jupiter, Saturn, Uranus, Neptune, Pluto).

    Constructor with default arguments:
    ```
    lodeproblem(
        q₀ = q₀,
        p₀ = p₀;
        timespan   = $(DEFAULT_TIMESPAN),
        timestep   = $(DEFAULT_TIMESTEP),
        parameters = $(default_parameters()),
        symbolic   = false
    )
    ```

    With `symbolic = true` the equations of motion are generated with EulerLagrange via
    [`lagrangian_system`](@ref) instead of using the hand-written vector fields. The two agree to
    round-off; the symbolic route is kept for cross-checking and costs a few seconds to build.
    """
    function lodeproblem(q₀ = q₀, p₀ = p₀;
                         timespan   = DEFAULT_TIMESPAN,
                         timestep   = DEFAULT_TIMESTEP,
                         parameters = default_parameters(),
                         symbolic   = false)
        if symbolic
            LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀;
                        v̄ = v̄, parameters = parameters)
        else
            LODEProblem(outer_solar_system_ϑ, outer_solar_system_f, outer_solar_system_g,
                        ω!, lagrangian, timespan, timestep, q₀, p₀;
                        v̄ = v̄, parameters = parameters)
        end
    end

end
