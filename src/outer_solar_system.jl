@doc raw"""
    OuterSolarSystem

Gravitational N-body model of the Sun and the five outer planets (Jupiter, Saturn,
Uranus, Neptune) and Pluto. The system has ``N = 6`` bodies in ``d = 3`` spatial
dimensions, giving 18 degrees of freedom.

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

    function hamiltonian(t, q, p, params)
        @unpack G, m = params
        Q  = reshape(q, N_DIM, N_BODIES)
        P  = reshape(p, N_DIM, N_BODIES)
        KE = sum(dot(P[:,i], P[:,i]) / m[i] for i in 1:N_BODIES) / 2
        V  = sum(G * m[i] * m[j] / norm(Q[:,i] - Q[:,j])
                 for i in 2:N_BODIES for j in 1:i-1)
        KE - V
    end

    function lagrangian(t, q, q̇, params)
        @unpack G, m = params
        Q  = reshape(q,  N_DIM, N_BODIES)
        Q̇  = reshape(q̇, N_DIM, N_BODIES)
        KE = sum(m[i] * dot(Q̇[:,i], Q̇[:,i]) for i in 1:N_BODIES) / 2
        V  = sum(G * m[i] * m[j] / norm(Q[:,i] - Q[:,j])
                 for i in 2:N_BODIES for j in 1:i-1)
        KE + V
    end

    function hamiltonian_system(parameters::NamedTuple)
        t, q, p = hamiltonian_variables(N_BODIES * N_DIM)
        sparams  = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(N_BODIES * N_DIM)
        sparams  = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
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
        parameters = $(default_parameters())
    )
    ```
    """
    function hodeproblem(q₀ = q₀, p₀ = p₀;
                         timespan   = DEFAULT_TIMESPAN,
                         timestep   = DEFAULT_TIMESTEP,
                         parameters = default_parameters())
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀;
                    parameters = parameters)
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
        parameters = $(default_parameters())
    )
    ```
    """
    function lodeproblem(q₀ = q₀, p₀ = p₀;
                         timespan   = DEFAULT_TIMESPAN,
                         timestep   = DEFAULT_TIMESTEP,
                         parameters = default_parameters())
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀;
                    v̄ = v̄, parameters = parameters)
    end

end
