@doc raw"""
    DoublePendulum

The `DoublePendulum` module provides functions `hodeproblem` and `lodeproblem`
each returning a Hamiltonian or Lagrangian problem, respectively, to be solved
in the GeometricIntegrators.jl ecosystem.
The actual code is generated with EulerLagrange.jl.

The double pendulum consists of two pendula, one attached to the origin at
``(x,y) = (0,0)``, and the second attached to the first. Each pendulum consists
of a point mass ``m_i`` attached to a massless rod of length ``l_i`` with ``i \in (1,2)``.
The dynamics of the system is described in terms of the angles ``\theta_i``
between the rods ``l_i`` and the vertical axis ``y``.
All motion is assumed to be frictionless.

System parameters:
* `l₁`: length of rod 1
* `l₂`: length of rod 2
* `m₁`: mass of pendulum 1
* `m₂`: mass of pendulum 2
* `g`: gravitational constant
"""
module DoublePendulum

    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem


    ϑ₁(t, q, q̇, params) = (params.m₁ + params.m₂) * params.l₁^2 * q̇[1] + 
        params.m₂ * params.l₁ * params.l₂ * q̇[2] * cos(q[1] - q[2])
    
    ϑ₂(t, q, q̇, params) = params.m₂ * params.l₂^2 * q̇[2] +
        params.m₂ * params.l₁ * params.l₂ * q̇[1] * cos(q[1] - q[2])

    ϑ(t, q, q̇, params) = [ϑ₁(t, q, q̇, params), ϑ₂(t, q, q̇, params)]


    function θ̇₁(t, q, p, params)
        @unpack l₁, l₂, m₁, m₂, g = params
        
        ( l₂ * p[1] - l₁ * p[2] * cos(q[1] - q[2]) ) /
        ( l₁^2 * l₂ * ( m₁ + m₂ * sin(q[1] - q[2])^2 ) )
    end

    function θ̇₂(t, q, p, params)
        @unpack l₁, l₂, m₁, m₂, g = params
        
        ( (m₁ + m₂) * l₁ * p[2] - m₂ * l₂ * p[1] * cos(q[1] - q[2]) ) /
        ( m₂ * l₁ * l₂^2 * ( m₁ + m₂ * sin(q[1] - q[2])^2 ) )
    end

    θ̇(t, q, p, params) = [θ̇₁(t, q, p, params), θ̇₂(t, q, p, params)]

    function θ̇(v, t, q, p, params)
        v[1] = θ̇₁(t, q, p, params)
        v[2] = θ̇₂(t, q, p, params)
        nothing
    end


    const tstep = 0.01
    const tspan = (0.0, 10.0)

    const default_parameters = (
        l₁ = 2.0,
        l₂ = 3.0,
        m₁ = 1.0,
        m₂ = 2.0,
        g = 9.80665,
    )

    const θ₀ = [π/4, π/2]
    const ω₀ = [0.0, π/8]
    const p₀ = ϑ(tspan[begin], θ₀, ω₀, default_parameters)


    function hamiltonian(t, q, p, params)
        @unpack l₁, l₂, m₁, m₂, g = params

        nom = (m₁ + m₂) * l₁^2 * p[2]^2 / 2+ 
                    m₂  * l₂^2 * p[1]^2 / 2 - 
              m₂ * l₁ * l₂ * p[1] * p[2] * cos(q[1] - q[2])
        
        den = m₂ * l₁^2 * l₂^2 * ( m₁ + m₂ * sin(q[1] - q[2])^2 )
        
        nom/den - g * (m₁ + m₂) * l₁ * cos(q[1]) - g * m₂ * l₂ * cos(q[2])
    end


    function lagrangian(t, q, q̇, params)
        @unpack l₁, l₂, m₁, m₂, g = params

        (m₁ + m₂) * l₁^2 * q̇[1]^2 / 2 + 
              m₂  * l₂^2 * q̇[2]^2 / 2 +
         m₂ * l₁  * l₂ * q̇[1] * q̇[2] * cos(q[1] - q[2]) +
        (m₁ + m₂) * l₁ * g * cos(q[1]) +
              m₂  * l₂ * g * cos(q[2])
    end


    """
        Hamiltonian problem for the double pendulum

    Constructor with default arguments:
    ```
    hodeproblem(
        q₀ = [π/4, π/2],
        p₀ = $(p₀);
        tspan = $(tspan),
        tstep = $(tstep),
        params = $(default_parameters)
    )
    ```
    """
    function hodeproblem(q₀ = θ₀, p₀ = p₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(2)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    """
        Lagrangian problem for the double pendulum

    Constructor with default arguments:
    ```
    lodeproblem(
        q₀ = [π/4, π/2],
        p₀ = $(p₀);
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
