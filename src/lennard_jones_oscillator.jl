@doc raw"""
# Lennard-Jones Oscillator

A single particle oscillating in a one-dimensional Lennard-Jones potential,
```math
V(q) = 4 \varepsilon \left[ \left( \frac{\sigma}{q} \right)^{12} - \left( \frac{\sigma}{q} \right)^{6} \right] ,
```
with Hamiltonian and Lagrangian
```math
H(q, p) = \frac{p^2}{2m} + V(q) , \qquad L(q, \dot{q}) = \frac{m}{2} \dot{q}^2 - V(q) .
```
The potential has its minimum at ``q = 2^{1/6}\sigma``; small oscillations about this point are
bounded. Only ``q > 0`` is physical.

System parameters:
* `m`: mass
* `ε`: depth of the potential well
* `σ`: finite distance at which the potential is zero
"""
module LennardJonesOscillator

    using EulerLagrange
    using LinearAlgebra
    using Parameters
    using GeometricEquations: HODEEnsemble, LODEEnsemble

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hodeensemble, lodeensemble
    export hamiltonian_system, lagrangian_system

    const DEFAULT_TIMESPAN = (0.0, 10.0)
    const DEFAULT_TIMESTEP = 0.01

    default_parameters(::Type{T}=Float64) where {T} = (m = T(1.0), ε = T(1.0), σ = T(1.0))

    # start near the potential minimum q = 2^(1/6) σ, at rest
    const q₀ = [1.2]
    const p₀ = [0.0]

    potential(q, ε, σ) = 4 * ε * ((σ / q)^12 - (σ / q)^6)

    function hamiltonian(t, q, p, parameters)
        @unpack m, ε, σ = parameters
        p[1]^2 / (2m) + potential(q[1], ε, σ)
    end

    function lagrangian(t, q, q̇, parameters)
        @unpack m, ε, σ = parameters
        m * q̇[1]^2 / 2 - potential(q[1], ε, σ)
    end

    function v̄(v, t, q, p, parameters)
        v[1] = p[1] / parameters.m
        nothing
    end

    function hamiltonian_system(parameters::NamedTuple)
        t, q, p = hamiltonian_variables(1)
        sparams = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(1)
        sparams = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
    end

    # Build the symbolic system from a single parameter set, while a vector of parameter
    # sets is passed on to the ensemble unchanged (see issue #64).
    _parameters(p::NamedTuple) = p
    _parameters(p::AbstractVector) = p[begin]

    "Hamiltonian problem for the Lennard-Jones oscillator."
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Hamiltonian ensemble for the Lennard-Jones oscillator (varying initial conditions and/or parameters)."
    function hodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        eqs = functions(hamiltonian_system(_parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Lagrangian problem for the Lennard-Jones oscillator."
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

    "Lagrangian ensemble for the Lennard-Jones oscillator (varying initial conditions and/or parameters)."
    function lodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        eqs = functions(lagrangian_system(_parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
