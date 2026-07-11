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

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const timespan = (0.0, 10.0)
    const timestep = 0.01

    const default_parameters = (m = 1.0, ε = 1.0, σ = 1.0)

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

    "Hamiltonian problem for the Lennard-Jones oscillator."
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(1)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Lagrangian problem for the Lennard-Jones oscillator."
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, x, v = lagrangian_variables(1)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
