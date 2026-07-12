@doc raw"""
# Morse Oscillator

The Morse oscillator models an anharmonic bond with the potential
```math
V(q) = D \left( 1 - e^{-a (q - r_0)} \right)^2 ,
```
with Hamiltonian and Lagrangian
```math
H(q, p) = \frac{p^2}{2m} + V(q) , \qquad L(q, \dot{q}) = \frac{m}{2} \dot{q}^2 - V(q) .
```
Trajectories with energy below the dissociation energy ``D`` are bound and oscillate anharmonically
about the equilibrium ``q = r_0``.

System parameters:
* `m`: reduced mass
* `D`: well depth (dissociation energy)
* `a`: controls the width of the well
* `r₀`: equilibrium position
"""
module MorseOscillator

    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const timespan = (0.0, 20.0)
    const timestep = 0.01

    const default_parameters = (m = 1.0, D = 1.0, a = 1.0, r₀ = 0.0)

    const q₀ = [0.5]
    const p₀ = [0.0]

    potential(q, D, a, r₀) = D * (1 - exp(-a * (q - r₀)))^2

    function hamiltonian(t, q, p, parameters)
        @unpack m, D, a, r₀ = parameters
        p[1]^2 / (2m) + potential(q[1], D, a, r₀)
    end

    function lagrangian(t, q, q̇, parameters)
        @unpack m, D, a, r₀ = parameters
        m * q̇[1]^2 / 2 - potential(q[1], D, a, r₀)
    end

    function v̄(v, t, q, p, parameters)
        v[1] = p[1] / parameters.m
        nothing
    end

    "Hamiltonian problem for the Morse oscillator."
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(1)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Lagrangian problem for the Morse oscillator."
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, x, v = lagrangian_variables(1)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
