@doc raw"""
# Duffing Oscillator

The (undamped, unforced) Duffing oscillator is a one-dimensional anharmonic oscillator with a
quartic potential. In this conservative form it is an autonomous Hamiltonian system with

```math
H(q, p) = \frac{p^2}{2m} + \frac{\alpha}{2} q^2 + \frac{\beta}{4} q^4 ,
```
and Lagrangian
```math
L(q, \dot{q}) = \frac{m}{2} \dot{q}^2 - \frac{\alpha}{2} q^2 - \frac{\beta}{4} q^4 ,
```
giving the equation of motion ``m \ddot{q} + \alpha q + \beta q^3 = 0``.

System parameters:
* `m`: mass
* `α`: linear stiffness
* `β`: coefficient of the cubic force (quartic potential) nonlinearity
"""
module DuffingOscillator

    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const timespan = (0.0, 100.0)
    const timestep = 0.1

    const default_parameters = (m = 1.0, α = 1.0, β = 1.0)

    const q₀ = [1.0]
    const p₀ = [0.0]

    function hamiltonian(t, q, p, parameters)
        @unpack m, α, β = parameters
        p[1]^2 / (2m) + α * q[1]^2 / 2 + β * q[1]^4 / 4
    end

    function lagrangian(t, q, q̇, parameters)
        @unpack m, α, β = parameters
        m * q̇[1]^2 / 2 - α * q[1]^2 / 2 - β * q[1]^4 / 4
    end

    function v̄(v, t, q, p, parameters)
        v[1] = p[1] / parameters.m
        nothing
    end

    "Hamiltonian problem for the Duffing oscillator."
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(1)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Lagrangian problem for the Duffing oscillator."
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, x, v = lagrangian_variables(1)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
