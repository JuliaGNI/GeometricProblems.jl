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
    using GeometricEquations: HODEEnsemble, LODEEnsemble

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hodeensemble, lodeensemble
    export hamiltonian_system, lagrangian_system

    const DEFAULT_TIMESPAN = (0.0, 100.0)
    const DEFAULT_TIMESTEP = 0.1

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

    "Hamiltonian problem for the Duffing oscillator."
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Hamiltonian ensemble for the Duffing oscillator (varying initial conditions and/or parameters)."
    function hodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        eqs = functions(hamiltonian_system(_parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Lagrangian problem for the Duffing oscillator."
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

    "Lagrangian ensemble for the Duffing oscillator (varying initial conditions and/or parameters)."
    function lodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters)
        eqs = functions(lagrangian_system(_parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
