@doc raw"""
# Mathews-Lakshmanan Oscillator

The Mathews-Lakshmanan oscillator is a nonlinear oscillator with a position-dependent effective
mass, defined by the Lagrangian
```math
L(q, \dot{q}) = \frac{1}{2} \frac{\dot{q}^2 - \omega^2 q^2}{1 + \lambda q^2} .
```
The conjugate momentum is ``p = \dot{q} / (1 + \lambda q^2)``, and the corresponding Hamiltonian is
```math
H(q, p) = \frac{1}{2} \big( 1 + \lambda q^2 \big) p^2 + \frac{1}{2} \frac{\omega^2 q^2}{1 + \lambda q^2} .
```
It exhibits amplitude-dependent (anisochronous for ``\lambda \neq 0``) oscillations and is a
standard example of a quadratic Liénard-type nonlinear oscillator (M. Lakshmanan and
S. Rajasekar, *Nonlinear Dynamics*, Springer, 2003).

System parameters:
* `ω`: (bare) angular frequency
* `λ`: nonlinearity parameter
"""
module MathewsLakshmananOscillator

    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const timespan = (0.0, 20.0)
    const timestep = 0.01

    const default_parameters = (ω = 1.0, λ = 1.0)

    const q₀ = [0.5]
    const p₀ = [0.0]

    function hamiltonian(t, q, p, parameters)
        @unpack ω, λ = parameters
        (1 + λ * q[1]^2) * p[1]^2 / 2 + ω^2 * q[1]^2 / (2 * (1 + λ * q[1]^2))
    end

    function lagrangian(t, q, q̇, parameters)
        @unpack ω, λ = parameters
        (q̇[1]^2 - ω^2 * q[1]^2) / (2 * (1 + λ * q[1]^2))
    end

    function v̄(v, t, q, p, parameters)
        v[1] = p[1] * (1 + parameters.λ * q[1]^2)
        nothing
    end

    "Hamiltonian problem for the Mathews-Lakshmanan oscillator."
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(1)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Lagrangian problem for the Mathews-Lakshmanan oscillator."
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, x, v = lagrangian_variables(1)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
