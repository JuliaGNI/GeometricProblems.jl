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
    using GeometricEquations: HODEEnsemble, LODEEnsemble

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hodeensemble, lodeensemble
    export hamiltonian_system, lagrangian_system

    const DEFAULT_TIMESPAN = (0.0, 20.0)
    const DEFAULT_TIMESTEP = 0.01

    default_parameters(::Type{T}=Float64) where {T} = (ω = T(1.0), λ = T(1.0))

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

    "Hamiltonian problem for the Mathews-Lakshmanan oscillator."
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Hamiltonian ensemble for the Mathews-Lakshmanan oscillator (varying initial conditions and/or parameters)."
    function hodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        eqs = functions(hamiltonian_system(_parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    "Lagrangian problem for the Mathews-Lakshmanan oscillator."
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

    "Lagrangian ensemble for the Mathews-Lakshmanan oscillator (varying initial conditions and/or parameters)."
    function lodeensemble(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN, timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
        eqs = functions(lagrangian_system(_parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
