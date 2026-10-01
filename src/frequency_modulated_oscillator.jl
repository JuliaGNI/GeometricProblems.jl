@doc raw"""
# Frequency-Modulated Oscillator

A one-dimensional harmonic oscillator whose frequency varies slowly in time,
```math
L(t, q, \dot{q}) = \frac{1}{2} \dot{q}^2 - \frac{1}{2} \omega(\varepsilon t)^2 q^2 , \qquad
\omega(s) = \omega_0 \big( 1 + \delta \sin s \big) ,
```
with Hamiltonian
```math
H(t, q, p) = \frac{1}{2} p^2 + \frac{1}{2} \omega(\varepsilon t)^2 q^2 .
```

The system is non-autonomous, so ``H`` is not conserved. For ``\varepsilon \ll 1`` the action
``J = H / \omega`` is an adiabatic invariant: it stays within ``O(\varepsilon)`` of its initial
value over the time scale ``1/\varepsilon`` of the modulation.

System parameters:
* `ε`: modulation rate
* `ω₀`: mean frequency
* `δ`: relative modulation amplitude (``|\delta| < 1``)
"""
module FrequencyModulatedOscillator

using EulerLagrange
using Parameters

export hamiltonian, lagrangian, frequency, adiabatic_invariant
export hodeproblem, lodeproblem
export hamiltonian_system, lagrangian_system
export default_parameters

default_parameters(::Type{T} = Float64) where {T} = (ε = T(1E-2), ω₀ = T(1.0), δ = T(0.5))

const q₀ = [1.0]
const p₀ = [0.0]

# one period 2π/ε of the modulation for the default ε
const DEFAULT_TIMESPAN = (0.0, 2π / default_parameters().ε)
const DEFAULT_TIMESTEP = 0.1

function frequency(t, params)
    @unpack ε, ω₀, δ = params
    ω₀ * (1 + δ * sin(ε * t))
end

function hamiltonian(t, q, p, params)
    p[1]^2 / 2 + frequency(t, params)^2 * q[1]^2 / 2
end

function lagrangian(t, q, q̇, params)
    q̇[1]^2 / 2 - frequency(t, params)^2 * q[1]^2 / 2
end

adiabatic_invariant(t, q, p, params) = hamiltonian(t, q, p, params) / frequency(t, params)

function v̄(v, t, q, p, params)
    v[1] = p[1]
    nothing
end

function hamiltonian_system(parameters::NamedTuple)
    t, q, p = hamiltonian_variables(1)
    sparams = symbolize(parameters)
    HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams; nanmath = true)
end

function lagrangian_system(parameters::NamedTuple)
    t, x, v = lagrangian_variables(1)
    sparams = symbolize(parameters)
    LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams; nanmath = true)
end

"Hamiltonian problem for the frequency-modulated oscillator."
function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    HODEProblem(
        hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
end

"Lagrangian problem for the frequency-modulated oscillator."
function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    LODEProblem(lagrangian_system(parameters), timespan,
        timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
end

end
