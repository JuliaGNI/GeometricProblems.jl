@doc raw"""
# Charged particle in 2D

A charged particle of mass ``m`` (and unit charge) moving in the plane, in the static magnetic
field ``B \, e_z`` and electrostatic potential ``\phi`` of
[`GeometricProblems.MasslessChargedParticle`](@ref). The Lagrangian is
```math
L(x, \dot{x}) = \frac{m}{2} |\dot{x}|^2 + A(x) \cdot \dot{x} - \phi(x) ,
```
with ``A`` and ``\phi`` as in `MasslessChargedParticle`, so that ``B(x) = A_0 (1 + 2 |x|^2)``.
The canonical momentum is ``p = m \dot{x} + A(x)`` and the Hamiltonian reads
```math
H(x, p) = \frac{|p - A(x)|^2}{2m} + \phi(x) .
```

The particle gyrates with frequency ``\Omega(x) = B(x) / m`` about a guiding centre that, for
``m \to 0``, follows the ``E \times B`` drift ``v_E(x)`` of `MasslessChargedParticle`. The
magnetic moment of the gyration, measured in the drift frame,
```math
\mu = \frac{m |\dot{x} - v_E(x)|^2}{2 B(x)} ,
```
is an adiabatic invariant for ``m \to 0``.

The default initial velocity is that drift velocity at ``q_0``, see [`initial_momentum`](@ref).

System parameters:
* `m`: mass
* `A₀`: magnetic field strength
* `E₀`: electric field strength
"""
module ChargedParticle2d

using EulerLagrange
using Parameters

using ..MasslessChargedParticle: A, B, ϕ, v₁, v₂

export hamiltonian, lagrangian, magnetic_moment, initial_momentum, drift_velocity
export hodeproblem, lodeproblem
export hamiltonian_system, lagrangian_system
export default_parameters

default_parameters(::Type{T} = Float64) where {T} = (m = T(1E-2), A₀ = T(1.0), E₀ = T(1.0))

const DEFAULT_TIMESPAN = (0.0, 10.0)
const DEFAULT_TIMESTEP = 1E-3

"Canonical momentum ``p = m v + A(q)`` of a particle at `q` with velocity `v`."
initial_momentum(q, v, params) = params.m .* v .+ A(q, params)

"``E \\times B`` drift velocity of the massless particle at `q`."
drift_velocity(q, params) = [v₁(0, q, params), v₂(0, q, params)]

const q₀ = [1.0, 1.0]

function hamiltonian(t, q, p, params)
    @unpack m = params
    a = A(q, params)
    ((p[1] - a[1])^2 + (p[2] - a[2])^2) / (2m) + ϕ(q, params)
end

function lagrangian(t, q, q̇, params)
    @unpack m = params
    a = A(q, params)
    m * (q̇[1]^2 + q̇[2]^2) / 2 + a[1] * q̇[1] + a[2] * q̇[2] - ϕ(q, params)
end

function magnetic_moment(t, q, p, params)
    @unpack m = params
    a = A(q, params)
    w₁ = (p[1] - a[1]) / m - v₁(t, q, params)
    w₂ = (p[2] - a[2]) / m - v₂(t, q, params)
    m * (w₁^2 + w₂^2) / (2 * B(q, params))
end

function v̄(v, t, q, p, params)
    @unpack m = params
    a = A(q, params)
    v[1] = (p[1] - a[1]) / m
    v[2] = (p[2] - a[2]) / m
    nothing
end

function hamiltonian_system(parameters::NamedTuple)
    t, q, p = hamiltonian_variables(2)
    sparams = symbolize(parameters)
    HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams; nanmath = true)
end

function lagrangian_system(parameters::NamedTuple)
    t, x, v = lagrangian_variables(2)
    sparams = symbolize(parameters)
    LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams; nanmath = true)
end

# `p₀ = nothing` resolves to the drift velocity at `q₀` for the *given* parameters, as a
# positional default cannot reference the `parameters` keyword.
_momentum(q₀, p₀, parameters) = p₀
_momentum(q₀, ::Nothing, parameters) = initial_momentum(
    q₀, drift_velocity(q₀, parameters), parameters)

"Hamiltonian problem for the charged particle in 2D."
function hodeproblem(q₀ = q₀, p₀ = nothing; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀,
        _momentum(q₀, p₀, parameters); parameters = parameters)
end

"Lagrangian problem for the charged particle in 2D."
function lodeproblem(q₀ = q₀, p₀ = nothing; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀,
        _momentum(q₀, p₀, parameters); v̄ = v̄, parameters = parameters)
end

end
