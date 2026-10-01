@doc raw"""
# Kepler Problem

The planar two-body problem in relative coordinates,
```math
L(q, \dot{q}) = \frac{1}{2} |\dot{q}|^2 + \frac{\mu}{|q|} , \qquad
H(q, p) = \frac{1}{2} |p|^2 - \frac{\mu}{|q|} .
```
Besides the energy ``H``, the angular momentum ``\ell = q_1 p_2 - q_2 p_1`` and the
Runge–Lenz vector
```math
R = \ell \, (p_2, \, -p_1) - \mu \, \frac{q}{|q|}
```
are conserved. ``R`` points to the perihelion and has length ``\mu e``.

The default initial condition is the perihelion of the orbit with semi-major axis `a = 1` and
eccentricity `e = 0.9`, see [`initial_condition`](@ref). Every bound orbit (``H < 0``) with
``0 < e < 1`` has the closed-form solution [`exact_solution`](@ref), obtained from Kepler's
equation ``M = E - e \sin E``.

System parameters:
* `μ`: gravitational parameter
"""
module KeplerProblem

using EulerLagrange
using Parameters

export hamiltonian, lagrangian, angular_momentum, runge_lenz_vector
export initial_condition, exact_solution
export hodeproblem, lodeproblem
export hamiltonian_system, lagrangian_system
export default_parameters

default_parameters(::Type{T} = Float64) where {T} = (μ = T(1.0),)

"Perihelion ``q = (a (1 - e), 0)``, ``p = (0, \\sqrt{\\mu (1 + e) / (a (1 - e))})`` of the prograde orbit."
function initial_condition(e, a = 1.0; parameters = default_parameters())
    @unpack μ = parameters
    ([a * (1 - e), 0.0], [0.0, sqrt(μ * (1 + e) / (a * (1 - e)))])
end

const q₀, p₀ = initial_condition(0.9)

# ten orbital periods 2π √(a³/μ) for a = μ = 1
const DEFAULT_TIMESPAN = (0.0, 20π)
const DEFAULT_TIMESTEP = 0.01

function hamiltonian(t, q, p, params)
    @unpack μ = params
    (p[1]^2 + p[2]^2) / 2 - μ / sqrt(q[1]^2 + q[2]^2)
end

function lagrangian(t, q, q̇, params)
    @unpack μ = params
    (q̇[1]^2 + q̇[2]^2) / 2 + μ / sqrt(q[1]^2 + q[2]^2)
end

angular_momentum(t, q, p, params) = q[1] * p[2] - q[2] * p[1]

function runge_lenz_vector(t, q, p, params)
    @unpack μ = params
    ℓ = angular_momentum(t, q, p, params)
    r = sqrt(q[1]^2 + q[2]^2)
    [ℓ * p[2] - μ * q[1] / r, -ℓ * p[1] - μ * q[2] / r]
end

function v̄(v, t, q, p, params)
    v[1] = p[1]
    v[2] = p[2]
    nothing
end

@doc raw"""
    eccentric_anomaly(M, e)

Solve Kepler's equation ``M = E - e \sin E`` for the eccentric anomaly ``E``. The left-hand side
is increasing in ``E`` and the root lies in ``[M - e, M + e]``, so Newton's method is safeguarded
by bisection on that bracket; this keeps it convergent for ``e \to 1``, where ``1 - e \cos E``
nearly vanishes at the perihelion.
"""
function eccentric_anomaly(M, e)
    a, b = M - e, M + e
    E = M
    for _ in 1:100
        f = E - e * sin(E) - M
        f == 0 && return E
        f < 0 ? (a = E) : (b = E)
        Eₙ = E - f / (1 - e * cos(E))
        Eₙ = a < Eₙ < b ? Eₙ : (a + b) / 2
        abs(Eₙ - E) ≤ 2eps(abs(E) + 1) && return Eₙ
        E = Eₙ
    end
    error("Kepler's equation did not converge for M = $M, e = $e.")
end

@doc raw"""
    exact_solution(t, t₀, q₀, p₀, params)

Position and momentum at time `t` of the orbit through `(q₀, p₀)` at time `t₀`. Requires a bound
orbit with eccentricity ``0 < e < 1``.
"""
function exact_solution(t, t₀, q₀, p₀, params)
    @unpack μ = params
    r₀ = sqrt(q₀[1]^2 + q₀[2]^2)
    a = -μ / (2 * hamiltonian(t₀, q₀, p₀, params))
    R = runge_lenz_vector(t₀, q₀, p₀, params)
    e = sqrt(R[1]^2 + R[2]^2) / μ
    @assert a > 0 && 0 < e < 1 "exact_solution requires a bound orbit with 0 < e < 1."

    s = sign(angular_momentum(t₀, q₀, p₀, params))      # orientation of the orbit
    b = a * sqrt(1 - e^2)
    n = sqrt(μ / a^3)                                      # mean motion
    cosϖ, sinϖ = R[1] / (μ * e), R[2] / (μ * e)             # direction of the perihelion

    # eccentric and mean anomaly at t₀
    E₀ = atan((q₀[1] * p₀[1] + q₀[2] * p₀[2]) / (e * sqrt(μ * a)), (1 - r₀ / a) / e)
    E = eccentric_anomaly(rem2pi(E₀ - e * sin(E₀) + n * (t - t₀), RoundNearest), e)

    # position and velocity in the perifocal frame, then rotated by the perihelion angle ϖ
    x, y = a * (cos(E) - e), s * b * sin(E)
    Ė = n / (1 - e * cos(E))
    ẋ, ẏ = -a * sin(E) * Ė, s * b * cos(E) * Ė

    ([cosϖ * x - sinϖ * y, sinϖ * x + cosϖ * y], [cosϖ * ẋ - sinϖ * ẏ, sinϖ * ẋ + cosϖ * ẏ])
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

"Hamiltonian problem for the Kepler problem."
function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    HODEProblem(
        hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
end

"Lagrangian problem for the Kepler problem."
function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    LODEProblem(lagrangian_system(parameters), timespan,
        timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
end

end
