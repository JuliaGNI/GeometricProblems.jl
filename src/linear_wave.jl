@doc raw"""
The discretized version of the 1d linear wave equation with homogeneous Dirichlet boundary
conditions.

It is a prime example of a non-trivial completely integrable system.

The system is discretized on `N` equidistant interior points ``ξ_i = iΔξ - 1/2``,
``i = 1, …, N``, with ``Δξ = 1/(N+1)``; `N` is the ``\tilde{N}`` of the documentation page. The two
boundary values ``q_0 = q_{N+1} = 0`` are fixed and are not part of the state, which therefore has
exactly `N` components. This is the full-order model of [buchfink2023symplectic](@cite): the
Hamiltonian is ``H = \tfrac12 p^Tp - \tfrac{μ^2}{2} q^T D_{ξξ} q`` with the central-difference
matrix ``D_{ξξ}`` and the canonical Poisson matrix.

Like the number of lattice sites of `TodaLattice`, `N` fixes the *size* of the system: it sets the
number of degrees of freedom and the summation bounds, so it cannot survive `symbolize` and hence
cannot be a system parameter. It is instead a leading positional argument of the problem
constructors, defaulting to `Ñ = 256`, and `hamiltonian`/`lagrangian` take it as a trailing
argument. The only system parameter proper is ``\mu``.

`hodeproblem`, `lodeproblem`, `hodeensemble` and `lodeensemble` use hand-written vector fields by
default: in-place `v`, `f`, `ϑ`, `g` and `ω!` over a shared [`∇V!`](@ref) kernel. Passing
`symbolic = true` generates the equations of motion with EulerLagrange via [`hamiltonian_system`](@ref)
or [`lagrangian_system`](@ref) instead, which agrees to round-off and is kept for cross-checking.

The symbolic route is unusable at the default size. EulerLagrange builds the ``2n × 2n`` Lagrange
two-form by symbolically differentiating a *dense* matrix, which at ``n = 258`` is 266 k entries:
`lagrangian_system(256, …)` takes 155 s and emits 14 MB of code for an `ω` whose first evaluation
costs a further 147 s — all for a two-form that no integrator in GeometricIntegrators ever evaluates
and that `GeometricEquations.check_methods` skips. Construction grows as ``n^{2.5}`` and the
generated `ω` as ``n^2``, while `hodeproblem`/`lodeproblem` in their default form build in a tenth of
a millisecond at any size. Per call the hand-written code is the faster one too, by 1.4–1.9× for the
force and 6–8× for `H` and `L`. See `benchmark/linear_wave.jl`. (These figures were measured
before 0.10, with the two boundary points still in the state; the scaling is unchanged.)
"""
module LinearWave

using EulerLagrange
using LinearAlgebra
using Parameters
using GeometricEquations: HODEEnsemble, LODEEnsemble

export hamiltonian, lagrangian
export hodeproblem, lodeproblem
export hodeensemble, lodeensemble
export hamiltonian_system, lagrangian_system

include("bump_initial_condition.jl")

const μ̃ = 0.6
const Ñ = 256

default_parameters(::Type{T} = Float64) where {T} = (μ = T(μ̃),)

"""
    interior_domain(N, T = Float64)

The `N` interior grid points ``ξ_i = i/(N+1) - 1/2``, ``i = 1, …, N``, at which the state is
defined. `compute_domain(N + 2)` is the same grid with the two boundary points included.
"""
interior_domain(N::Integer, T = Float64) = compute_domain(N + 2, T)[(begin + 1):(end - 1)]

@doc raw"""
    initial_condition(μ, N)

The initial condition ``q(0, ξ_i) = u_0(ξ_i; μ)``, ``p(0, ξ_i) = -μ ∂_ξ u_0(ξ_i; μ)`` on the `N`
[`interior_domain`](@ref) points, as a named tuple `(q, p)`. Before it reaches a boundary, the
solution of the continuous problem from this initial condition is the travelling wave
``u_0(ξ - μt; μ)``.
"""
function initial_condition(μ::T, N::Integer) where {T}
    Ω = interior_domain(N, T)
    (q = u₀(Ω, μ), p = compute_p₀(Ω, μ))
end

@doc raw"""
    potential(q, parameters, N)

The discretized gradient energy on `N` interior points with ``q_0 = q_{N+1} = 0``,

```math
V(q) = \frac{\mu^2}{2Δξ^2} \sum_{i = 0}^{N} (q_{i+1} - q_i)^2
     = -\frac{\mu^2}{2} q^T D_{ξξ} q,
\qquad Δξ = \frac{1}{N+1},
```

where ``D_{ξξ}`` is the ``N × N`` central-difference matrix ``\mathrm{tridiag}(1, -2, 1)/Δξ^2``.

Both the Hamiltonian and the Lagrangian take `N` as a trailing argument (as in `TodaLattice`)
rather than reading it from `parameters`: it fixes the number of degrees of freedom and the
summation bounds, so it has to be a plain integer and cannot survive `symbolize`.
"""
function potential(q, parameters, N)
    @unpack μ = parameters

    Δx² = (one(μ) / (N + 1)) ^ 2
    μ ^ 2 / 2Δx² * (q[1] ^ 2 + sum((q[i + 1] - q[i]) ^ 2 for i in 1:(N - 1); init = zero(μ)) +
                    q[N] ^ 2)
end

@doc raw"""
    ∇V!(dV, q, parameters, N, α = 1)

``α`` times the gradient ``∂V/∂q = -μ^2 D_{ξξ} q`` of the discretized gradient energy
[`potential`](@ref), written into `dV`. `α = -1` gives the force ``μ^2 D_{ξξ} q`` directly, which
is what the force functions want and saves a second pass over the state compared with negating
afterwards.

The boundary values ``q_0 = q_{N+1} = 0`` enter the stencil at the first and the last point, so the
two end rows drop one neighbor and the loop over the interior rows is branch-free.
`test/linear_wave.jl` checks this against a central finite difference of [`potential`](@ref) at
several `N`, including `N = 1`, where the two end rows coincide.
"""
function ∇V!(dV, q, parameters, N, α = 1)
    @unpack μ = parameters

    Δx² = (one(μ) / (N + 1)) ^ 2
    c = α * μ ^ 2 / Δx²

    @inbounds begin
        if N == 1
            dV[1] = 2c * q[1]
        else
            dV[1] = c * (2q[1] - q[2])
            for j in 2:(N - 1)
                dV[j] = c * (2q[j] - q[j - 1] - q[j + 1])
            end
            dV[N] = c * (2q[N] - q[N - 1])
        end
    end
    nothing
end

function hamiltonian(t, q, p, parameters, N)
    sum(p[n] ^ 2 for n in 1:N) / 2 + potential(q, parameters, N)
end

function lagrangian(t, q, q̇, parameters, N)
    sum(q̇[n] ^ 2 for n in 1:N) / 2 - potential(q, parameters, N)
end

# GeometricEquations requires a HODE's `H` and a LODE's `l` to be callable as (t, q, p, params);
# see `GeometricEquations.check_methods`. The state has N components, so the lattice size can be
# read off it and need not be closed over. `hamiltonian_system`/`lagrangian_system` keep calling the
# five-argument methods, so that the summation bounds stay a plain integer, as `symbolize` requires.
hamiltonian(t, q, p, parameters) = hamiltonian(t, q, p, parameters, length(q))
lagrangian(t, q, q̇, parameters) = lagrangian(t, q, q̇, parameters, length(q))

_timestep(timespan::Tuple, n_time_steps::Integer) = (timespan[2] - timespan[1]) /
                                                    (n_time_steps-1)

const DEFAULT_TIMESPAN = (0, 1)
const n_time_steps = 200
const DEFAULT_TIMESTEP = _timestep(DEFAULT_TIMESPAN, n_time_steps)

const q₀ = initial_condition(μ̃, Ñ).q
const p₀ = initial_condition(μ̃, Ñ).p

# ---------------------------------------------------------------------------------------------
# Hand-written vector fields
#
# These are what the four problem and ensemble constructors use by default. They follow the
# pattern of the other hand-written problems in this package (see `outer_solar_system.jl` and
# `harmonic_oscillator.jl`), including the extra arity methods GeometricEquations calls.
#
# None of them takes `N`: the state has N components, so N = length(q). That makes them
# size-agnostic, which is the whole point — the symbolic route has to be rebuilt per size.
# ---------------------------------------------------------------------------------------------

# In-place velocity estimate from momentum. The mass matrix is the identity, so q̇ = p — the same
# expression as ∂H/∂p, hence also the initial guess the LODE wants.
function v̄(v, t, q, p, params)
    v .= p
    nothing
end

# ∂H/∂p = p
linear_wave_v(v, t, q, p, params) = v̄(v, t, q, p, params)

# -∂H/∂q = -∂V/∂q, and ∂L/∂q = -∂V/∂q as well, so the Hamiltonian and Lagrangian force functions
# coincide and differ only in whether the third slot holds p or q̇.
function linear_wave_f(f, t, q, w, params)
    ∇V!(f, q, params, length(q), -1)
    nothing
end

# ϑ = ∂L/∂q̇ = q̇
function linear_wave_ϑ(Θ, t, q, w, params)
    Θ .= w
    nothing
end

# Projection direction for enforcing p = ϑ(q, q̇). The Lagrangian is regular, so the constraint is
# non-degenerate and the projection acts in every direction; as in `harmonic_oscillator.jl` and
# `outer_solar_system.jl` the multiplier absorbs the (here identity) mass matrix, giving g = λ.
function linear_wave_g(g, t, q, λ, params)
    g .= λ
    nothing
end

linear_wave_g(g, t, q, w, λ, params) = linear_wave_g(g, t, q, λ, params)

@doc raw"""
    ω!(Ω, t, q, params)

The Lagrange two-form ``ω = dθ`` of this regular Lagrangian, where ``θ = ϑ_i \, dq^i``.

A regular Lagrangian is a second-order system of ``n`` equations, equivalent to a first-order
system of ``2n``, so `Ω` is the ``2n × 2n`` form on ``(q, \dot q)`` with ``n = N``. In block
form it is ``[∂ϑ_i/∂q_j - ∂ϑ_j/∂q_i \; -M^T; \; M \; 0]``; since ``ϑ = ∂L/∂\dot q = \dot q``
depends on the velocities alone and the mass matrix is the identity, the upper-left block
vanishes and what remains is the canonical ``[0 \; -I; \; I \; 0]``.

This is the convention EulerLagrange's `LagrangianSystem` produces, so `lodeproblem` agrees
whether or not `symbolic = true`. Note that no integrator in GeometricIntegrators evaluates `ω`
and `GeometricEquations.check_methods` skips it — which is what makes the 14 MB of code
EulerLagrange emits for it at the default size, and the 147 s its first evaluation then costs,
pure overhead.
"""
function ω!(Ω, t, q, params)
    n = length(q)
    fill!(Ω, zero(eltype(Ω)))
    @inbounds for i in 1:n
        Ω[i, n + i] = -one(eltype(Ω))
        Ω[n + i, i] = +one(eltype(Ω))
    end
    nothing
end

# LODE/LDAE evaluate the two-form with an extra velocity slot; ω depends on neither q nor q̇.
ω!(Ω, t, q, w, params) = ω!(Ω, t, q, params)

# ---------------------------------------------------------------------------------------------
# Symbolic formulation (EulerLagrange)
#
# Kept for cross-checking the hand-written functions above, and reachable via `symbolic = true`.
# `hamiltonian_system` is cheap (1.9 s at the default size); `lagrangian_system` is not, because
# EulerLagrange differentiates a dense (2n)×(2n) matrix for `ω` and runs `build_function` over
# all 266 k entries of it.
# ---------------------------------------------------------------------------------------------

"""
    hamiltonian_system(N, parameters)

The EulerLagrange `HamiltonianSystem` for the linear wave equation on `N` interior points, from
which `hodeproblem(…; symbolic = true)` takes its `v`, `f` and `H`.

Cheap at every size — 1.9 s at the default `Ñ = 256` — since a Hamiltonian system carries no
two-form. Contrast [`lagrangian_system`](@ref).
"""
function hamiltonian_system(N::Int, parameters::NamedTuple)
    t, q, p = hamiltonian_variables(N)
    sparams = symbolize(parameters)
    HamiltonianSystem(hamiltonian(t, q, p, sparams, N), t, q, p, sparams; nanmath = true)
end

"""
    lagrangian_system(N, parameters)

The EulerLagrange `LagrangianSystem` for the linear wave equation on `N` interior points, from
which `lodeproblem(…; symbolic = true)` takes its `ϑ`, `f`, `g`, `ω` and `L`.

This is the expensive one: it builds the ``2n × 2n`` two-form by differentiating a dense matrix,
so its cost grows as ``n^{2.5}`` and reaches 155 s at the default `Ñ = 256`. See the module
docstring and `benchmark/linear_wave.jl`.
"""
function lagrangian_system(N::Int, parameters::NamedTuple)
    t, x, v = lagrangian_variables(N)
    sparams = symbolize(parameters)
    LagrangianSystem(lagrangian(t, x, v, sparams, N), t, x, v, sparams; nanmath = true)
end

# Build the symbolic system from a single parameter set, while a vector of parameter
# sets is passed on to the ensemble unchanged (see issue #64).
_parameters(p::NamedTuple) = p
_parameters(p::AbstractVector) = p[begin]

# The state has N components, one per interior point, so N is the length of an initial condition.
_nint(q::AbstractArray{<:Number}) = length(q)
_nint(q::AbstractVector{<:AbstractArray}) = length(q[begin])

# The hand-written vector fields read the lattice size off the state while the symbolic ones have
# it baked in, so an `N` that disagrees with the initial condition would silently mean two
# different systems depending on `symbolic`. Reject it up front.
#
# *Every* sample of an ensemble is checked, not just the first one `_nint` reports on: the
# hand-written fields read the size off each state individually, so a ragged ensemble would
# integrate its samples as lattices of different sizes, while the generated ones have one size
# baked in and would leave the components past `N` of an oversized sample without a force.
# Neither errors on its own, and neither is what was asked for.
function _check_size(N, q₀, p₀)
    _check_components(N, q₀, "q₀")
    _check_components(N, p₀, "p₀")
    nothing
end

_check_components(N, x::AbstractArray{<:Number}, name) = @assert length(x) == N "$name has $(length(x)) components, expected N = $N"

function _check_components(N, x::AbstractVector{<:AbstractArray}, name)
    for (i, sample) in pairs(x)
        @assert length(sample) == N "$name sample $i has $(length(sample)) components, expected N = $N"
    end
    nothing
end

"""
    hodeproblem(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Hamiltonian problem for the linear wave equation on `N` interior points.

Constructor with default arguments:
```
hodeproblem(
    N  = $(Ñ),
    q₀ = initial_condition(μ̃, N).q,
    p₀ = initial_condition(μ̃, N).p;
    timespan   = $(DEFAULT_TIMESPAN),
    timestep   = $(DEFAULT_TIMESTEP),
    parameters = $(default_parameters()),
    symbolic   = false
)
```

With `symbolic = true` the equations of motion are generated with EulerLagrange via
[`hamiltonian_system`](@ref) instead of using the hand-written vector fields. The two agree to
round-off; the symbolic route is kept for cross-checking and costs 1.9 s to build at `N = 256`.
"""
function hodeproblem(N::Int = Ñ, q₀ = initial_condition(μ̃, N).q,
        p₀ = initial_condition(μ̃, N).p;
        timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP,
        parameters = default_parameters(),
        symbolic = false)
    _check_size(N, q₀, p₀)
    if symbolic
        HODEProblem(hamiltonian_system(N, parameters), timespan, timestep, q₀, p₀;
            parameters = parameters)
    else
        HODEProblem(linear_wave_v, linear_wave_f, hamiltonian,
            timespan, timestep, q₀, p₀; parameters = parameters)
    end
end

function hodeproblem(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    hodeproblem(_nint(q₀), q₀, p₀; kwargs...)
end

"""
    hodeensemble(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Hamiltonian ensemble for the linear wave equation (varying initial conditions and/or parameters).

Takes the same arguments as [`hodeproblem`](@ref), with `q₀`, `p₀` and/or `parameters` given as
vectors of samples.
"""
function hodeensemble(N::Int = Ñ, q₀ = initial_condition(μ̃, N).q,
        p₀ = initial_condition(μ̃, N).p;
        timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP,
        parameters = default_parameters(),
        symbolic = false)
    _check_size(N, q₀, p₀)
    if symbolic
        eqs = functions(hamiltonian_system(N, _parameters(parameters)))
        HODEEnsemble(
            eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    else
        HODEEnsemble(linear_wave_v, linear_wave_f, hamiltonian,
            timespan, timestep, q₀, p₀; parameters = parameters)
    end
end

function hodeensemble(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    hodeensemble(_nint(q₀), q₀, p₀; kwargs...)
end

"""
    lodeproblem(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Lagrangian problem for the linear wave equation on `N` interior points.

Constructor with default arguments:
```
lodeproblem(
    N  = $(Ñ),
    q₀ = initial_condition(μ̃, N).q,
    p₀ = initial_condition(μ̃, N).p;
    timespan   = $(DEFAULT_TIMESPAN),
    timestep   = $(DEFAULT_TIMESTEP),
    parameters = $(default_parameters()),
    symbolic   = false
)
```

With `symbolic = true` the equations of motion are generated with EulerLagrange via
[`lagrangian_system`](@ref) instead of using the hand-written vector fields. The two agree to
round-off, but the symbolic route takes 155 s to build at `N = 256` — see the module docstring.
"""
function lodeproblem(N::Int = Ñ, q₀ = initial_condition(μ̃, N).q,
        p₀ = initial_condition(μ̃, N).p;
        timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP,
        parameters = default_parameters(),
        symbolic = false)
    _check_size(N, q₀, p₀)
    if symbolic
        LODEProblem(lagrangian_system(N, parameters), timespan, timestep, q₀, p₀;
            v̄ = v̄, parameters = parameters)
    else
        LODEProblem(linear_wave_ϑ, linear_wave_f, linear_wave_g, ω!, lagrangian,
            timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end
end

function lodeproblem(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    lodeproblem(_nint(q₀), q₀, p₀; kwargs...)
end

"""
    lodeensemble(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Lagrangian ensemble for the linear wave equation (varying initial conditions and/or parameters).

Takes the same arguments as [`lodeproblem`](@ref), with `q₀`, `p₀` and/or `parameters` given as
vectors of samples.
"""
function lodeensemble(N::Int = Ñ, q₀ = initial_condition(μ̃, N).q,
        p₀ = initial_condition(μ̃, N).p;
        timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP,
        parameters = default_parameters(),
        symbolic = false)
    _check_size(N, q₀, p₀)
    if symbolic
        eqs = functions(lagrangian_system(N, _parameters(parameters)))
        LODEEnsemble(eqs.ϑ, eqs.f, eqs.g, eqs.ω, eqs.L, timespan, timestep, q₀, p₀;
            v̄ = v̄, parameters = parameters)
    else
        LODEEnsemble(linear_wave_ϑ, linear_wave_f, linear_wave_g, ω!, lagrangian,
            timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end
end

function lodeensemble(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    lodeensemble(_nint(q₀), q₀, p₀; kwargs...)
end

end
