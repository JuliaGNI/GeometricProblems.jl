@doc raw"""
The Toda lattice is a model for a one-dimensional crystal named after its discoverer Morikazu Toda [toda1967vibration](@cite).

It is a prime example of a non-trivial completely integrable system.

It is configured by the *number of points* ``N`` in the periodic lattice and by ``\alpha``, which
adjusts the strength of the interactions in the lattice.

Like the number of interior points of `LinearWave`, ``N`` fixes the *size* of the system: it sets
the number of degrees of freedom and the summation bounds, so it cannot survive `symbolize` and
hence cannot be a system parameter. It is instead a leading positional argument of the problem
constructors, defaulting to `Ñ = 200`, and `hamiltonian`/`lagrangian` take it as a trailing
argument. The only system parameter proper is ``\alpha``.

`hodeproblem`, `lodeproblem`, `hodeensemble` and `lodeensemble` use hand-written vector fields by
default: in-place `v`, `f`, `ϑ`, `g` and `ω!` over a shared [`∇V!`](@ref) kernel. Passing
`symbolic = true` generates the equations of motion with EulerLagrange via [`hamiltonian_system`](@ref)
or [`lagrangian_system`](@ref) instead, which agrees to round-off and is kept for cross-checking.

The symbolic route is unusable at the default size. EulerLagrange builds the ``2N × 2N`` Lagrange
two-form by symbolically differentiating a *dense* matrix, which at ``N = 200`` is 160 k entries:
`lagrangian_system(200, …)` takes 73 s and emits 8.4 MB of code whose first `ω` evaluation costs a
further 82 s — all for a two-form that no integrator in GeometricIntegrators ever evaluates and that
`GeometricEquations.check_methods` skips. Construction grows as ``N^{2.4}`` and the generated `ω` as
``N^2``, while `hodeproblem`/`lodeproblem` in their default form build in a tenth of a millisecond at
any size.

Per call the two are much closer than for the linear wave, because both force functions are dominated
by the ``N`` exponentials rather than by arithmetic: the generated force is in fact 5–12%
*faster* over ``N = 8`` to ``200``, while the hand-written `H` and `L` are 1.3–2× faster and `v`, `ϑ`
and `g` change places with size. So the 155 s of setup a `symbolic = true` LODE pays at ``N = 200``
would take 6.9 × 10⁹ force evaluations to earn back — around six orders of magnitude more than a
default run performs. See `benchmark/toda_lattice.jl`.
"""
module TodaLattice

using EulerLagrange
using LinearAlgebra
using Parameters
using GeometricEquations: HODEEnsemble, HODEProblem, LODEEnsemble, LODEProblem

export hodeproblem, hodeensemble, hamiltonian, hamiltonian_system
export lodeproblem, lodeensemble, lagrangian, lagrangian_system

include("bump_initial_condition.jl")

default_parameters(::Type{T}=Float64) where {T} = (
    α=T(0.64),
)

@doc raw"""
    potential(q, parameters, N)

The exponential nearest-neighbour interaction of the periodic Toda lattice on `N` points,

```math
V(q) = \alpha \sum_{n = 1}^{N} e^{q_n - q_{n+1}}, \qquad q_{N+1} \equiv q_1 .
```

Both the Hamiltonian and the Lagrangian take `N` as a trailing argument (as in `LinearWave`) rather
than reading it from `parameters`: it fixes the number of degrees of freedom and the summation
bounds, so it has to be a plain integer and cannot survive `symbolize`.
"""
function potential(q, params, N)
    params.α * sum(exp(q[n] - q[n%N+1]) for n in 1:N)
end

@doc raw"""
    ∇V!(dV, q, parameters, N, scale = 1)

`scale` times the gradient ``∂V/∂q`` of the periodic Toda [`potential`](@ref), written into `dV`.
`scale = -1` gives the force ``-∂V/∂q`` directly, which is what the force functions want and saves
the extra sweep over the output that negating afterwards would cost.

Writing ``E_n = e^{q_n - q_{n+1}}`` with cyclic indices, ``q_j`` enters ``V = \alpha \sum_n E_n``
with ``+1`` in the term ``n = j`` and with ``-1`` in the term ``n = j - 1``, so

```math
\frac{∂V}{∂q_j} = \alpha \, (E_j - E_{j-1}), \qquad E_0 \equiv E_N .
```

``E_N`` is both the last term and ``E_0``, so it is computed once, before anything is written, and
kept in a scalar; the loop then carries ``E_{j-1}`` in a second scalar and writes ``∂V/∂q_j`` as it
goes. So the kernel makes a single pass, allocates nothing, evaluates each exponential exactly once
— the expensive part — and never reads an entry of `dV` back, which is what makes it safe even when
`dV` aliases `q`.

The one entry that has to be peeled off is the last, and that is also the only case analysis: the
expression is otherwise correct for every ``N ≥ 1``. At ``N = 1`` the loop body never runs and
``E_1 = e^0 = 1`` leaves ``∂V/∂q_1 = c \, (E_1 - E_1) = 0``, which is right, since ``V`` is then
constant.

The cyclic wrap-around is the one thing here that is easy to get wrong and that nothing else checks,
so `test/toda_lattice_tests.jl` verifies this against a central finite difference of
[`potential`](@ref), which shares no code with it, at several `N` including `N = 1` and `N = 2` — and
at each of them with `dV` aliased onto `q`, which is the property the single pass buys.
"""
function ∇V!(dV, q, parameters, N, scale = 1)
    @unpack α = parameters

    c = scale * α

    @inbounds begin
        # the periodic wrap-around term, which is both E_N and E₀, computed before the first write so
        # that q[1] is still intact even if dV aliases q
        Eₙ = exp(q[N] - q[1])

        # ∂V/∂q_j = α (E_j - E_{j-1}), carrying E_{j-1} in a scalar
        Eₚ = Eₙ                         # E₀ ≡ E_N
        for j in 1 : (N - 1)
            Eⱼ = exp(q[j] - q[j + 1])
            dV[j] = c * (Eⱼ - Eₚ)
            Eₚ = Eⱼ
        end
        dV[N] = c * (Eₙ - Eₚ)           # Eₚ is E_{N-1} by now
    end
    nothing
end

hamiltonian(t, q, p, params, N) = p ⋅ p / 2 + potential(q, params, N)
lagrangian(t, q, q̇, params, N) = q̇ ⋅ q̇ / 2 - potential(q, params, N)

# GeometricEquations requires a HODE's `H` and a LODE's `l` to be callable as (t, q, p, params);
# see `GeometricEquations.check_methods`. The state has exactly N components, so the lattice size can
# be read off it and need not be closed over. `hamiltonian_system`/`lagrangian_system` keep calling
# the five-argument methods, so that the summation bounds stay a plain integer, as `symbolize`
# requires.
hamiltonian(t, q, p, params) = hamiltonian(t, q, p, params, length(q))
lagrangian(t, q, q̇, params) = lagrangian(t, q, q̇, params, length(q))

const DEFAULT_TIMESTEP = 0.1
const DEFAULT_TIMESPAN = (0.0, 120.0)

# parameter for the default initial conditions
const Ñ = 200
const μ = 0.3

const q₀ = compute_initial_q(μ, Ñ)
const p₀ = zero(q₀)


# -------------------------------------------------------------------------------------------------
# Hand-written vector fields
#
# These are what the four problem and ensemble constructors use by default. They follow the pattern
# of the other hand-written problems in this package (see `linear_wave.jl` and
# `outer_solar_system.jl`), including the extra arity methods GeometricEquations calls.
#
# None of them takes `N`: the state has exactly N components, so N = length(q). That makes them
# size-agnostic, which is the whole point — the symbolic route has to be rebuilt per size.
# -------------------------------------------------------------------------------------------------

# In-place velocity estimate from momentum. The mass matrix is the identity, so q̇ = p — the same
# expression as ∂H/∂p, hence also the initial guess the LODE wants. This replaces the second, full
# `HamiltonianSystem` the LODE constructors used to build purely to obtain this one function.
function v̄(v, t, q, p, params)
    v .= p
    nothing
end

# ∂H/∂p = p
toda_lattice_v(v, t, q, p, params) = v̄(v, t, q, p, params)

# -∂H/∂q = -∂V/∂q, and ∂L/∂q = -∂V/∂q as well, so the Hamiltonian and Lagrangian force functions
# coincide and differ only in whether the third slot holds p or q̇.
function toda_lattice_f(f, t, q, w, params)
    ∇V!(f, q, params, length(q), -1)
    nothing
end

# ϑ = ∂L/∂q̇ = q̇
function toda_lattice_ϑ(Θ, t, q, w, params)
    Θ .= w
    nothing
end

# Projection direction for enforcing p = ϑ(q, q̇). The Lagrangian is regular, so the constraint is
# non-degenerate and the projection acts in every direction; as in `linear_wave.jl` and
# `outer_solar_system.jl` the multiplier absorbs the (here identity) mass matrix, giving g = λ.
function toda_lattice_g(g, t, q, λ, params)
    g .= λ
    nothing
end

toda_lattice_g(g, t, q, w, λ, params) = toda_lattice_g(g, t, q, λ, params)

@doc raw"""
    ω!(Ω, t, q, params)

The Lagrange two-form ``ω = dθ`` of this regular Lagrangian, where ``θ = ϑ_i \, dq^i``.

A regular Lagrangian is a second-order system of ``N`` equations, equivalent to a first-order system
of ``2N``, so `Ω` is the ``2N × 2N`` form on ``(q, \dot q)``. In block form it is
``[∂ϑ_i/∂q_j - ∂ϑ_j/∂q_i \; -M^T; \; M \; 0]``; since ``ϑ = ∂L/∂\dot q = \dot q`` depends on the
velocities alone and the mass matrix is the identity, the upper-left block vanishes and what remains
is the canonical ``[0 \; -I; \; I \; 0]``.

This is the convention EulerLagrange's `LagrangianSystem` produces, so `lodeproblem` agrees whether
or not `symbolic = true`. Note that no integrator in GeometricIntegrators evaluates `ω` and
`GeometricEquations.check_methods` skips it — which is what makes the 8.4 MB of code EulerLagrange
emits for it at the default size, and the 82 s its first evaluation then costs, pure overhead.
"""
function ω!(Ω, t, q, params)
    n = length(q)
    fill!(Ω, zero(eltype(Ω)))
    @inbounds for i in 1 : n
        Ω[i, n + i] = -one(eltype(Ω))
        Ω[n + i, i] = +one(eltype(Ω))
    end
    nothing
end

# LODE/LDAE evaluate the two-form with an extra velocity slot; ω depends on neither q nor q̇.
ω!(Ω, t, q, w, params) = ω!(Ω, t, q, params)


# -------------------------------------------------------------------------------------------------
# Symbolic formulation (EulerLagrange)
#
# Kept for cross-checking the hand-written functions above, and reachable via `symbolic = true`.
# `hamiltonian_system` is comparatively cheap; `lagrangian_system` is not, because EulerLagrange
# differentiates a dense (2N)×(2N) matrix for `ω` and runs `build_function` over all 160 k entries of
# it at the default size.
# -------------------------------------------------------------------------------------------------

"""
    hamiltonian_system(N, parameters)

The EulerLagrange `HamiltonianSystem` for the Toda lattice on `N` points, from which
`hodeproblem(…; symbolic = true)` takes its `v`, `f` and `H`.

Far cheaper than [`lagrangian_system`](@ref) at every size, since a Hamiltonian system carries no
two-form — 1.2 s at the default `Ñ = 200`, against a tenth of a millisecond for the hand-written
route.
"""
function hamiltonian_system(N::Int, parameters::NamedTuple)
    t, q, p = hamiltonian_variables(N)
    sparams = symbolize(parameters)
    HamiltonianSystem(hamiltonian(t, q, p, sparams, N), t, q, p, sparams; nanmath = true)
end

"""
    lagrangian_system(N, parameters)

The EulerLagrange `LagrangianSystem` for the Toda lattice on `N` points, from which
`lodeproblem(…; symbolic = true)` takes its `ϑ`, `f`, `g`, `ω` and `L`.

This is the expensive one: it builds the ``2N × 2N`` two-form by differentiating a dense matrix, so
its cost grows as ``N^{2.4}`` and reaches 73 s at the default `Ñ = 200` — followed by another 82 s
the first time the generated `ω` is evaluated. See the module docstring and
`benchmark/toda_lattice.jl`.
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

# The state has exactly N components, so an initial condition of length n corresponds to N = n. The
# ensemble constructors are handed a vector of samples instead of a single one.
_length(q::AbstractArray{<:Number}) = length(q)
_length(q::AbstractVector{<:AbstractArray}) = length(q[begin])

# The hand-written vector fields read the lattice size off the state while the symbolic ones have it
# baked in, so an `N` that disagrees with the initial condition would silently mean two different
# systems depending on `symbolic`. Reject it up front.
#
# *Every* sample of an ensemble is checked, not just the first one `_length` reports on: the
# hand-written fields read the size off each state individually, so a ragged ensemble would integrate
# its samples as lattices of different sizes, while the generated ones have one size baked in and
# would leave the components past `N` of an oversized sample without a force. Neither errors on its
# own, and neither is what was asked for.
function _check_size(N, q₀, p₀)
    _check_components(N, q₀, "q₀")
    _check_components(N, p₀, "p₀")
    nothing
end

_check_components(N, x::AbstractArray{<:Number}, name) =
    @assert length(x) == N "$name has $(length(x)) components, expected N = $N"

function _check_components(N, x::AbstractVector{<:AbstractArray}, name)
    for (i, sample) in pairs(x)
        @assert length(sample) == N "$name sample $i has $(length(sample)) components, expected N = $N"
    end
    nothing
end

"""
    hodeproblem(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Hamiltonian problem for the Toda lattice on `N` points.

Constructor with default arguments:
```
hodeproblem(
    N  = $(Ñ),
    q₀ = compute_initial_q(μ, N),
    p₀ = zero(q₀);
    timespan   = $(DEFAULT_TIMESPAN),
    timestep   = $(DEFAULT_TIMESTEP),
    parameters = $(default_parameters()),
    symbolic   = false
)
```

With `symbolic = true` the equations of motion are generated with EulerLagrange via
[`hamiltonian_system`](@ref) instead of using the hand-written vector fields. The two agree to
round-off; the symbolic route is kept for cross-checking and costs 1.2 s to build at `N = 200`.
"""
function hodeproblem(N::Int=Ñ, q₀=compute_initial_q(μ, N), p₀=zero(q₀);
                     timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP,
                     parameters=default_parameters(), symbolic=false)
    _check_size(N, q₀, p₀)
    if symbolic
        HODEProblem(hamiltonian_system(N, _parameters(parameters)), timespan, timestep, q₀, p₀;
                    parameters=parameters)
    else
        HODEProblem(toda_lattice_v, toda_lattice_f, hamiltonian,
                    timespan, timestep, q₀, p₀; parameters=parameters)
    end
end

function hodeproblem(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    hodeproblem(_length(q₀), q₀, p₀; kwargs...)
end

"""
    hodeensemble(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Hamiltonian ensemble for the Toda lattice (varying initial conditions and/or parameters).

Takes the same arguments as [`hodeproblem`](@ref), with `q₀`, `p₀` and/or `parameters` given as
vectors of samples.
"""
function hodeensemble(N::Int=Ñ, q₀=compute_initial_q(μ, N), p₀=zero(q₀);
                      timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP,
                      parameters=default_parameters(), symbolic=false)
    _check_size(N, q₀, p₀)
    if symbolic
        eqs = functions(hamiltonian_system(N, _parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters=parameters)
    else
        HODEEnsemble(toda_lattice_v, toda_lattice_f, hamiltonian,
                     timespan, timestep, q₀, p₀; parameters=parameters)
    end
end

function hodeensemble(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    hodeensemble(_length(q₀), q₀, p₀; kwargs...)
end

"""
    lodeproblem(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Lagrangian problem for the Toda lattice on `N` points.

Constructor with default arguments:
```
lodeproblem(
    N  = $(Ñ),
    q₀ = compute_initial_q(μ, N),
    p₀ = zero(q₀);
    timespan   = $(DEFAULT_TIMESPAN),
    timestep   = $(DEFAULT_TIMESTEP),
    parameters = $(default_parameters()),
    symbolic   = false
)
```

With `symbolic = true` the equations of motion are generated with EulerLagrange via
[`lagrangian_system`](@ref) instead of using the hand-written vector fields. The two agree to
round-off, but the symbolic route takes 73 s to build at `N = 200`, and 155 s in all before the
first step — see the module docstring.
"""
function lodeproblem(N::Int=Ñ, q₀=compute_initial_q(μ, N), p₀=zero(q₀);
                     timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP,
                     parameters=default_parameters(), symbolic=false)
    _check_size(N, q₀, p₀)
    if symbolic
        LODEProblem(lagrangian_system(N, _parameters(parameters)), timespan, timestep, q₀, p₀;
                    v̄=v̄, parameters=parameters)
    else
        LODEProblem(toda_lattice_ϑ, toda_lattice_f, toda_lattice_g, ω!, lagrangian,
                    timespan, timestep, q₀, p₀; v̄=v̄, parameters=parameters)
    end
end

function lodeproblem(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    lodeproblem(_length(q₀), q₀, p₀; kwargs...)
end

"""
    lodeensemble(N, q₀, p₀; timespan, timestep, parameters, symbolic)

Lagrangian ensemble for the Toda lattice (varying initial conditions and/or parameters).

Takes the same arguments as [`lodeproblem`](@ref), with `q₀`, `p₀` and/or `parameters` given as
vectors of samples.
"""
function lodeensemble(N::Int=Ñ, q₀=compute_initial_q(μ, N), p₀=zero(q₀);
                      timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP,
                      parameters=default_parameters(), symbolic=false)
    _check_size(N, q₀, p₀)
    if symbolic
        leqs = functions(lagrangian_system(N, _parameters(parameters)))
        LODEEnsemble(leqs.ϑ, leqs.f, leqs.g, leqs.ω, leqs.L, timespan, timestep, q₀, p₀;
                     v̄=v̄, parameters=parameters)
    else
        LODEEnsemble(toda_lattice_ϑ, toda_lattice_f, toda_lattice_g, ω!, lagrangian,
                     timespan, timestep, q₀, p₀; v̄=v̄, parameters=parameters)
    end
end

function lodeensemble(q₀, p₀; kwargs...)
    @assert length(q₀) == length(p₀)
    lodeensemble(_length(q₀), q₀, p₀; kwargs...)
end

end
