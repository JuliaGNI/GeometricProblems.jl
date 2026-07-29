@doc raw"""
# Harmonic Oscillator

"""
module HarmonicOscillator

using GeometricEquations
using GeometricSolutions
using Parameters

export odeproblem, podeproblem, hodeproblem, iodeproblem, lodeproblem, sodeproblem,
    daeproblem, pdaeproblem, hdaeproblem, idaeproblem, ldaeproblem,
    degenerate_iodeproblem, degenerate_lodeproblem
export deleproblem_midpoint, deleproblem_trapezoidal

export odeensemble, podeensemble, hodeensemble

export hamiltonian, lagrangian

export compute_energy_error, exact_solution

export default_parameters

export plot_solution, plot_spring, plot_phase_portrait, plot_traces, plot_hamiltonian
export labels_ode, labels_hamiltonian


const t₀ = 0.0
const Δt = 0.1
const nt = 10

const DEFAULT_TIMESPAN = (t₀, Δt * nt)
const DEFAULT_TIMESTEP = Δt

const m = 1.0
const k = 0.5
const ω = √(k / m)

default_parameters(::Type{T}=Float64) where {T} = (m=T(m), k=T(k), ω=T(ω))


# Components of the one-form of the *degenerate* (phase-space) formulation, in which the state is
# q = (x, ẋ) and ϑ = (m ẋ, 0). The mass has to appear here for the same reason it appears in
# `degenerate_oscillator_iode_ϑ`: both are ∂L/∂v of `degenerate_lagrangian`.
ϑ₁(t, q, params) = params.m * q[2]
ϑ₂(t, q, params) = zero(eltype(q))

function ϑ(q, params)
    p = zero(q)
    p[1] = ϑ₁(0, q, params)
    p[2] = ϑ₂(0, q, params)
    return p
end

# The Lagrange two-form ω = dθ of the *regular* Lagrangian L = m v²/2 - k q²/2, where θ = ϑᵢ dqⁱ.
#
# A regular Lagrangian is a second-order system of n equations, equivalent to a first-order system
# of 2n, so its ω is the 2n×2n form on (q, q̇) — here 2×2 for the single degree of freedom. In block
# form it is [∂ϑᵢ/∂qⱼ - ∂ϑⱼ/∂qᵢ   -Mᵀ; M   0]; since ϑ = ∂L/∂v = m v depends on the velocity alone,
# the upper-left block vanishes and what remains is the canonical [0 -m; m 0].
#
# This is the same convention EulerLagrange's `LagrangianSystem` produces. Degenerate formulations
# are first-order systems of n equations and carry an n×n ω instead — see `degenerate_ω!` below and
# LotkaVolterra2d/4d, PointVortices and the massless charged particle.
function ω!(Ω, t, q, params)
    @unpack m = params
    Ω[1, 1] = 0
    Ω[1, 2] = -m
    Ω[2, 1] = +m
    Ω[2, 2] = 0
    nothing
end

# LODE/LDAE evaluate the two-form with an extra velocity slot; ω depends only on q.
ω!(Ω, t, q, v, params) = ω!(Ω, t, q, params)

# The two-form of the degenerate formulation, in which the state is q = (x, ẋ) and the system is
# first order, so ω is n×n = 2×2. With ϑ = (m q₂, 0) and the convention Ωᵢⱼ = ∂ϑᵢ/∂qⱼ - ∂ϑⱼ/∂qᵢ
# (shared with LotkaVolterra2d/4d, PointVortices and the massless charged particle) it reads
# Ω = [0 m; -m 0], so that the Euler-Lagrange equations Ω v = -∇H reproduce
# `degenerate_oscillator_iode_v`.
function degenerate_ω!(Ω, t, q, params)
    @unpack m = params
    Ω[1, 1] = 0
    Ω[1, 2] = +m
    Ω[2, 1] = -m
    Ω[2, 2] = 0
    nothing
end

degenerate_ω!(Ω, t, q, v, params) = degenerate_ω!(Ω, t, q, params)

function hamiltonian(t, q::AbstractArray, params)
    @unpack m, k = params
    m * q[2]^2 / 2 + k * q[1]^2 / 2
end

function hamiltonian(t, q::Number, p::Number, params)
    @unpack m, k = params
    p^2 / (2m) + k * q^2 / 2
end

function hamiltonian(t, q::AbstractArray, p::AbstractArray, params)
    hamiltonian(t, q[1], p[1], params)
end

function lagrangian(t, q::Number, v::Number, params)
    @unpack m, k = params
    m * v^2 / 2 - k * q^2 / 2
end

function lagrangian(t, q::AbstractArray, v::AbstractArray, params)
    lagrangian(t, q[1], v[1], params)
end

function degenerate_lagrangian(t, q, v, params)
    ϑ₁(t, q, params) * v[1] + ϑ₂(t, q, params) * v[2] - hamiltonian(t, q, params)
end


A(q, p, params) = q * sqrt(1 + p^2 / q^2 / (params.m * params.k))
ϕ(q, p, params) = atan(p / q / params.ω / params.m)

exact_solution_q(t, q₀, p₀, t₀, params) = A(q₀, p₀, params) * cos(params.ω * (t - t₀) - ϕ(q₀, p₀, params))
exact_solution_p(t, q₀, p₀, t₀, params) = -params.m * params.ω * A(q₀, p₀, params) * sin(params.ω * (t - t₀) - ϕ(q₀, p₀, params))

exact_solution_q(t, q₀::AbstractVector, p₀::AbstractVector, t₀, params) = exact_solution_q(t, q₀[1], p₀[1], t₀, params)
exact_solution_p(t, q₀::AbstractVector, p₀::AbstractVector, t₀, params) = exact_solution_p(t, q₀[1], p₀[1], t₀, params)

exact_solution_q(t, x₀::AbstractVector, t₀, params) = exact_solution_q(t, x₀[1], params.m * x₀[2], t₀, params)
exact_solution_p(t, x₀::AbstractVector, t₀, params) = exact_solution_p(t, x₀[1], params.m * x₀[2], t₀, params) / params.m
exact_solution(t, x₀::AbstractVector, t₀, params) = [exact_solution_q(t, x₀, t₀, params), exact_solution_p(t, x₀, t₀, params)]


const q₀ = [0.5]
const p₀ = [0.0]
const x₀ = vcat(q₀, p₀)

const xmin = [-2.0, -2.0]
const xmax = [+2.0, +2.0]
const nsamples = [10, 10]

const reference_solution_q = exact_solution_q(Δt * nt, q₀[1], p₀[1], t₀, default_parameters())
const reference_solution_p = exact_solution_p(Δt * nt, q₀[1], p₀[1], t₀, default_parameters())

const reference_solution = [reference_solution_q, reference_solution_p]


function _ode_samples(qmin, qmax, nsamples)
    qs = [range(qmin[i], qmax[i]; length=nsamples[i]) for i in eachindex(qmin, qmax, nsamples)]

    samples = vec(collect.(collect(Base.Iterators.product(qs...))))

    (q=samples,)
end

function _pode_samples(qmin, qmax, pmin, pmax, qsamples, psamples)
    qs = [range(qmin[i], qmax[i]; length=qsamples[i]) for i in eachindex(qmin, qmax, qsamples)]
    ps = [range(pmin[i], pmax[i]; length=psamples[i]) for i in eachindex(pmin, pmax, psamples)]

    qsamples = vec(collect.(collect(Base.Iterators.product(qs...))))
    psamples = vec(collect.(collect(Base.Iterators.product(ps...))))
    zsamples = Base.Iterators.product(qsamples, psamples)

    return (
        q=vec([zs[1] for zs in zsamples]),
        p=vec([zs[2] for zs in zsamples]),
    )
end


function oscillator_ode_v(v, t, x, params)
    @unpack m, k = params
    v[1] = x[2]
    v[2] = -k * x[1] / m
    nothing
end

function odeproblem(x₀=x₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(x₀) == 2
    ODEProblem(oscillator_ode_v, timespan, timestep, x₀; invariants=(h=hamiltonian,), parameters=parameters)
end

function odeensemble(qmin=xmin, qmax=xmax, nsamples=nsamples; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    samples = _ode_samples(qmin, qmax, nsamples)
    ODEEnsemble(oscillator_ode_v, timespan, timestep, samples...; invariants=(h=hamiltonian,), parameters=parameters)
end

function exact_solution!(sol::GeometricSolution, prob::ODEProblem)
    for n in eachtimestep(sol)
        sol[n].q .= exact_solution(sol[n].t, sol[0].q, sol[0].t, parameters(prob))
    end
    return sol
end

function exact_solution(prob::ODEProblem)
    exact_solution!(GeometricSolution(prob), prob)
end


function oscillator_pode_v(v, t, q, p, params)
    @unpack m = params
    v[1] = p[1] / m
    nothing
end

function oscillator_pode_f(f, t, q, p, params)
    @unpack k = params
    f[1] = -k * q[1]
    nothing
end

function podeproblem(q₀=q₀, p₀=p₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == 1
    PODEProblem(oscillator_pode_v, oscillator_pode_f, timespan, timestep, q₀, p₀; invariants=(h=hamiltonian,), parameters=parameters)
end

function hodeproblem(q₀=q₀, p₀=p₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == 1
    HODEProblem(oscillator_pode_v, oscillator_pode_f, hamiltonian, timespan, timestep, q₀, p₀; parameters=parameters)
end


function podeensemble(
    qmin=[xmin[1]],
    qmax=[xmax[1]],
    pmin=[xmin[2]],
    pmax=[xmax[2]],
    qsamples=[nsamples[1]],
    psamples=[nsamples[2]],
    ;
    timespan=DEFAULT_TIMESPAN,
    timestep=DEFAULT_TIMESTEP,
    parameters=default_parameters())
    samples = _pode_samples(qmin, qmax, pmin, pmax, qsamples, psamples)
    PODEEnsemble(oscillator_pode_v, oscillator_pode_f, timespan, timestep, samples...; invariants=(h=hamiltonian,), parameters=parameters)
end

function hodeensemble(
    qmin=[xmin[1]],
    qmax=[xmax[1]],
    pmin=[xmin[2]],
    pmax=[xmax[2]],
    qsamples=[nsamples[1]],
    psamples=[nsamples[2]],
    ;
    timespan=DEFAULT_TIMESPAN,
    timestep=DEFAULT_TIMESTEP,
    parameters=default_parameters())
    samples = _pode_samples(qmin, qmax, pmin, pmax, qsamples, psamples)
    HODEEnsemble(oscillator_pode_v, oscillator_pode_f, hamiltonian, timespan, timestep, samples...; parameters=parameters)
end

function exact_solution!(sol::GeometricSolution, prob::Union{PODEProblem,HODEProblem})
    for n in eachtimestep(sol)
        sol[n].q = [exact_solution_q(sol[n].t, sol[0].q, sol[0].p, sol[0].t, parameters(prob))]
        sol[n].p = [exact_solution_p(sol[n].t, sol[0].q, sol[0].p, sol[0].t, parameters(prob))]
    end
    return sol
end

function exact_solution(prob::Union{PODEProblem,HODEProblem})
    exact_solution!(GeometricSolution(prob), prob)
end


function oscillator_sode_v_1(v, t, q, params)
    v[1] = q[2]
    v[2] = 0
    nothing
end

function oscillator_sode_v_2(v, t, q, params)
    @unpack m, k = params
    v[1] = 0
    v[2] = -k * q[1] / m
    nothing
end

function oscillator_sode_q_1(q₁, t₁, q₀, t₀, params)
    q₁[1] = q₀[1] + (t₁ - t₀) * q₀[2]
    q₁[2] = q₀[2]
    nothing
end

function oscillator_sode_q_2(q₁, t₁, q₀, t₀, params)
    @unpack m, k = params
    q₁[1] = q₀[1]
    q₁[2] = q₀[2] - (t₁ - t₀) * k * q₀[1] / m
    nothing
end

function sodeproblem(x₀=x₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    SODEProblem((oscillator_sode_v_1, oscillator_sode_v_2),
        (oscillator_sode_q_1, oscillator_sode_q_2),
        timespan, timestep, x₀; v̄=oscillator_ode_v, parameters=parameters)
end


function oscillator_iode_ϑ(p, t, q, v, params)
    @unpack m = params
    p[1] = m * v[1]
    nothing
end

function oscillator_iode_f(f, t, q, v, params)
    @unpack k = params
    f[1] = -k * q[1]
    nothing
end

function oscillator_iode_g(g, t, q, v, λ, params)
    g[1] = λ[1]
    nothing
end

function oscillator_iode_v(v, t, q, p, params)
    @unpack m = params
    v[1] = p[1] / m
    nothing
end

function iodeproblem(q₀=q₀, p₀=p₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == 1
    IODEProblem(oscillator_iode_ϑ, oscillator_iode_f,
        oscillator_iode_g, timespan, timestep, q₀, p₀;
        invariants=(h=hamiltonian,), parameters=parameters,
        v̄=oscillator_iode_v)
end

function lodeproblem(q₀=q₀, p₀=p₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == 1
    LODEProblem(oscillator_iode_ϑ, oscillator_iode_f,
        oscillator_iode_g, ω!, lagrangian,
        timespan, timestep, q₀, p₀;
        invariants=(h=hamiltonian,), parameters=parameters,
        v̄=oscillator_iode_v)
end



function degenerate_oscillator_iode_ϑ(p, t, q, v, params)
    @unpack m = params
    p[1] = m * q[2]
    p[2] = 0
    nothing
end

function degenerate_oscillator_iode_f(f, t, q, v, params)
    @unpack m, k = params
    f[1] = -k * q[1]
    f[2] = m * (v[1] - q[2])
    nothing
end

function degenerate_oscillator_iode_g(g, t, q, v, λ, params)
    g[1] = 0
    g[2] = λ[1]
    nothing
end

function degenerate_oscillator_iode_v(v, t, q, p, params)
    @unpack m, k = params
    v[1] = q[2]
    v[2] = -k * q[1] / m
    nothing
end

# `p₀ = nothing` resolves to ϑ(q₀) evaluated with the *given* parameters. A positional default
# cannot reference the `parameters` keyword, and hard-coding `ϑ(q₀)` with the default parameters
# would silently produce a momentum off by a factor `m` whenever `m ≠ 1`.
function degenerate_iodeproblem(q₀=x₀, p₀=nothing; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    p₀ = p₀ === nothing ? ϑ(q₀, parameters) : p₀
    @assert length(q₀) == length(p₀) == 2
    IODEProblem(degenerate_oscillator_iode_ϑ, degenerate_oscillator_iode_f,
        degenerate_oscillator_iode_g, timespan, timestep, q₀, p₀;
        invariants=(h=hamiltonian,), parameters=parameters,
        v̄=degenerate_oscillator_iode_v)
end

function degenerate_lodeproblem(q₀=x₀, p₀=nothing; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    p₀ = p₀ === nothing ? ϑ(q₀, parameters) : p₀
    @assert length(q₀) == length(p₀) == 2
    LODEProblem(degenerate_oscillator_iode_ϑ, degenerate_oscillator_iode_f,
        degenerate_oscillator_iode_g, degenerate_ω!, degenerate_lagrangian,
        timespan, timestep, q₀, p₀;
        invariants=(h=hamiltonian,), parameters=parameters,
        v̄=degenerate_oscillator_iode_v)
end



function oscillator_dae_u(u, t, x, λ, params)
    @unpack m, k = params
    u[1] = k * x[1] * λ[1]
    u[2] = m * x[2] * λ[1]
end

@doc raw"""
Build the energy constraint ``\phi = H(t, \cdot) - H_0`` of the differential-algebraic
formulations, where ``H_0`` is the energy of the *problem's own* initial condition.

`H₀` is captured in a closure per problem rather than read from the module-level defaults, so
that the constraint vanishes at ``t_0`` for any initial data. A shared constraint closing over
the default initial condition would be satisfied only for that one initial condition.

`_dae_energy_constraint` covers the plain DAE, whose state `x = (q, v)` carries the velocity, so
that `H` is evaluated in its velocity form.
"""
_dae_energy_constraint(H₀) = function (ϕ, t, x, params)
    ϕ[1] = hamiltonian(t, x, params) - H₀
    nothing
end

@doc raw"""
Energy constraint for the partitioned, Hamiltonian, implicit and variational DAEs, which keep
position and momentum separate. The implicit/variational forms evaluate it with an additional
velocity slot, which the energy does not depend on.

See `_dae_energy_constraint` for why the reference energy is captured per problem.
"""
function _pdae_energy_constraint(H₀)
    function constraint(ϕ, t, q, p, params)
        ϕ[1] = hamiltonian(t, q, p, params) - H₀
        nothing
    end
    constraint(ϕ, t, q, v, p, params) = constraint(ϕ, t, q, p, params)
    return constraint
end

function daeproblem(x₀=x₀, λ₀=[zero(eltype(x₀))]; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    constraint = _dae_energy_constraint(hamiltonian(timespan[begin], x₀, parameters))
    DAEProblem(oscillator_ode_v, oscillator_dae_u, constraint, timespan, timestep, x₀, λ₀;
        v̄=oscillator_ode_v, invariants=(h=hamiltonian,), parameters=parameters)
end


function oscillator_pdae_v(v, t, q, p, params)
    @unpack m = params
    v[1] = p[1] / m
    nothing
end

function oscillator_pdae_f(f, t, q, p, params)
    @unpack k = params
    f[1] = -k * q[1]
    nothing
end

function oscillator_pdae_u(u, t, q, p, λ, params)
    @unpack k = params
    u[1] = k * q[1] * λ[1]
    nothing
end

function oscillator_pdae_g(g, t, q, p, λ, params)
    @unpack m = params
    g[1] = p[1] * λ[1] / m
    nothing
end

function oscillator_pdae_ū(u, t, q, p, λ, params)
    @unpack k = params
    u[1] = k * q[1] * λ[1]
    nothing
end

function oscillator_pdae_ḡ(g, t, q, p, λ, params)
    g[1] = p[1] * λ[1] / params.m
    nothing
end

# The energy constraint is built per problem, see `_pdae_energy_constraint`.

function oscillator_pdae_ψ(ψ, t, q, p, q̇, ṗ, params)
    @unpack k = params
    ψ[1] = p[1] * ṗ[1] / params.m + k * q[1] * q̇[1]
    nothing
end

function pdaeproblem(q₀=q₀, p₀=p₀, λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == 1
    constraint = _pdae_energy_constraint(hamiltonian(timespan[begin], q₀, p₀, parameters))
    PDAEProblem(oscillator_pdae_v, oscillator_pdae_f,
        oscillator_pdae_u, oscillator_pdae_g, constraint,
        timespan, timestep, q₀, p₀, λ₀; invariants=(h=hamiltonian,), parameters=parameters)
end

function hdaeproblem(q₀=q₀, p₀=p₀, λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == 1
    constraint = _pdae_energy_constraint(hamiltonian(timespan[begin], q₀, p₀, parameters))
    HDAEProblem(oscillator_pdae_v, oscillator_pdae_f,
        oscillator_pdae_u, oscillator_pdae_g, constraint,
        oscillator_pdae_ū, oscillator_pdae_ḡ, oscillator_pdae_ψ,
        hamiltonian, timespan, timestep, q₀, p₀, λ₀; parameters=parameters)
end


oscillator_idae_u(u, t, q, v, p, λ, params) = oscillator_pdae_u(u, t, q, p, λ, params)
oscillator_idae_g(g, t, q, v, p, λ, params) = oscillator_pdae_g(g, t, q, p, λ, params)
oscillator_idae_ū(u, t, q, v, p, λ, params) = oscillator_pdae_ū(u, t, q, p, λ, params)
oscillator_idae_ḡ(g, t, q, v, p, λ, params) = oscillator_pdae_ḡ(g, t, q, p, λ, params)
# The energy constraint of the implicit/variational forms is the same closure as for the
# partitioned ones; it already accepts the extra velocity slot.
oscillator_idae_ψ(ψ, t, q, v, p, q̇, ṗ, params) = oscillator_pdae_ψ(ψ, t, q, p, q̇, ṗ, params)

function idaeproblem(q₀=q₀, p₀=p₀, λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == length(λ₀) == 1
    constraint = _pdae_energy_constraint(hamiltonian(timespan[begin], q₀, p₀, parameters))
    IDAEProblem(oscillator_iode_ϑ, oscillator_iode_f,
        oscillator_idae_u, oscillator_idae_g, constraint,
        timespan, timestep, q₀, p₀, λ₀; v̄=oscillator_iode_v, invariants=(h=hamiltonian,), parameters=parameters)
end

function ldaeproblem(q₀=q₀, p₀=p₀, λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == length(p₀) == length(λ₀) == 1
    constraint = _pdae_energy_constraint(hamiltonian(timespan[begin], q₀, p₀, parameters))
    LDAEProblem(oscillator_iode_ϑ, oscillator_iode_f,
        oscillator_idae_u, oscillator_idae_g, constraint, ω!, lagrangian,
        timespan, timestep, q₀, p₀, λ₀; v̄=oscillator_iode_v, invariants=(h=hamiltonian,), parameters=parameters)
end


function oscillator_dele_midpoint_Ld(t₀, t₁, q₀, q₁, params)
    h = (t₁ - t₀)
    t = (t₀ + t₁) / 2
    q = (q₀[1] + q₁[1]) / 2
    v = (q₁[1] - q₀[1]) / h
    return h * lagrangian(t, q, v, params)
end

function oscillator_dele_midpoint_D1Ld(d, t₀, t₁, q₀, q₁, params)
    @unpack m, k = params
    h = (t₁ - t₀)
    q = (q₀[1] + q₁[1]) / 2
    v = (q₁[1] - q₀[1]) / h
    d[1] = -m * v - h * k * q / 2
    return nothing
end

function oscillator_dele_midpoint_D2Ld(d, t₀, t₁, q₀, q₁, params)
    @unpack m, k = params
    h = (t₁ - t₀)
    q = (q₀[1] + q₁[1]) / 2
    v = (q₁[1] - q₀[1]) / h
    d[1] = +m * v - h * k * q / 2
    return nothing
end

function deleproblem_midpoint(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == 1

    q₁ = [exact_solution_q(timespan[begin] - Δt, q₀, zero(q₀), timespan[begin], parameters)]

    DELEProblem(oscillator_dele_midpoint_Ld,
        oscillator_dele_midpoint_D1Ld,
        oscillator_dele_midpoint_D2Ld,
        timespan, timestep, q₁, q₀; invariants=(h=hamiltonian,), parameters=parameters)
end


function oscillator_dele_trapezoidal_Ld(t₀, t₁, q₀, q₁, params)
    h = (t₁ - t₀)
    v = (q₁[1] - q₀[1]) / h
    return h * (lagrangian(t₀, q₀[1], v, params) + lagrangian(t₁, q₁[1], v, params)) / 2
end

function oscillator_dele_trapezoidal_D1Ld(d, t₀, t₁, q₀, q₁, params)
    @unpack m, k = params
    h = (t₁ - t₀)
    v = (q₁[1] - q₀[1]) / h
    d[1] = -m * v - h * k * q₀[1] / 2
    return nothing
end

function oscillator_dele_trapezoidal_D2Ld(d, t₀, t₁, q₀, q₁, params)
    @unpack m, k = params
    h = (t₁ - t₀)
    v = (q₁[1] - q₀[1]) / h
    d[1] = +m * v - h * k * q₁[1] / 2
    return nothing
end

function deleproblem_trapezoidal(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    @assert length(q₀) == 1

    q₁ = [exact_solution_q(timespan[begin] - Δt, q₀, zero(q₀), timespan[begin], parameters)]

    DELEProblem(oscillator_dele_trapezoidal_Ld,
        oscillator_dele_trapezoidal_D1Ld,
        oscillator_dele_trapezoidal_D2Ld,
        timespan, timestep, q₁, q₀; invariants=(h=hamiltonian,), parameters=parameters)
end


function exact_solution(probs::Union{ODEEnsemble,PODEEnsemble,HODEEnsemble})
    sols = EnsembleSolution(probs)
    for (sol, prob) in zip(sols, probs)
        exact_solution!(sol, prob)
    end
    return sols
end


function compute_energy_error(t, q::DataSeries{T}, params) where {T}
    h = DataSeries(T, q.nt)
    e = DataSeries(T, q.nt)

    for i in axes(q, 2)
        h[i] = hamiltonian(t[i], q[:, i], params)
        e[i] = (h[i] - h[0]) / h[0]
    end

    (h, e)
end


const labels_ode = (t = "t", q = "x", p = "ẋ", h = "E")
const labels_hamiltonian = (t = "t", q = "q", p = "p", h = "H")

function plot_spring end
function plot_solution end
function plot_phase_portrait end
function plot_traces end
function plot_hamiltonian end

end
