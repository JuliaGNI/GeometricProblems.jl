@doc raw"""
# Kubo Oscillator

The Kubo oscillator is a unit-frequency harmonic oscillator, $\dot{q}_1 = q_2$, $\dot{q}_2 = -q_1$,
driven by multiplicative (Stratonovich) noise. It is a standard test problem for stochastic
geometric integrators. The module provides it as a stochastic differential equation
(`sdeproblem`), a partitioned SDE (`psdeproblem`), a split partitioned SDE (`spsdeproblem`), and
the underlying deterministic ODE (`odeproblem`). The `*problem` builders create a single problem
from one initial condition; the `*ensemble` builders (`sdeensemble`, `psdeensemble`,
`spsdeensemble`) build an ensemble from several initial conditions. The noise is a
one-dimensional `WienerProcess`.

Because the diffusion is proportional to the drift, the solution is the deterministic one
evaluated at the random time $\theta(t) = t + \nu W(t)$, so the energy is conserved exactly along
every sample path. That is what makes this problem useful for checking a stochastic geometric
integrator: any energy drift is the scheme's, not the problem's. `exact_solution` gives the
closed form.

The `damped_*` builders add linear damping, following Kraus & Tyranowski, *Variational
integrators for stochastic dissipative Hamiltonian systems*. There the system is
$H = (p^2+q^2)/2$, $h = \nu (p^2+q^2)/2$ with forcing $F(q,p) = -\gamma p$ and
$f(q,p) = -\nu \gamma p$. The energy then decays rather than being conserved, and
`exact_mean_energy` gives $E(H)$ in closed form. Unlike the undamped problem, the damped one has
genuinely non-zero $f_2$ and $G_2$, so it exercises the split half of the SPSDE.

System parameters: `ν` — the noise intensity; `γ` — the damping coefficient (`0` for the
undamped problems). Note the paper writes these as $\beta$ and $\nu$ respectively; the names
here follow this module's existing convention, in which `ν` has always been the noise intensity.
"""
module KuboOscillator

using GeometricEquations
using Parameters

export sdeproblem, psdeproblem, spsdeproblem
export sdeensemble, psdeensemble, spsdeensemble
export damped_psdeproblem, damped_spsdeproblem
export odeproblem
export exact_solution, exact_solution_q, exact_solution_p, exact_mean_energy

const q_init_A=[0.5, 0.0]
const q_init_B=[[0.5, 0.0],
    [0.0, 0.5],
    [-0.5, 0.0]]

const noise_intensity = 0.1

const Δt = 0.01
const nt = 10
const DEFAULT_TIMESPAN = (0.0, Δt*nt)
const DEFAULT_TIMESTEP = Δt

default_parameters(::Type{T} = Float64) where {T} = (ν = T(noise_intensity),)

# Parameters of the damped oscillator, from Kraus & Tyranowski §4.1. Their β is this module's ν.
const damped_noise_intensity = 0.5
const damping_coefficient = 0.001

# The damping is weak — with γ = 0.001 the energy decays on a timescale 1/γ = 1000 — so over the
# undamped problems' timespan of 0.1 it moves E(H) by 5e-06 relative and the damped problem is
# numerically indistinguishable from the undamped one. The damped builders therefore carry their
# own default, long enough for the decay to be visible (about 10% by the final time).
const damped_Δt = 0.1
const damped_nt = 1000
const DEFAULT_DAMPED_TIMESPAN = (0.0, damped_Δt*damped_nt)
const DEFAULT_DAMPED_TIMESTEP = damped_Δt

const q_init_damped = [2.0]
const p_init_damped = [0.0]

function damped_parameters(::Type{T} = Float64;
        ν = damped_noise_intensity, γ = damping_coefficient) where {T}
    (ν = T(ν), γ = T(γ))
end

# The closed forms below are the *underdamped* solution: ω is real only for |γ| < 2, and at
# |γ| = 2 the sin(ωθ)/ω terms are 0/0 and would return NaN rather than fail. Rejecting the whole
# range keeps the one silent case out.
function damped_frequency(γ)
    abs(γ) < 2 || throw(DomainError(γ,
        "the Kubo oscillator's closed-form solution is the underdamped one and needs |γ| < 2"))
    sqrt(4 - γ^2) / 2
end

@doc raw"""
    exact_solution_q(t, W, q₀, p₀, t₀, params)
    exact_solution_p(t, W, q₀, p₀, t₀, params)
    exact_solution(t, W, q₀, p₀, t₀, params)
    exact_solution(t, W, x₀, t₀, params)

Exact solution of the (possibly damped) Kubo oscillator at time `t`, along the sample path whose
Wiener increment $W(t) - W(t_0)$ is `W`.

Because the diffusion is proportional to the drift, the solution is the deterministic damped
oscillator evaluated at the random time $\theta = (t - t_0) + \nu (W(t) - W(t_0))$:
```math
\begin{aligned}
q(t) &= e^{-\gamma \theta / 2} \left[ q_0 \cos \omega \theta
        + \tfrac{1}{\omega} \left( p_0 + \tfrac{\gamma}{2} q_0 \right) \sin \omega \theta \right], \\
p(t) &= e^{-\gamma \theta / 2} \left[ p_0 \cos \omega \theta
        - \tfrac{1}{\omega} \left( q_0 + \tfrac{\gamma}{2} p_0 \right) \sin \omega \theta \right],
\end{aligned}
```
with $\omega = \tfrac{1}{2} \sqrt{4 - \gamma^2}$. `params` needs a noise intensity `ν`; a
damping coefficient `γ` is taken as zero when absent, which reduces the above to a rotation by
$\theta$. This is the underdamped solution, so $|\gamma| < 2$; anything else is a `DomainError`.

`q₀` and `p₀` are numbers or one-element vectors, and the result matches: a tuple of numbers or a
tuple of vectors. The form taking a single state vector `x₀` serves the `sdeproblem` formulation,
where $q = (q_1, q_2)$ plays the role of $(q, p)$, and returns `[q, p]`.
"""
function exact_solution(t, W, q₀::Number, p₀::Number, t₀, params)
    @unpack ν = params
    γ = hasproperty(params, :γ) ? params.γ : zero(ν)
    θ = (t - t₀) + ν * W
    ω = damped_frequency(γ)
    decay = exp(-γ * θ / 2)

    q = decay * (q₀ * cos(ω * θ) + (p₀ + γ * q₀ / 2) * sin(ω * θ) / ω)
    p = decay * (p₀ * cos(ω * θ) - (q₀ + γ * p₀ / 2) * sin(ω * θ) / ω)

    (q, p)
end

"""
    exact_solution_q(t, W, q₀, p₀, t₀, params)

Position component of [`exact_solution`](@ref).
"""
exact_solution_q(t, W, q₀, p₀, t₀, params) = exact_solution(t, W, q₀, p₀, t₀, params)[1]

"""
    exact_solution_p(t, W, q₀, p₀, t₀, params)

Momentum component of [`exact_solution`](@ref).
"""
exact_solution_p(t, W, q₀, p₀, t₀, params) = exact_solution(t, W, q₀, p₀, t₀, params)[2]

function exact_solution(t, W, q₀::AbstractVector, p₀::AbstractVector, t₀, params)
    q, p = exact_solution(t, W, q₀[begin], p₀[begin], t₀, params)
    ([q], [p])
end

function exact_solution(t, W, x₀::AbstractVector, t₀, params)
    q, p = exact_solution(t, W, x₀[begin], x₀[begin + 1], t₀, params)
    [q, p]
end

@doc raw"""
    exact_mean_energy(t, q₀, p₀, t₀, params)

Expected value of the Hamiltonian $H = (p^2 + q^2)/2$ of the damped Kubo oscillator at time `t`,
in closed form (Kraus & Tyranowski §4.1), written in the elapsed time $s = t - t_0$:
```math
E(H) = a \, e^{-\frac{\gamma (2 - \nu^2 \gamma)}{2} s}
     + e^{-((2 - \gamma^2)\nu^2 + \gamma) s}
       \Big[ b \cos \big( 2 (1 - \nu^2 \gamma) \omega s \big)
           + c \sin \big( 2 (1 - \nu^2 \gamma) \omega s \big) \Big],
```
with $\omega = \tfrac{1}{2}\sqrt{4 - \gamma^2}$ and

```math
a = \frac{2 (p_0^2 + q_0^2 + \gamma p_0 q_0)}{4 - \gamma^2}, \qquad
b = -\frac{\gamma^2 (p_0^2 + q_0^2) + 4 \gamma p_0 q_0}{2 (4 - \gamma^2)}, \qquad
c = \frac{\gamma (q_0^2 - p_0^2)}{2 \sqrt{4 - \gamma^2}} .
```

This averages over the noise, so unlike [`exact_solution`](@ref) it needs no sample path. As
there, $|\gamma| < 2$. Without damping it collapses to the constant $(p_0^2 + q_0^2)/2$, which is
the exact conservation the undamped problems exhibit pathwise.
"""
function exact_mean_energy(t, q₀::Number, p₀::Number, t₀, params)
    @unpack ν = params
    γ = hasproperty(params, :γ) ? params.γ : zero(ν)
    ω = damped_frequency(γ)
    s = t - t₀

    a = 2 * (p₀^2 + q₀^2 + γ * p₀ * q₀) / (4 - γ^2)
    b = -(γ^2 * (p₀^2 + q₀^2) + 4 * γ * p₀ * q₀) / (2 * (4 - γ^2))
    c = γ * (q₀^2 - p₀^2) / (2 * sqrt(4 - γ^2))

    a * exp(-γ * (2 - ν^2 * γ) * s / 2) +
    exp(-((2 - γ^2) * ν^2 + γ) * s) *
    (b * cos(2 * (1 - ν^2 * γ) * ω * s) + c * sin(2 * (1 - ν^2 * γ) * ω * s))
end

function exact_mean_energy(t, q₀::AbstractVector, p₀::AbstractVector, t₀, params)
    exact_mean_energy(t, q₀[begin], p₀[begin], t₀, params)
end

function kubo_oscillator_sde_v(v, t, q, params)
    v[1] = q[2]
    v[2] = -q[1]
end

function kubo_oscillator_sde_B(B::AbstractVector, t, q, params)
    @unpack ν = params
    B[1] = +ν * q[2]
    B[2] = -ν * q[1]
end

function kubo_oscillator_sde_B(B::AbstractMatrix, t, q, params)
    @unpack ν = params
    for j in axes(B, 2)
        B[1, j] = +ν * q[2]
        B[2, j] = -ν * q[1]
    end
end

function sdeproblem(q₀ = q_init_A; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    # single initial condition (q_init_A); the number of sample paths is chosen by the integrator
    # 1-dimensional noise
    SDEProblem(kubo_oscillator_sde_v, kubo_oscillator_sde_B, WienerProcess(1),
        timespan, timestep, q₀; parameters = parameters)
end

# NOTE: GeometricEquations exports `SDEEnsemble`/`PSDEEnsemble`/`SPSDEEnsemble`, but only as
# type aliases — unlike `ODEEnsemble` and friends they have no convenience constructor. The
# ensembles below therefore have to assemble the `EnsembleProblem` by hand, which also means
# reaching for the non-exported `GeometricEquations.parameter_types`. Once upstream adds the
# constructors these three functions collapse into one-liners.

function sdeensemble(q₀ = q_init_B; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    # q_init_B holds several initial conditions -> ensemble problem
    equ = SDE(kubo_oscillator_sde_v, kubo_oscillator_sde_B, WienerProcess(1);
        parameters = GeometricEquations.parameter_types(parameters))
    ics = [(q = StateVariable(x),) for x in q₀]
    GeometricEquations.EnsembleProblem(equ, timespan, timestep, ics, parameters)
end

# ODE

function odeproblem(q₀ = q_init_A; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    ODEProblem(kubo_oscillator_sde_v, timespan, timestep, q₀; parameters = parameters)
end

# PSDE

const q_init_C=[0.5]
const p_init_C=[0.0]

const q_init_D=[[0.5], [0.0], [-0.5]]
const p_init_D=[[0.0], [0.5], [0.0]]

function kubo_oscillator_psde_v(v, t, q, p, params)
    v[1] = p[1]
end

function kubo_oscillator_psde_f(f, t, q, p, params)
    f[1] = -q[1]
end

function kubo_oscillator_psde_B(B, t, q, p, params)
    @unpack ν = params
    B[1, 1] = +ν * p[1]
end

function kubo_oscillator_psde_G(G, t, q, p, params)
    @unpack ν = params
    G[1, 1] = -ν * q[1]
end

function psdeproblem(q₀ = q_init_C, p₀ = p_init_C; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    # single initial condition (q_init_C, p_init_C); the number of sample paths is chosen by the integrator
    # 1-dimensional noise
    PSDEProblem(kubo_oscillator_psde_v, kubo_oscillator_psde_f,
        kubo_oscillator_psde_B, kubo_oscillator_psde_G, WienerProcess(1),
        timespan, timestep, q₀, p₀; parameters = parameters)
end

function psdeensemble(q₀ = q_init_D, p₀ = p_init_D; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    # q_init_D / p_init_D hold several initial conditions -> ensemble problem
    equ = PSDE(kubo_oscillator_psde_v, kubo_oscillator_psde_f,
        kubo_oscillator_psde_B, kubo_oscillator_psde_G, WienerProcess(1);
        parameters = GeometricEquations.parameter_types(parameters))
    ics = [(q = StateVariable(x), p = StateVariable(y)) for (x, y) in zip(q₀, p₀)]
    GeometricEquations.EnsembleProblem(equ, timespan, timestep, ics, parameters)
end

# SPSDE

function kubo_oscillator_spsde_v(v, t, q, p, params)
    v[1] = p[1]
end

function kubo_oscillator_spsde_f1(f, t, q, p, params)
    f[1] = -q[1]
end

function kubo_oscillator_spsde_f2(f, t, q, p, params)
    f[1] = 0
end

function kubo_oscillator_spsde_B(B, t, q, p, params)
    @unpack ν = params
    B[1, 1] = +ν * p[1]
end

function kubo_oscillator_spsde_G1(G, t, q, p, params)
    @unpack ν = params
    G[1, 1] = -ν * q[1]
end

function kubo_oscillator_spsde_G2(G, t, q, p, params)
    G[1, 1] = 0
end

function spsdeproblem(q₀ = q_init_C, p₀ = p_init_C; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    # single initial condition (q_init_C, p_init_C); the number of sample paths is chosen by the integrator
    # 1-dimensional noise
    SPSDEProblem(
        kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2,
        kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, WienerProcess(1),
        timespan, timestep, q₀, p₀; parameters = parameters)
end

function spsdeensemble(q₀ = q_init_D, p₀ = p_init_D; timespan = DEFAULT_TIMESPAN,
        timestep = DEFAULT_TIMESTEP, parameters = default_parameters())
    # q_init_D / p_init_D hold several initial conditions -> ensemble problem
    equ = SPSDE(
        kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2,
        kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, WienerProcess(1);
        parameters = GeometricEquations.parameter_types(parameters))
    ics = [(q = StateVariable(x), p = StateVariable(y)) for (x, y) in zip(q₀, p₀)]
    GeometricEquations.EnsembleProblem(equ, timespan, timestep, ics, parameters)
end

# Damped Kubo oscillator
#
# The forcing terms F(q,p) = -γ p and f(q,p) = -ν γ p of Kraus & Tyranowski §4.1. In the split
# formulation these are exactly f2 and G2, which are identically zero in the undamped problems
# above — so this is the only variant here that exercises the split half of an SPSDE.

function kubo_oscillator_damped_f(f, t, q, p, params)
    @unpack γ = params
    f[1] = -q[1] - γ * p[1]
end

function kubo_oscillator_damped_G(G, t, q, p, params)
    @unpack ν, γ = params
    G[1, 1] = -ν * q[1] - ν * γ * p[1]
end

function kubo_oscillator_damped_f2(f, t, q, p, params)
    @unpack γ = params
    f[1] = -γ * p[1]
end

function kubo_oscillator_damped_G2(G, t, q, p, params)
    @unpack ν, γ = params
    G[1, 1] = -ν * γ * p[1]
end

"""
Damped Kubo oscillator as a partitioned SDE, with the damping folded into `f` and `G`.

Pairs with [`damped_spsdeproblem`](@ref), which splits the same dynamics; the two must give the
same trajectory on a common sample path.
"""
function damped_psdeproblem(q₀ = q_init_damped, p₀ = p_init_damped;
        timespan = DEFAULT_DAMPED_TIMESPAN, timestep = DEFAULT_DAMPED_TIMESTEP,
        parameters = damped_parameters())
    PSDEProblem(kubo_oscillator_psde_v, kubo_oscillator_damped_f,
        kubo_oscillator_psde_B, kubo_oscillator_damped_G, WienerProcess(1),
        timespan, timestep, q₀, p₀; parameters = parameters)
end

"""
Damped Kubo oscillator as a split partitioned SDE, with the Hamiltonian part in `f1`/`G1` and the
damping in `f2`/`G2`.
"""
function damped_spsdeproblem(q₀ = q_init_damped, p₀ = p_init_damped;
        timespan = DEFAULT_DAMPED_TIMESPAN, timestep = DEFAULT_DAMPED_TIMESTEP,
        parameters = damped_parameters())
    SPSDEProblem(
        kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_damped_f2,
        kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_damped_G2,
        WienerProcess(1), timespan, timestep, q₀, p₀; parameters = parameters)
end

end
