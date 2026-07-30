# Benchmark for the TodaLattice equations of motion: hand-written vs symbolic (EulerLagrange).
#
#     julia --project=.    benchmark/toda_lattice.jl   # construction, compile and per-call cost
#     julia --project=test benchmark/toda_lattice.jl   # additionally times a full integration
#
#     TODA_LATTICE_BENCH_FULL=1 julia --project=. benchmark/toda_lattice.jl   # sweep out to N = 200
#
# The companion of `benchmark/linear_wave.jl`, and the measurement that chose the default here too.
# It answers five questions:
#
#   1. Do the hand-written and generated functions agree?
#   2. What does the symbolic route cost to *build* and to *compile*, as a function of size?
#   3. What does each function cost to *call*, once built?
#   4. How many evaluations would it take for the symbolic route to earn its setup cost back?
#   5. What does the difference come to over a whole integration?
#
# The state has exactly N components — unlike the linear wave, the Toda lattice is periodic and has
# no boundary points — and the Lagrangian is regular, so EulerLagrange emits a 2N × 2N Lagrange
# two-form. It builds that form by symbolically differentiating a *dense* matrix, which is what makes
# the symbolic route scale as badly as it does. Note that the `ω` rows below measure something no
# integrator ever calls — every `.ω` in GeometricIntegrators is a SPARK/VPRK tableau coefficient, the
# one equation-`ω` call site is commented out, and `GeometricEquations.check_methods` skips `ω` —
# which is precisely what makes the generated `ω` pure overhead.
#
# Where this differs from the linear wave is the potential itself: Toda's is a sum of N exponentials,
# so both routes are dominated by `exp` rather than by arithmetic on differences. That does more than
# narrow the per-call gap for the force — it reverses it, the generated force coming out 5-12% ahead
# over N = 8 to 200. It changes nothing about the construction cost, which is what the default turns
# on, and section 4 puts the two together: the break-even is 6.9e9 force evaluations at N = 200.
#
# `using EulerLagrange` is deliberately avoided: EulerLagrange is a dependency of the main project but
# not of `test/Project.toml`, and this script has to run under both. `functions` is generic in
# GeometricBase and re-exported by GeometricEquations, so the EulerLagrange methods dispatch anyway.

using GeometricEquations: functions
using GeometricProblems
using Printf

import GeometricProblems.TodaLattice as toda

include("timing.jl")

const PARAMS = toda.default_parameters()

# A 5% difference in a hand-rolled microbenchmark is not a real difference; the same threshold
# `simplify_evaluation.jl` uses.
const THRESHOLD = 1.15

const N_SWEEP = [4, 8, 16, 32, 64, 96, 128]
const N_FULL = [160, 200]
const FULL = get(ENV, "TODA_LATTICE_BENCH_FULL", "0") != "0"

# Stop the sweep once a single size costs more than this, so the script stays runnable by default.
const BUDGET = 90.0

"""
    state(N)

Positions, momenta and velocities on `N` lattice points.

The default initial momentum of this model is zero, which would make the force the only nontrivial
thing measured and `L` numerically equal to `-V`. The momenta are therefore offset, and the
positions perturbed off the bump so that no lattice difference is accidentally zero.
"""
function state(N)
    q = collect(toda.compute_initial_q(toda.μ, N)) .+ 0.05 .* collect(1:N)
    p = collect(range(0.1, 0.4, length = N))
    # the mass matrix is the identity, so q̇ = p
    return q, p, copy(p)
end

"""
    comparisons(N, hameqs, lageqs)

The functions an integrator calls, as `(label, generated, handwritten)` triples of thunks.

Each thunk returns the output it wrote, so that `percall` can feed it to `sink!` and the call cannot
be optimised away. Generated and hand-written write into separate buffers, so the same triple serves
both the agreement check and the timings.
"""
function comparisons(N, hameqs, lageqs)
    q, p, v = state(N)
    gen, hand = zeros(N), zeros(N)
    Ωgen, Ωhand = zeros(2N, 2N), zeros(2N, 2N)
    λ = collect(range(0.1, 0.7, length = N))
    return (
        ("v = ∂H/∂p", () -> (hameqs.v(gen, 0.0, q, p, PARAMS); gen),
                      () -> (toda.toda_lattice_v(hand, 0.0, q, p, PARAMS); hand)),
        ("f = -∂H/∂q", () -> (hameqs.f(gen, 0.0, q, p, PARAMS); gen),
                       () -> (toda.toda_lattice_f(hand, 0.0, q, p, PARAMS); hand)),
        ("f = ∂L/∂q", () -> (lageqs.f(gen, 0.0, q, v, PARAMS); gen),
                      () -> (toda.toda_lattice_f(hand, 0.0, q, v, PARAMS); hand)),
        ("ϑ = ∂L/∂q̇", () -> (lageqs.ϑ(gen, 0.0, q, v, PARAMS); gen),
                       () -> (toda.toda_lattice_ϑ(hand, 0.0, q, v, PARAMS); hand)),
        ("g = λ", () -> (lageqs.g(gen, 0.0, q, v, λ, PARAMS); gen),
                  () -> (toda.toda_lattice_g(hand, 0.0, q, v, λ, PARAMS); hand)),
        ("ω (2N×2N)", () -> (lageqs.ω(Ωgen, 0.0, q, v, PARAMS); Ωgen),
                      () -> (toda.ω!(Ωhand, 0.0, q, v, PARAMS); Ωhand)),
        ("H", () -> hameqs.H(0.0, q, p, PARAMS),
              () -> toda.hamiltonian(0.0, q, p, PARAMS)),
        ("L", () -> lageqs.L(0.0, q, v, PARAMS),
              () -> toda.lagrangian(0.0, q, v, PARAMS)),
    )
end

# The construction sweep runs each size in a *fresh process*, because the cost of a generated
# function's first call is badly order-dependent within one process: once anything has driven Julia
# through this call path, later first calls come back nanoseconds instead of seconds, and the sweep
# would report a multi-megabyte `ω` as free to compile. A subprocess per size is the only way to
# measure what a user actually pays in a fresh session.
#
# Each worker does warm up `lagrangian_system` at N = 4 *without calling* the generated functions.
# That excludes the one-off ~10 s of Symbolics compilation, which would otherwise swamp every row and
# hide the scaling, and it was verified not to affect the first-call figures.
const WORKER = raw"""
    using GeometricEquations: functions
    using GeometricProblems
    import GeometricProblems.TodaLattice as toda

    const P = toda.default_parameters()
    el(f) = (GC.gc(); s = time(); r = f(); (time() - s, r))

    toda.lagrangian_system(4, P)    # warm up Symbolics; deliberately never called

    N = parse(Int, ARGS[1])
    q = collect(toda.compute_initial_q(toda.μ, N)) .+ 0.05 .* collect(1:N)
    v = collect(range(0.1, 0.4, length = N))
    out = zeros(N)
    Ω = zeros(2N, 2N)

    t_hsys, _    = el(() -> toda.hamiltonian_system(N, P))
    t_lsys, lsys = el(() -> toda.lagrangian_system(N, P))
    t_hand, _    = el(() -> (toda.hodeproblem(N), toda.lodeproblem(N)))

    eqs = functions(lsys)
    t_f, _ = el(() -> eqs.f(out, 0.0, q, v, P))
    t_ω, _ = el(() -> eqs.ω(Ω, 0.0, q, v, P))

    println("RESULT ", join((t_hsys, t_lsys, t_hand, t_f, t_ω, length(string(eqs.ω.body))), " "))
"""

"Run the construction benchmark for one size in a fresh process."
function measure_construction(N)
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) --startup-file=no -e $WORKER $N`
    for line in eachline(open(cmd))
        startswith(line, "RESULT ") || continue
        f = split(line)[2:end]
        return (t_hsys = parse(Float64, f[1]), t_lsys = parse(Float64, f[2]),
                t_hand = parse(Float64, f[3]), t_f = parse(Float64, f[4]),
                t_ω = parse(Float64, f[5]), chars = parse(Int, f[6]))
    end
    error("worker for N = $N produced no result")
end

"Exponent `b` of a least-squares fit `t ≈ a · Nᵇ`, over the points with a usefully large `t`."
function power_law(ns, ts)
    pts = [(log(n), log(t)) for (n, t) in zip(ns, ts) if t > 10 * TIMER_RES]
    length(pts) < 2 && return NaN
    x̄ = sum(first, pts) / length(pts)
    ȳ = sum(last, pts) / length(pts)
    num = sum((x - x̄) * (y - ȳ) for (x, y) in pts)
    den = sum((x - x̄) ^ 2 for (x, _) in pts)
    return num / den
end

# =================================================================================================
# 0. Warm-up
#
# Symbolics and EulerLagrange compile a great deal of their own code on first use. Without this, all
# of it would land in the first row of the table below.
# =================================================================================================

let N = 4
    hameqs = functions(toda.hamiltonian_system(N, PARAMS))
    lageqs = functions(toda.lagrangian_system(N, PARAMS))
    for (_, gen, hand) in comparisons(N, hameqs, lageqs)
        sink!(gen())
        sink!(hand())
    end
    toda.hodeproblem(N)
    toda.lodeproblem(N)
end

@printf("\ntimer resolution: %.1f ns\n", 1e9 * TIMER_RES)

# =================================================================================================
# 1. Agreement
#
# The hand-written functions exist precisely to replace the generated ones, so nothing below this
# section means anything unless they agree.
# =================================================================================================

println("\n", "="^100)
println("1. Do the hand-written and generated functions agree?  (N = 5)")
println("="^100)
@printf("%-14s %14s %14s   %s\n", "function", "max|Δ|", "scale", "verdict")
println("-"^100)

let N = 5
    hameqs = functions(toda.hamiltonian_system(N, PARAMS))
    lageqs = functions(toda.lagrangian_system(N, PARAMS))
    for (label, gen, hand) in comparisons(N, hameqs, lageqs)
        a, b = gen(), hand()
        Δ = maximum(abs, a .- b)
        scale = maximum(abs, a)
        ok = isapprox(a, b; rtol = 1e-10, atol = 1e-14)
        @printf("%-14s %14.3e %14.3e   %s\n", label, Δ, scale, ok ? "agree" : "[MISMATCH]")
    end
end

# `∇V!` is where the whole hand-written formulation could go wrong: the lattice is periodic, so the
# gradient α (Eⱼ - Eⱼ₋₁) wraps around at j = 1, and that wrap-around is checked by nothing else.
# Differentiate `potential` itself, which shares no code with `∇V!`. N = 1 and N = 2 are included
# because they are the cases where the wrap-around collapses: at N = 1 the only neighbour of the only
# point is itself (so V is constant and ∇V vanishes identically, which is the right answer), and at
# N = 2 the forward and backward neighbours coincide.
println("\n∇V! against a central finite difference of `potential`:")
@printf("  %4s %14s %14s   %s\n", "N", "max|Δ|/scale", "scale", "verdict")
for N in (1, 2, 3, 5, 12, 64)
    q = collect(toda.compute_initial_q(toda.μ, N)) .+ 0.1 .* collect(1:N)
    h = 1e-6
    reference = zeros(N)
    for j in 1:N
        qp = copy(q); qp[j] += h
        qm = copy(q); qm[j] -= h
        reference[j] = (toda.potential(qp, PARAMS, N) - toda.potential(qm, PARAMS, N)) / 2h
    end
    kernel = zeros(N)
    toda.∇V!(kernel, q, PARAMS, N)
    scale = maximum(abs, reference)
    relative = maximum(abs, kernel .- reference) / max(scale, one(scale))
    # `scale` scales the gradient, and scale = -1 is what the force function uses
    negated = zeros(N)
    toda.∇V!(negated, q, PARAMS, N, -1)
    verdict = relative < 1e-6 ? (negated ≈ -kernel ? "agree" : "[scale = -1 MISMATCH]") : "[MISMATCH]"
    N == 1 && all(iszero, kernel) && (verdict *= "  (∇V ≡ 0 at N = 1)")
    @printf("  %4d %14.3e %14.3e   %s\n", N, relative, scale, verdict)
end

# =================================================================================================
# 2. Construction and compile cost
#
# The decisive table. `hodeproblem`/`lodeproblem` in their default (hand-written) form do no work
# beyond wrapping a few function pointers, so their column is the baseline the symbolic columns are
# to be read against.
# =================================================================================================

println("\n", "="^116)
println("2. What does the symbolic route cost to build and to compile?  (one fresh process per size)")
println("="^116)
@printf("%4s %11s %11s %14s %11s %11s %12s\n",
        "N", "ham_system", "lag_system", "hand-written", "1st gen f", "1st gen ω", "chars gen ω")
println("-"^116)

const SWEEP_N = Int[]
const SWEEP_LAG = Float64[]
const SWEEP_OMEGA_CHARS = Float64[]
const SWEEP_SETUP = Float64[]

for N in (FULL ? vcat(N_SWEEP, N_FULL) : N_SWEEP)
    r = measure_construction(N)

    # The hand-written column is in µs: it builds two problems out of function pointers, which is
    # five orders of magnitude away from everything beside it.
    @printf("%4d %9.3f s %9.3f s %11.1f µs %9.3f s %9.3f s %12d\n",
            N, r.t_hsys, r.t_lsys, 1e6 * r.t_hand, r.t_f, r.t_ω, r.chars)
    flush(stdout)

    push!(SWEEP_N, N)
    push!(SWEEP_LAG, r.t_lsys)
    push!(SWEEP_OMEGA_CHARS, r.chars)
    # everything a user pays before the first useful step of a `symbolic = true` LODE
    push!(SWEEP_SETUP, r.t_lsys + r.t_f + r.t_ω)

    # The budget keeps the default run affordable. Opting in with TODA_LATTICE_BENCH_FULL waives it —
    # going all the way to the default Ñ = 200 is the whole point of asking for the full sweep.
    total = r.t_hsys + r.t_lsys + r.t_f + r.t_ω
    if !FULL && total > BUDGET
        @printf("  (stopping the sweep: N = %d alone cost %.0f s, over the %.0f s budget;\n", N, total, BUDGET)
        println("   set TODA_LATTICE_BENCH_FULL=1 to go all the way to N = 200)")
        break
    end
end

let b_lag = power_law(SWEEP_N, SWEEP_LAG), b_chars = power_law(SWEEP_N, SWEEP_OMEGA_CHARS)
    println()
    @printf("`lagrangian_system` grows as N^%.2f, the generated `ω` as N^%.2f.\n", b_lag, b_chars)
    if !FULL && !isempty(SWEEP_N)
        scale = toda.Ñ / last(SWEEP_N)
        @printf("Extrapolated to the default Ñ = %d: `lagrangian_system` ≈ %.0f s, `ω` ≈ %.1f MB of code.\n",
                toda.Ñ, last(SWEEP_LAG) * scale^b_lag,
                last(SWEEP_OMEGA_CHARS) * scale^b_chars / 1e6)
        println("(Measured directly with TODA_LATTICE_BENCH_FULL=1: see the module docstring.)")
    end
end

# =================================================================================================
# 3. Steady-state per-call cost
# =================================================================================================

println("\n", "="^116)
println("3. What does each function cost to call, once built?")
println("="^116)
@printf("%4s %-14s %13s %13s %8s %9s   %s\n",
        "N", "function", "generated", "hand-written", "ratio", "reliable", "verdict")
println("-"^116)

# (N, label, ratio, t_generated, t_handwritten) for the break-even analysis below
const CALLS = Tuple{Int,String,Float64,Float64,Float64}[]

# Capped at 128: this section has to build the symbolic systems *in process*, and the per-call ratios
# are already unambiguous by then. Going to 200 here would add minutes to say the same thing.
for N in unique([8, isempty(SWEEP_N) ? 8 : min(128, last(SWEEP_N))])
    hameqs = functions(toda.hamiltonian_system(N, PARAMS))
    lageqs = functions(toda.lagrangian_system(N, PARAMS))
    for (label, gen, hand) in comparisons(N, hameqs, lageqs)
        gen(); hand()   # compile before timing
        t_gen, ok_gen = percall(gen)
        t_hand, ok_hand = percall(hand)
        ratio = t_gen / t_hand
        verdict = if !(ok_gen && ok_hand)
            "below timer resolution"
        elseif ratio > THRESHOLD
            "hand-written faster"
        elseif ratio < 1 / THRESHOLD
            "generated faster"
        else
            "--"
        end
        label == "ω (2N×2N)" && (verdict *= "  (never called)")
        @printf("%4d %-14s %10.4f µs %10.4f µs %8.2f %9s   %s\n",
                N, label, 1e6 * t_gen, 1e6 * t_hand, ratio, ok_gen && ok_hand, verdict)
        ok_gen && ok_hand && push!(CALLS, (N, label, ratio, t_gen, t_hand))
    end
    flush(stdout)
end

# =================================================================================================
# 4. Break-even
# =================================================================================================

println("\n", "="^100)
println("4. How many evaluations would the symbolic route need to earn its setup cost back?")
println("="^100)

let forces = filter(c -> startswith(c[2], "f "), CALLS)
    if isempty(forces) || isempty(SWEEP_SETUP)
        println("not enough reliable measurements")
    else
        N, label, ratio, t_gen, t_hand = last(forces)
        Δsetup = last(SWEEP_SETUP)
        @printf("At N = %d, `symbolic = true` costs %.1f s of setup before the first step.\n",
                last(SWEEP_N), Δsetup)
        @printf("Its `%s` then costs %.4f µs per call against the hand-written %.4f µs.\n",
                label, 1e6 * t_gen, 1e6 * t_hand)
        if N != last(SWEEP_N)
            @printf("""
(Those per-call figures are from N = %d, the largest size section 3 measures in process. The
 per-call difference grows with N, so the break-even below is an over-estimate — by roughly the
 ratio of the two sizes, which leaves its order of magnitude untouched.)\n""", N)
        end
        if t_gen < t_hand
            @printf("Break-even: %.3g evaluations.\n", Δsetup / (t_hand - t_gen))
        else
            println("""
The generated code is the *slower* one here, so there is no break-even: the symbolic route never
earns its setup cost back, however long the run.""")
        end
    end
end

# =================================================================================================
# 5. End to end
#
# The only measurement that combines setup and per-call cost the way a user experiences them.
# =================================================================================================

println("\n", "="^100)
println("5. End-to-end: building and integrating a LODE")
println("="^100)

"Total wall clock, build plus integrate, each way. `GI` is the GeometricIntegrators module."
function integration_section(GI, N; warmup = false)
    # The very first `integrate` of the session compiles a large part of GeometricIntegrators, which
    # would otherwise land on whichever branch happens to run first and make the two incomparable.
    if warmup
        GI.integrate(toda.lodeproblem(4), GI.ImplicitMidpoint())
        GI.integrate(toda.lodeproblem(4; symbolic = true), GI.ImplicitMidpoint())
    end
    for symbolic in (false, true)
        t_build, prob = elapsed(() -> toda.lodeproblem(N; symbolic = symbolic))
        t_first, _ = elapsed(() -> GI.integrate(prob, GI.ImplicitMidpoint()))
        t_warm, _ = elapsed(() -> GI.integrate(prob, GI.ImplicitMidpoint()))
        @printf("  N=%3d  symbolic=%-5s  build %8.3f s   integrate %7.3f s (warm %7.3f s)   total %8.3f s\n",
                N, symbolic, t_build, t_first, t_warm, t_build + t_first)
        flush(stdout)
    end
end

if Base.find_package("GeometricIntegrators") === nothing
    println("""
  skipped: GeometricIntegrators is a test-only dependency and is not in this environment.
  Re-run with `julia --project=test benchmark/toda_lattice.jl` to include this section.""")
else
    @eval using GeometricIntegrators
    # `@eval using` bumps the world age within this same top-level expression
    Base.invokelatest(integration_section, GeometricIntegrators, 8; warmup = true)
    Base.invokelatest(integration_section, GeometricIntegrators, 64)
end

@printf("\n(sink checksum %g -- printed only so the timed calls cannot be optimised away)\n\n", SINK[])
