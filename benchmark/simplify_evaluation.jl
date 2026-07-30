# Does `simplify=true` or `simplify=false` give faster-to-*evaluate* equations?
#
#     julia --project=. benchmark/simplify_evaluation.jl
#
# `simplify` is a keyword of `EulerLagrange.LagrangianSystem`/`HamiltonianSystem`. Its effect on
# *construction* time is already known and large (see `benchmark/outer_solar_system.jl`). The
# question here is the opposite one: given that the equations have been built, is the generated code
# cheaper to call one way or the other? Construction happens once, evaluation happens once per
# nonlinear-solver iteration per timestep, so if `simplify=true` produced meaningfully faster code
# it could still be worth paying for.
#
# For every problem in this package whose equations come from EulerLagrange, this builds the system
# both ways and times the functions an integrator actually calls — `ϑ, f, ω, L` for a Lagrangian
# system, `v, f, H` for a Hamiltonian one — reporting seconds per call.

using EulerLagrange
using GeometricProblems
using LinearAlgebra
using Printf

import GeometricProblems.CoupledHarmonicOscillator as cho
import GeometricProblems.DoublePendulum as dp
import GeometricProblems.DuffingOscillator as duff
import GeometricProblems.HenonHeilesPotential as hh
import GeometricProblems.LennardJonesOscillator as lj
import GeometricProblems.LinearWave as lw
import GeometricProblems.LotkaVolterra4dLagrangian as lv4l
import GeometricProblems.MathewsLakshmananOscillator as ml
import GeometricProblems.MorseOscillator as morse
import GeometricProblems.OuterSolarSystem as oss
import GeometricProblems.PerturbedPendulum as pp
import GeometricProblems.ThreeBody as tb
import GeometricProblems.TodaLattice as toda

include("timing.jl")

# ------------------------------------------------------------------------------------------------
# Builders. Every module wraps `LagrangianSystem`/`HamiltonianSystem` the same way; `TodaLattice`
# and `LinearWave` take the lattice size as an extra argument to their Lagrangian.

build_lagrangian(mod, d, params, extra...) = simplify -> begin
    t, x, v = lagrangian_variables(d)
    sp = symbolize(params)
    LagrangianSystem(mod.lagrangian(t, x, v, sp, extra...), t, x, v, sp; simplify=simplify)
end

build_hamiltonian(mod, d, params, extra...) = simplify -> begin
    t, q, p = hamiltonian_variables(d)
    sp = symbolize(params)
    HamiltonianSystem(mod.hamiltonian(t, q, p, sp, extra...), t, q, p, sp; simplify=simplify)
end

# `LotkaVolterra4dLagrangian` is the one `DegenerateLagrangianSystem` here, built from a kinetic and
# a Hamiltonian part rather than from a single Lagrangian.
build_degenerate(params) = simplify -> begin
    t, x, v = lagrangian_variables(4)
    sp = symbolize(params)
    Ks = lv4l.K(x, v, lv4l.A_default, lv4l.B_default)
    Hs = lv4l.H(x, lv4l.get_parameters(sp)...)
    DegenerateLagrangianSystem(Ks, Hs, t, x, v, sp; simplify=simplify)
end

# `simplify=false` has been the EulerLagrange default since 0.5, and these two are the only problems
# here whose size is a free parameter. Both are measured on a small lattice, because construction
# cost grows steeply with it — see `benchmark/linear_wave.jl` and `benchmark/toda_lattice.jl` for how
# steeply.
const TODA_N = 8
const WAVE_N = 8

# Admissible states are taken from each problem's own initial conditions, since several of these
# have restricted domains (Lennard-Jones and Morse are singular or flat at the wrong radius,
# Mathews-Lakshmanan carries a 1/(1+λq²) factor).
ics(prob) = (collect(prob.ics.q), collect(prob.ics.p))

const toda_q, toda_p = ics(toda.hodeproblem(TODA_N))
const wave_q, wave_p = ics(lw.hodeproblem(WAVE_N))

# (label, kind, dimension, builder, params, q, w)  --  `w` is the momentum for :ham, velocity for :lag
const PROBLEMS = Any[
    ("CoupledHarmonicOscillator", :lag, 2, build_lagrangian(cho, 2, cho.default_parameters()), cho.default_parameters(), [0.5, 0.7], [0.3, -0.2]),
    ("CoupledHarmonicOscillator", :ham, 2, build_hamiltonian(cho, 2, cho.default_parameters()), cho.default_parameters(), [0.5, 0.7], [0.3, -0.2]),
    ("DoublePendulum", :lag, 2, build_lagrangian(dp, 2, dp.default_parameters()), dp.default_parameters(), [0.4, 0.6], [0.3, -0.2]),
    ("DoublePendulum", :ham, 2, build_hamiltonian(dp, 2, dp.default_parameters()), dp.default_parameters(), [0.4, 0.6], [0.3, -0.2]),
    ("DuffingOscillator", :lag, 1, build_lagrangian(duff, 1, duff.default_parameters()), duff.default_parameters(), [0.7], [0.3]),
    ("DuffingOscillator", :ham, 1, build_hamiltonian(duff, 1, duff.default_parameters()), duff.default_parameters(), [0.7], [0.3]),
    ("HenonHeilesPotential", :lag, 2, build_lagrangian(hh, 2, hh.default_parameters()), hh.default_parameters(), [0.2, 0.3], [0.1, -0.1]),
    ("HenonHeilesPotential", :ham, 2, build_hamiltonian(hh, 2, hh.default_parameters()), hh.default_parameters(), [0.2, 0.3], [0.1, -0.1]),
    ("LennardJonesOscillator", :lag, 1, build_lagrangian(lj, 1, lj.default_parameters()), lj.default_parameters(), [1.2], [0.3]),
    ("LennardJonesOscillator", :ham, 1, build_hamiltonian(lj, 1, lj.default_parameters()), lj.default_parameters(), [1.2], [0.3]),
    ("MathewsLakshmananOsc.", :lag, 1, build_lagrangian(ml, 1, ml.default_parameters()), ml.default_parameters(), [0.4], [0.3]),
    ("MathewsLakshmananOsc.", :ham, 1, build_hamiltonian(ml, 1, ml.default_parameters()), ml.default_parameters(), [0.4], [0.3]),
    ("MorseOscillator", :lag, 1, build_lagrangian(morse, 1, morse.default_parameters()), morse.default_parameters(), [1.1], [0.3]),
    ("MorseOscillator", :ham, 1, build_hamiltonian(morse, 1, morse.default_parameters()), morse.default_parameters(), [1.1], [0.3]),
    ("PerturbedPendulum", :lag, 1, build_lagrangian(pp, 1, pp.default_parameters()), pp.default_parameters(), [0.5], [0.3]),
    ("PerturbedPendulum", :ham, 1, build_hamiltonian(pp, 1, pp.default_parameters()), pp.default_parameters(), [0.5], [0.3]),
    ("ThreeBody", :lag, 6, build_lagrangian(tb, 6, tb.default_parameters()), tb.default_parameters(), tb.initial_condition.q, tb.initial_condition.p),
    ("ThreeBody", :ham, 6, build_hamiltonian(tb, 6, tb.default_parameters()), tb.default_parameters(), tb.initial_condition.q, tb.initial_condition.p),
    ("TodaLattice (N=$TODA_N)", :lag, TODA_N, build_lagrangian(toda, TODA_N, toda.default_parameters(), TODA_N), toda.default_parameters(), toda_q, toda_p .+ 0.1),
    ("TodaLattice (N=$TODA_N)", :ham, TODA_N, build_hamiltonian(toda, TODA_N, toda.default_parameters(), TODA_N), toda.default_parameters(), toda_q, toda_p .+ 0.1),
    ("LinearWave (N=$WAVE_N)", :lag, WAVE_N + 2, build_lagrangian(lw, WAVE_N + 2, lw.default_parameters(), WAVE_N), lw.default_parameters(), wave_q, wave_p),
    ("LinearWave (N=$WAVE_N)", :ham, WAVE_N + 2, build_hamiltonian(lw, WAVE_N + 2, lw.default_parameters(), WAVE_N), lw.default_parameters(), wave_q, wave_p),
    ("LotkaVolterra4d (degen.)", :deg, 4, build_degenerate(lv4l.default_parameters()), lv4l.default_parameters(), lv4l.q₀, [0.3, -0.2, 0.1, 0.4]),
    ("OuterSolarSystem", :lag, 18, build_lagrangian(oss, 18, oss.default_parameters()), oss.default_parameters(), oss.q₀, oss.p₀),
    ("OuterSolarSystem", :ham, 18, build_hamiltonian(oss, 18, oss.default_parameters()), oss.default_parameters(), oss.q₀, oss.p₀),
]

"""
    callables(kind, eqs, d, q, w, params)

The generated functions an integrator actually calls, as (name, thunk) pairs.

Each thunk returns the output it wrote, so that `percall` can feed it to `sink!` and the call cannot
be optimised away.
"""
function callables(kind, eqs, d, q, w, params)
    t = 0.0
    out = zeros(d)
    if kind === :ham
        return [("v", () -> (eqs.v(out, t, q, w, params); out)),
            ("f", () -> (eqs.f(out, t, q, w, params); out)),
            ("H", () -> eqs.H(t, q, w, params))]
    end
    # A regular Lagrangian carries a 2n×2n two-form, a degenerate one an n×n form.
    Ω = kind === :deg ? zeros(d, d) : zeros(2d, 2d)
    return [("ϑ", () -> (eqs.ϑ(out, t, q, w, params); out)),
        ("f", () -> (eqs.f(out, t, q, w, params); out)),
        ("ω", () -> (eqs.ω(Ω, t, q, w, params); Ω)),
        ("L", () -> eqs.L(t, q, w, params))]
end

# ------------------------------------------------------------------------------------------------

@printf("\ntimer resolution: %.1f ns\n", 1e9 * TIMER_RES)
println("\n", "="^112)
println("Evaluation cost of the generated equations: simplify=true vs simplify=false")
println("="^112)
@printf("%-26s %-5s %-3s %11s %11s %8s %6s   %s\n",
    "problem", "kind", "fn", "simp=true", "simp=false", "ratio", "chars", "verdict")
println("-"^112)

# A 5% difference in a hand-rolled microbenchmark is not a real difference. Only call something a
# win when it is well clear of that.
const THRESHOLD = 1.15

# (label, kind, fn, ratio) for the summary; only reliable measurements are recorded
const RATIOS = Tuple{String,String,String,Float64}[]
const BUILD = Tuple{String,String,Float64,Float64}[]
const SKIPPED = String[]

for (label, kind, d, builder, params, q, w) in PROBLEMS
    systems = Dict{Bool,Any}()
    times = Dict{Bool,Float64}()
    failed = false
    for simplify in (true, false)
        try
            start = time()
            sys = builder(simplify)
            times[simplify] = time() - start
            systems[simplify] = functions(sys)
        catch err
            @printf("%-26s %-5s     BUILD FAILED (simplify=%s): %s\n",
                label, string(kind), simplify, first(split(sprint(showerror, err), "\n")))
            failed = true
        end
    end
    failed && continue
    push!(BUILD, (label, string(kind), times[true], times[false]))

    with = callables(kind, systems[true], d, q, w, params)
    without = callables(kind, systems[false], d, q, w, params)

    for ((fname, fyes), (_, fno)) in zip(with, without)
        ryes, rno = fyes(), fno()   # also warms the compiler
        agree = isapprox(ryes, rno; rtol=1e-8, atol=1e-14)
        tyes, ok_yes = percall(fyes)
        tno, ok_no = percall(fno)
        ratio = tyes / tno

        # size of the emitted code, as a noise-free cross-check on the timings
        nchars = length(string(getproperty(systems[true], Symbol(fname)).body)) /
                 length(string(getproperty(systems[false], Symbol(fname)).body))

        verdict = if !(ok_yes && ok_no)
            push!(SKIPPED, "$label/$kind/$fname")
            "below timer resolution"
        elseif ratio < 1 / THRESHOLD
            "simplify=true faster"
        elseif ratio > THRESHOLD
            "simplify=false faster"
        else
            "--"
        end
        agree || (verdict *= "  [VALUE MISMATCH]")
        ok_yes && ok_no && push!(RATIOS, (label, string(kind), fname, ratio))

        @printf("%-26s %-5s %-3s %8.4f µs %8.4f µs %8.2f %6.2f   %s\n",
            label, string(kind), fname, 1e6 * tyes, 1e6 * tno, ratio, nchars, verdict)
    end
    flush(stdout)
end

println("\n", "="^112)
println("Summary")
println("="^112)
println("`ratio` and `chars` are simplify=true relative to simplify=false; > 1 means",
    "\n`simplify=true` is slower / larger. A win needs to exceed $(THRESHOLD)x.")
println()

let r = [x[4] for x in RATIOS]
    faster_true = count(<(1 / THRESHOLD), r)
    faster_false = count(>(THRESHOLD), r)
    neutral = length(r) - faster_true - faster_false
    @printf("%d of %d comparisons were above the timer resolution (%d skipped)\n",
        length(r), length(r) + length(SKIPPED), length(SKIPPED))
    @printf("  simplify=true faster:  %3d\n", faster_true)
    @printf("  simplify=false faster: %3d\n", faster_false)
    @printf("  no difference:         %3d\n", neutral)
    if !isempty(r)
        s = sort(r)
        @printf("  median ratio: %.3f   min %.3f   max %.3f\n",
            s[cld(length(s), 2)], first(s), last(s))
    end
    println()
    println("Construction time, summed over every problem:")
    @printf("  simplify=true:  %8.2f s\n", sum(x[3] for x in BUILD))
    @printf("  simplify=false: %8.2f s\n", sum(x[4] for x in BUILD))
    println()
    @printf("(sink checksum %g -- printed only so the timed calls cannot be optimised away)\n", SINK[])
end

println()
