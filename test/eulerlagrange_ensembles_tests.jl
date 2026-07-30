using GeometricIntegrators: ImplicitMidpoint, integrate
using Logging: Warn, with_logger
using Test

# Tests for issue #64: parameter-varying (and initial-condition-varying) ensembles for the
# EulerLagrange-generated HODE/LODE problems. These complement the existing ensemble tests
# for CoupledHarmonicOscillator and TodaLattice, extending coverage to the remaining
# EulerLagrange-based problems.

import GeometricProblems.DuffingOscillator as duffing
import GeometricProblems.MorseOscillator as morse
import GeometricProblems.LennardJonesOscillator as lj
import GeometricProblems.MathewsLakshmananOscillator as ml
import GeometricProblems.DoublePendulum as dp
import GeometricProblems.ThreeBody as tb
import GeometricProblems.LinearWave as lw

# Integrate and assert that the nonlinear solver stayed quiet. An implicit step whose Newton
# iteration cannot reach a root — as every ThreeBody step used to do when the default initial
# condition was still a member of the collisional `initial_conditions` grid, see
# `ThreeBody.sympnets_initial_condition` — still returns a solution, so it fails no assertion here
# while flooding the log with line-search warnings. Warnings are attributed by originating module
# rather than caught with `@test_nowarn`, both to ignore the unrelated `GeometricEquations` warning
# about single-sample ensembles and because SimpleSolvers is only an indirect dependency (via
# GeometricIntegrators), so it must not be imported.
function integrate_quietly(problem, method)
    logger = TestLogger(min_level = Warn)
    solution = with_logger(() -> integrate(problem, method), logger)
    # Assert on the distinct messages rather than their count, so that a failure reports what the
    # solver actually complained about.
    @test isempty(unique(l.message for l in logger.logs if nameof(l._module) === :SimpleSolvers))
    solution
end

# Perturb a single parameter of the module's default parameter set.
_perturbed(M, field, δ) = merge(M.default_parameters(),
    NamedTuple{(field,)}((getfield(M.default_parameters(), field) + δ,)))

# Run both an initial-condition ensemble and a parameter ensemble through a HODE and a LODE
# integration, and assert that distinct samples yield distinct trajectories. The default
# initial conditions are passed explicitly because the modules name them differently
# (e.g. `q₀`, `θ₀`, `initial_condition.q`). Every module is integrated over its own default
# timespan and timestep, which also keeps those defaults covered.
function test_ensembles(M, field, δ, q₀, p₀)
    param_vec = [M.default_parameters(), _perturbed(M, field, δ)]

    # Varying parameters only.
    hsol = integrate_quietly(M.hodeensemble(; parameters = param_vec), ImplicitMidpoint())
    lsol = integrate_quietly(M.lodeensemble(; parameters = param_vec), ImplicitMidpoint())
    @test length(hsol.s) == 2
    @test length(lsol.s) == 2
    @test hsol.s[2].q.d.parent ≉ hsol.s[1].q.d.parent
    @test lsol.s[2].q.d.parent ≉ lsol.s[1].q.d.parent

    # Varying initial conditions only.
    q₀_vec = [q₀, q₀ .+ 0.1]
    p₀_vec = [p₀, p₀ .+ 0.1]
    hsol_ic = integrate_quietly(M.hodeensemble(q₀_vec, p₀_vec), ImplicitMidpoint())
    lsol_ic = integrate_quietly(M.lodeensemble(q₀_vec, p₀_vec), ImplicitMidpoint())
    @test length(hsol_ic.s) == 2
    @test length(lsol_ic.s) == 2
    @test hsol_ic.s[2].q.d.parent ≉ hsol_ic.s[1].q.d.parent
    @test lsol_ic.s[2].q.d.parent ≉ lsol_ic.s[1].q.d.parent
end

@testset "Duffing oscillator ensembles"            begin test_ensembles(duffing, :β, 0.5, duffing.q₀, duffing.p₀) end
@testset "Morse oscillator ensembles"              begin test_ensembles(morse,   :D, 0.3, morse.q₀,   morse.p₀)   end
@testset "Lennard-Jones oscillator ensembles"      begin test_ensembles(lj,      :ε, 0.2, lj.q₀,      lj.p₀)      end
@testset "Mathews-Lakshmanan oscillator ensembles" begin test_ensembles(ml,      :λ, 0.5, ml.q₀,      ml.p₀)      end
@testset "Double pendulum ensembles"               begin test_ensembles(dp,      :g, 0.5, dp.θ₀,      dp.p₀)      end
@testset "Three-body ensembles"                    begin test_ensembles(tb,      :G, 0.2, tb.figure_eight.q, tb.figure_eight.p)          end

# `LinearWave` builds its equations of motion by hand *by default*, so `symbolic = true` is required
# here — without it these three calls would not touch EulerLagrange at all and this testset would
# silently stop testing what the file is about.
#
# It still runs on a small lattice: `lagrangian_system` grows as n^2.5 and the generated `ω` as n^2,
# so the default Ñ = 256 (258 degrees of freedom) costs minutes. Construction is all this testset
# checks, as in `toda_lattice_tests.jl`.
const N_wave = 20

@testset "Linear wave ensembles (construction)" begin
    param_vec = [lw.default_parameters(), merge(lw.default_parameters(), (μ = lw.default_parameters().μ + 0.1,))]
    q₀ = lw.compute_initial_condition2(lw.μ̃, N_wave + 2).q
    p₀ = lw.compute_initial_condition2(lw.μ̃, N_wave + 2).p
    @test lw.hodeensemble(N_wave; parameters = param_vec, symbolic = true) !== nothing
    @test lw.lodeensemble(N_wave; parameters = param_vec, symbolic = true) !== nothing
    @test lw.hodeensemble([q₀, q₀ .+ 0.01], [p₀, p₀]; symbolic = true) !== nothing
end

# The hand-written ensembles are cheap enough to actually integrate, at the same lattice size.
@testset "Linear wave ensembles (hand-written)" begin
    param_vec = [lw.default_parameters(), merge(lw.default_parameters(), (μ = lw.default_parameters().μ + 0.1,))]
    q₀ = lw.compute_initial_condition2(lw.μ̃, N_wave + 2).q
    p₀ = lw.compute_initial_condition2(lw.μ̃, N_wave + 2).p

    hsol = integrate_quietly(lw.hodeensemble(N_wave; parameters = param_vec), ImplicitMidpoint())
    lsol = integrate_quietly(lw.lodeensemble(N_wave; parameters = param_vec), ImplicitMidpoint())
    @test length(hsol.s) == 2
    @test length(lsol.s) == 2
    @test hsol.s[2].q.d.parent ≉ hsol.s[1].q.d.parent
    @test lsol.s[2].q.d.parent ≉ lsol.s[1].q.d.parent

    hsol_ic = integrate_quietly(lw.hodeensemble([q₀, q₀ .+ 0.01], [p₀, p₀]), ImplicitMidpoint())
    @test length(hsol_ic.s) == 2
    @test hsol_ic.s[2].q.d.parent ≉ hsol_ic.s[1].q.d.parent
end
