using GeometricIntegrators: ImplicitMidpoint, integrate
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

# Perturb a single parameter of the module's default parameter set.
_perturbed(M, field, δ) = merge(M.default_parameters(),
    NamedTuple{(field,)}((getfield(M.default_parameters(), field) + δ,)))

# Run both an initial-condition ensemble and a parameter ensemble through a HODE and a LODE
# integration, and assert that distinct samples yield distinct trajectories. The default
# initial conditions are passed explicitly because the modules name them differently
# (e.g. `q₀`, `θ₀`, `initial_condition.q`).
function test_ensembles(M, field, δ, q₀, p₀)
    param_vec = [M.default_parameters(), _perturbed(M, field, δ)]

    # Varying parameters only.
    hsol = integrate(M.hodeensemble(; parameters = param_vec), ImplicitMidpoint())
    lsol = integrate(M.lodeensemble(; parameters = param_vec), ImplicitMidpoint())
    @test length(hsol.s) == 2
    @test length(lsol.s) == 2
    @test hsol.s[2].q.d.parent ≉ hsol.s[1].q.d.parent
    @test lsol.s[2].q.d.parent ≉ lsol.s[1].q.d.parent

    # Varying initial conditions only.
    q₀_vec = [q₀, q₀ .+ 0.1]
    p₀_vec = [p₀, p₀ .+ 0.1]
    hsol_ic = integrate(M.hodeensemble(q₀_vec, p₀_vec), ImplicitMidpoint())
    @test length(hsol_ic.s) == 2
    @test hsol_ic.s[2].q.d.parent ≉ hsol_ic.s[1].q.d.parent
end

@testset "Duffing oscillator ensembles"            begin test_ensembles(duffing, :β, 0.5, duffing.q₀, duffing.p₀) end
@testset "Morse oscillator ensembles"              begin test_ensembles(morse,   :D, 0.3, morse.q₀,   morse.p₀)   end
@testset "Lennard-Jones oscillator ensembles"      begin test_ensembles(lj,      :ε, 0.2, lj.q₀,      lj.p₀)      end
@testset "Mathews-Lakshmanan oscillator ensembles" begin test_ensembles(ml,      :λ, 0.5, ml.q₀,      ml.p₀)      end
@testset "Double pendulum ensembles"               begin test_ensembles(dp,      :g, 0.5, dp.θ₀,      dp.p₀)      end
@testset "Three-body ensembles"                    begin test_ensembles(tb,      :G, 0.2, tb.initial_condition.q, tb.initial_condition.p) end

# The linear wave is a high-dimensional system; integrating an ensemble is expensive, so
# here we only assert that the parameter- and IC-varying ensembles can be constructed.
@testset "Linear wave ensembles (construction)" begin
    param_vec = [lw.default_parameters(), merge(lw.default_parameters(), (μ = lw.default_parameters().μ + 0.1,))]
    @test lw.hodeensemble(; parameters = param_vec) !== nothing
    @test lw.lodeensemble(; parameters = param_vec) !== nothing
    @test lw.hodeensemble([lw.q₀, lw.q₀ .+ 0.01], [lw.p₀, lw.p₀]) !== nothing
end
