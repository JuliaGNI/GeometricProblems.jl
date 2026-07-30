using GeometricEquations
using GeometricEquations: functions, parameters
using Test

import GeometricProblems.HarmonicOscillator as ho
import GeometricProblems.LinearWave as lw
import GeometricProblems.LotkaVolterra2d as lv2
import GeometricProblems.LotkaVolterra2dGauge as lv2gauge
import GeometricProblems.LotkaVolterra2dSingular as lv2sing
import GeometricProblems.LotkaVolterra2dSymmetric as lv2symm
import GeometricProblems.LotkaVolterra4d as lv4
import GeometricProblems.MasslessChargedParticle as mcp
import GeometricProblems.MasslessChargedParticleSingular as mcps
import GeometricProblems.PointVortices as pv
import GeometricProblems.PointVorticesLinear as pvl
import GeometricProblems.TodaLattice as toda

# `LODEProblem(ϑ, f, g, ω, l, …)` and `LDAEProblem(…, ϕ, ū, ḡ, ψ, ω, l, …)` take the symplectic
# matrix `ω` and the Lagrangian `l` as *adjacent positional* arguments. No integrator in
# GeometricIntegrators evaluates `ω`, and `check_methods` skips it, so swapping the two — or
# omitting the `(Ω, t, q, v, params)` method that the LODE/LDAE actually call — goes completely
# unnoticed at construction time and during every integration.
#
# Both mistakes were present in this package: `LotkaVolterra2d.ldaeproblem`/`ldaeproblem_slrk` had
# `ω` and `l` swapped, and `HarmonicOscillator.ω!` / `LotkaVolterra4d.lotka_volterra_4d_ω` lacked
# the five-argument method. This test therefore checks every LODE/LDAE constructor in the package
# at once: `ω` must be callable with a velocity slot and produce an antisymmetric matrix of the
# right size, and `l` must be callable and return a finite scalar.
#
# Any new `lodeproblem`/`ldaeproblem` belongs in the list below.
#
# `ω` also has to have the right *size*, and that depends on the nature of the system:
#
#   * a **regular** Lagrangian is a second-order system of n equations, equivalent to a first-order
#     system of 2n, so its two-form on (q, q̇) is `2n × 2n`;
#   * a **degenerate** Lagrangian is already a first-order system of n equations, so its two-form
#     is `n × n`.
#
# EulerLagrange follows the same split: `LagrangianSystem` emits a 2n×2n `ω`, and
# `DegenerateLagrangianSystem` an n×n one.

# (name, problem, degrees of freedom, regular?, a valid (q, v) at which to evaluate)
#
# The initial data has to be admissible for the model: the Lotka-Volterra Hamiltonians take
# `log(q)`, so `q > 0`, and the point-vortex Hamiltonians take the logarithm of the distance
# between the two vortices, so the two positions must differ.
#
# `LinearWave` is entered on a small lattice rather than at its default Ñ = 256: this list is a
# `const` evaluated at load time, and allocating a 516×516 `ones` per assertion would serve no
# purpose. Its state has N + 2 components, so N = 5 gives d = 7 and a 14×14 two-form.
const WAVE_N = 5
const wave_q = lw.compute_initial_condition2(lw.μ̃, WAVE_N + 2).q
const wave_v = lw.compute_initial_condition2(lw.μ̃, WAVE_N + 2).p

# `TodaLattice` is entered on a small lattice for the same reason, and its state has exactly N
# components, so N = 5 gives d = 5 and a 10×10 two-form. Its default initial momentum is zero, which
# would make `l` degenerate here, so the velocity is offset.
const TODA_N = 5
const toda_q = toda.compute_initial_q(toda.μ, TODA_N)
const toda_v = zero(toda_q) .+ 0.1

const LODE_PROBLEMS = (
    ("HarmonicOscillator.lodeproblem",              ho.lodeproblem(),               1, true,  ([0.7], [0.3])),
    ("HarmonicOscillator.ldaeproblem",              ho.ldaeproblem(),               1, true,  ([0.7], [0.3])),
    ("HarmonicOscillator.degenerate_lodeproblem",   ho.degenerate_lodeproblem(),    2, false, ([0.7, 0.3], [0.3, -0.2])),
    ("LinearWave.lodeproblem",                      lw.lodeproblem(WAVE_N),         WAVE_N + 2, true,  (wave_q, wave_v)),
    ("LotkaVolterra2d.lodeproblem",                 lv2.lodeproblem(),              2, false, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2d.ldaeproblem",                 lv2.ldaeproblem(),              2, false, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2d.ldaeproblem_slrk",            lv2.ldaeproblem_slrk(),         2, false, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2dGauge.lodeproblem",            lv2gauge.lodeproblem(),         2, false, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2dSingular.lodeproblem",         lv2sing.lodeproblem(),          2, false, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2dSymmetric.lodeproblem",        lv2symm.lodeproblem(),          2, false, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra4d.lodeproblem",                 lv4.lodeproblem(),              4, false, ([2.0, 1.0, 1.0, 1.0], [0.3, -0.2, 0.1, 0.4])),
    ("LotkaVolterra4d.ldaeproblem",                 lv4.ldaeproblem(),              4, false, ([2.0, 1.0, 1.0, 1.0], [0.3, -0.2, 0.1, 0.4])),
    ("LotkaVolterra4d.ldaeproblem_secondary",       lv4.ldaeproblem_secondary(),    4, false, ([2.0, 1.0, 1.0, 1.0], [0.3, -0.2, 0.1, 0.4])),
    ("MasslessChargedParticle.lodeproblem",         mcp.lodeproblem(),              2, false, ([1.3, 0.2], [0.7, -0.4])),
    ("MasslessChargedParticle.ldaeproblem",         mcp.ldaeproblem(),              2, false, ([1.3, 0.2], [0.7, -0.4])),
    ("MasslessChargedParticleSingular.lodeproblem", mcps.lodeproblem(),             2, false, ([1.3, 0.2], [0.7, -0.4])),
    ("MasslessChargedParticleSingular.ldaeproblem", mcps.ldaeproblem(),             2, false, ([1.3, 0.2], [0.7, -0.4])),
    ("PointVortices.lodeproblem_formal_lagrangian", pv.lodeproblem_formal_lagrangian(),  4, false, ([0.2, 0.4, 0.4, 0.2], [0.1, -0.1, 0.2, 0.3])),
    ("PointVorticesLinear.lodeproblem_formal_lagrangian", pvl.lodeproblem_formal_lagrangian(), 4, false, ([0.3, 0.0, -0.6, 0.1], [0.1, -0.1, 0.2, 0.3])),
    ("TodaLattice.lodeproblem",                     toda.lodeproblem(TODA_N),       TODA_N, true,  (toda_q, toda_v)),
)

@testset "$(rpad("LODE/LDAE ω and l wiring",80))" begin
    @testset "$(name)" for (name, prob, d, regular, (q, v)) in LODE_PROBLEMS
        equs = functions(prob)
        pars = parameters(prob)
        t = 0.0

        # 2n×2n for a regular (second-order) Lagrangian, n×n for a degenerate (first-order) one.
        n = regular ? 2d : d

        # `ω` must accept the velocity slot the LODE/LDAE evaluate it with. If `ω` and `l` are
        # swapped this raises a MethodError, because a Lagrangian takes no output matrix.
        Ω = zeros(n, n)
        @test_nowarn equs.ω(Ω, t, q, v, pars)

        # A symplectic two-form is antisymmetric, which also rules out a matrix left partly
        # unwritten because `ω` filled the wrong number of entries.
        equs.ω(Ω, t, q, v, pars)
        @test Ω ≈ -transpose(Ω) atol = 1e-12

        # A two-form that came out all zeros would satisfy antisymmetry vacuously. Every problem
        # here has a nondegenerate two-form, so require it to be nonzero — this is what would have
        # caught a `ω` left as an identically vanishing matrix.
        @test any(!iszero, Ω)

        # `l` must be the Lagrangian: callable with (t, q, v, params) and scalar-valued. If the two
        # were swapped this returns `nothing` (an in-place `ω` writes and returns nothing).
        l = equs.l(t, q, v, pars)
        @test l isa Number
        @test isfinite(l)
    end
end
