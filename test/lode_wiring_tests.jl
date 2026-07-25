using GeometricEquations
using GeometricEquations: functions, parameters
using Test

import GeometricProblems.HarmonicOscillator as ho
import GeometricProblems.LotkaVolterra2d as lv2
import GeometricProblems.LotkaVolterra2dGauge as lv2gauge
import GeometricProblems.LotkaVolterra2dSingular as lv2sing
import GeometricProblems.LotkaVolterra2dSymmetric as lv2symm
import GeometricProblems.LotkaVolterra4d as lv4
import GeometricProblems.MasslessChargedParticle as mcp
import GeometricProblems.MasslessChargedParticleSingular as mcps
import GeometricProblems.PointVortices as pv
import GeometricProblems.PointVorticesLinear as pvl

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

# (name, problem, degrees of freedom, a valid (q, v) at which to evaluate)
#
# The initial data has to be admissible for the model: the Lotka-Volterra Hamiltonians take
# `log(q)`, so `q > 0`, and the point-vortex Hamiltonians take the logarithm of the distance
# between the two vortices, so the two positions must differ.
const LODE_PROBLEMS = (
    ("HarmonicOscillator.lodeproblem",              ho.lodeproblem(),               1, ([0.7], [0.3])),
    ("HarmonicOscillator.ldaeproblem",              ho.ldaeproblem(),               1, ([0.7], [0.3])),
    ("HarmonicOscillator.degenerate_lodeproblem",   ho.degenerate_lodeproblem(),    2, ([0.7, 0.3], [0.3, -0.2])),
    ("LotkaVolterra2d.lodeproblem",                 lv2.lodeproblem(),              2, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2d.ldaeproblem",                 lv2.ldaeproblem(),              2, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2d.ldaeproblem_slrk",            lv2.ldaeproblem_slrk(),         2, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2dGauge.lodeproblem",            lv2gauge.lodeproblem(),         2, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2dSingular.lodeproblem",         lv2sing.lodeproblem(),          2, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra2dSymmetric.lodeproblem",        lv2symm.lodeproblem(),          2, ([2.0, 1.0], [0.3, -0.2])),
    ("LotkaVolterra4d.lodeproblem",                 lv4.lodeproblem(),              4, ([2.0, 1.0, 1.0, 1.0], [0.3, -0.2, 0.1, 0.4])),
    ("LotkaVolterra4d.ldaeproblem",                 lv4.ldaeproblem(),              4, ([2.0, 1.0, 1.0, 1.0], [0.3, -0.2, 0.1, 0.4])),
    ("LotkaVolterra4d.ldaeproblem_secondary",       lv4.ldaeproblem_secondary(),    4, ([2.0, 1.0, 1.0, 1.0], [0.3, -0.2, 0.1, 0.4])),
    ("MasslessChargedParticle.lodeproblem",         mcp.lodeproblem(),              2, ([1.3, 0.2], [0.7, -0.4])),
    ("MasslessChargedParticle.ldaeproblem",         mcp.ldaeproblem(),              2, ([1.3, 0.2], [0.7, -0.4])),
    ("MasslessChargedParticleSingular.lodeproblem", mcps.lodeproblem(),             2, ([1.3, 0.2], [0.7, -0.4])),
    ("MasslessChargedParticleSingular.ldaeproblem", mcps.ldaeproblem(),             2, ([1.3, 0.2], [0.7, -0.4])),
    ("PointVortices.lodeproblem_formal_lagrangian", pv.lodeproblem_formal_lagrangian(),  4, ([0.2, 0.4, 0.4, 0.2], [0.1, -0.1, 0.2, 0.3])),
    ("PointVorticesLinear.lodeproblem_formal_lagrangian", pvl.lodeproblem_formal_lagrangian(), 4, ([0.3, 0.0, -0.6, 0.1], [0.1, -0.1, 0.2, 0.3])),
)

@testset "$(rpad("LODE/LDAE ω and l wiring",80))" begin
    @testset "$(name)" for (name, prob, d, (q, v)) in LODE_PROBLEMS
        equs = functions(prob)
        pars = parameters(prob)
        t = 0.0

        # `ω` must accept the velocity slot the LODE/LDAE evaluate it with. If `ω` and `l` are
        # swapped this raises a MethodError, because a Lagrangian takes no output matrix.
        Ω = zeros(d, d)
        @test_nowarn equs.ω(Ω, t, q, v, pars)

        # A symplectic two-form is antisymmetric, which also rules out a matrix left partly
        # unwritten because `ω` filled the wrong number of entries.
        equs.ω(Ω, t, q, v, pars)
        @test Ω ≈ -transpose(Ω) atol = 1e-12

        # `l` must be the Lagrangian: callable with (t, q, v, params) and scalar-valued. If the two
        # were swapped this returns `nothing` (an in-place `ω` writes and returns nothing).
        l = equs.l(t, q, v, pars)
        @test l isa Number
        @test isfinite(l)
    end
end
