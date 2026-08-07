using Test
using PoincareInvariants

import GeometricProblems
import GeometricProblems.LotkaVolterra2d as lv2
import GeometricProblems.LotkaVolterra2dGauge as lv2g
import GeometricProblems.LotkaVolterra2dSingular as lv2si
import GeometricProblems.LotkaVolterra2dSymmetric as lv2sy

# Loading PoincareInvariants activates the `LotkaVolterra2dPoincareInvariants` extension. Its bodies
# are dead — they call `PoincareInvariant1st`, which PoincareInvariants 0.5 does not define, and
# `lotka_volterra_2d_ode`/`lotka_volterra_2d_iode`, which are now `odeproblem`/`iodeproblem` — and
# they stay that way until the Poincaré-invariant support gets its makeover.
#
# What is worth asserting in the meantime is that the extension still *resolves*: that it
# precompiles, that Pkg attaches it, and that it reaches all four Lotka-Volterra 2d modules. Under
# `Requires` none of that was checkable, and both call paths rotted through two renames unnoticed.
# These tests are what makes the move to an extension pay for itself; they fail on a typo in the
# `@eval` loop, on a renamed parent module, and on a `[weakdeps]`/`[extensions]` mismatch.

const MODULES = (lv2, lv2g, lv2si, lv2sy)

@testset "$(rpad("Poincaré invariants extension",80))" begin

    @testset "extension is loaded" begin
        ext = Base.get_extension(GeometricProblems, :LotkaVolterra2dPoincareInvariants)
        @test ext !== nothing
    end

    @testset "methods are attached to all four modules" begin
        for M in MODULES
            @test !isempty(methods(M.ode_loop))
            @test !isempty(methods(M.iode_loop))
            @test !isempty(methods(M.ode_poincare_invariant_1st))
            @test !isempty(methods(M.iode_poincare_invariant_1st))
        end
    end

    # `f_loop` and `initial_conditions_loop` are deliberately *not* in the extension: they need
    # nothing from PoincareInvariants, so they must work in a bare environment too. This asserts
    # they are live methods rather than empty stubs.
    @testset "loop sampling lives in the package, not the extension" begin
        for M in MODULES
            q₀ = M.initial_conditions_loop(8)
            @test size(q₀) == (2, 8)
            @test all(q₀[:, i] == M.f_loop(i, 8) for i in 1:8)
        end
    end

    # Pin the dead paths to the error they currently throw, so that a repair — or an accidental
    # revival against an interface that happens to fit — has to come past this test.
    @testset "invariant bodies are still dead" begin
        for M in MODULES
            @test_throws UndefVarError M.ode_loop(4)
            @test_throws UndefVarError M.iode_loop(4)
            @test_throws UndefVarError M.ode_poincare_invariant_1st(0.01, 4, 10, 1)
            @test_throws UndefVarError M.iode_poincare_invariant_1st(0.01, 4, 10, 1)
        end
    end

end
