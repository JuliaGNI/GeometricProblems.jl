using Test
using GeometricIntegrators
using PoincareInvariants

import GeometricProblems
import GeometricProblems.LotkaVolterra2d as lv2
import GeometricProblems.LotkaVolterra2dGauge as lv2g
import GeometricProblems.LotkaVolterra2dSingular as lv2si
import GeometricProblems.LotkaVolterra2dSymmetric as lv2sy
import GeometricProblems.MasslessChargedParticle as mcp
import GeometricProblems.MasslessChargedParticleSingular as mcpsi

# Loading PoincareInvariants activates the `LotkaVolterra2dPoincareInvariants` and
# `MasslessChargedParticlePoincareInvariants` extensions, which give the six degenerate Lagrangian
# problems their `poincare_invariant_1st`/`poincare_invariant_2nd` constructors.
#
# The tests below cover the wiring, which is where these can plausibly break: that both extensions
# precompile and that Pkg attaches them, that they reach every parent module, that each invariant
# picks up its own module's one- or two-form, and that the whole advection path
# (`PIEnsembleProblem` → `integrate` → `compute!`) runs and conserves what a variational integrator
# is supposed to conserve. Under `Requires` none of that was checkable, and the pre-0.4 call paths
# rotted through two renames unnoticed.

const LV2D = (lv2, lv2g, lv2si, lv2sy)
const MCP  = (mcp, mcpsi)
const MODULES = (LV2D..., MCP...)

# The loop takes any number of sample points; the surface's Chebyshev plan samples at Padua points
# and rounds up to the next Padua number, of which 45 = 9·10/2 is one.
const NLOOP = 32
const NSURFACE = 45

# Four steps is enough to see an invariant move if it is going to; the point here is the wiring,
# not the long-time behaviour.
const TIMESTEP = 0.01
const TIMESPAN = (0.0, 4TIMESTEP)


@testset "$(rpad("Poincaré invariants extensions",80))" begin

    @testset "extensions are loaded" begin
        @test Base.get_extension(GeometricProblems, :LotkaVolterra2dPoincareInvariants) !== nothing
        @test Base.get_extension(GeometricProblems, :MasslessChargedParticlePoincareInvariants) !== nothing
    end

    @testset "constructors are attached to all six modules" begin
        for M in MODULES
            @test !isempty(methods(M.poincare_invariant_1st))
            @test !isempty(methods(M.poincare_invariant_2nd))
        end
    end

    # Each module has its own gauge, so each invariant has to carry *its* form, not the one of the
    # module the `@eval` loop happened to close over last.
    @testset "each invariant carries its own module's form" begin
        for M in LV2D
            @test getform(M.poincare_invariant_1st(NLOOP)) === M.lotka_volterra_2d_ϑ
            @test getform(M.poincare_invariant_2nd(NSURFACE)) === M.lotka_volterra_2d_ω
        end
        for M in MCP
            @test getform(M.poincare_invariant_1st(NLOOP)) === M.massless_charged_particle_ϑ
            @test getform(M.poincare_invariant_2nd(NSURFACE)) === M.massless_charged_particle_ω
        end
    end

    # The Lagrangians are degenerate: the loop and the surface live in the two-dimensional
    # configuration space alone, with the momentum determined by ϑ(q).
    @testset "invariants live in configuration space" begin
        for M in MODULES
            @test getdim(M.poincare_invariant_1st(NLOOP)) == 2
            @test getdim(M.poincare_invariant_2nd(NSURFACE)) == 2
            @test getpointnum(M.poincare_invariant_1st(NLOOP)) == NLOOP
            @test getpointnum(M.poincare_invariant_2nd(NSURFACE)) ≥ NSURFACE
        end
    end

    # `f_loop`, `f_surface` and `initial_conditions_loop` are deliberately *not* in the extensions:
    # they need nothing from PoincareInvariants, so they must work in a bare environment too. This
    # asserts they are live methods rather than empty stubs.
    @testset "the parameterisations live in the package, not the extensions" begin
        for M in LV2D
            q₀ = M.initial_conditions_loop(8)
            @test size(q₀) == (2, 8)
            @test all(q₀[:, i] == M.f_loop(i, 8) for i in 1:8)
        end

        for M in MODULES
            @test length(M.f_loop(0.25)) == 2
            @test length(M.f_surface(0.25, 0.75)) == 2
            # The loop closes — the Fourier plan of the first invariant needs it periodic — and
            # the surface is centred on it.
            @test M.f_loop(0.0) ≈ M.f_loop(1.0)
            @test M.f_surface(0.5, 0.5) ≈ (M.f_loop(0.0) .+ M.f_loop(0.5)) ./ 2
        end
    end

    # The second invariant of a two-dimensional system is the ω-weighted area of the surface, which
    # for the singular Lotka-Volterra gauge (ω₁₂ = 1/q₁q₂ over a rectangle) integrates in closed
    # form. This pins the whole chain — parameterisation, two-form, quadrature — against a number
    # that owes nothing to the implementation.
    @testset "second invariant matches the analytic area integral" begin
        pinv = lv2si.poincare_invariant_2nd(NSURFACE)
        points = getpoints(lv2si.f_surface, pinv)
        I = compute!(pinv, points, 0.0, lv2si.default_parameters())

        # f_surface spans x ∈ 1 ± 0.1, y ∈ 1 ± 0.15
        exact = log(1.1 / 0.9) * log(1.15 / 0.85)

        @test abs(I) ≈ exact rtol = 1E-10
    end

    # The full path the companion packages use: sample the loop or surface, advect every sample
    # point with a variational integrator, read the invariant off the ensemble solution.
    @testset "$(nameof(M)) — advection conserves the invariants" for M in MODULES
        for (pinv, init) in ((M.poincare_invariant_1st(NLOOP), M.f_loop),
                             (M.poincare_invariant_2nd(NSURFACE), M.f_surface))
            prob = M.iodeproblem(; timespan = TIMESPAN, timestep = TIMESTEP)
            sol = integrate(PIEnsembleProblem(prob, pinv, init), VPRKGauss(2);
                            f_abstol = 1E-14, f_reltol = 1E-14)
            I = compute!(pinv, sol, parameters(prob))

            @test I isa Vector
            @test length(I) == 5
            @test all(isfinite, I)
            @test !iszero(I[begin])
            # A loose bound only; the sharp statement is the order test below. The quadrature error
            # is the same at every step and cancels in the relative variation, so what is left is
            # the integrator's own error, which a method that did *not* respect the structure would
            # blow through by orders of magnitude.
            @test maximum(abs, (I .- I[begin]) ./ I[begin]) < 1E-4
        end
    end

    # These Lagrangians are degenerate, so a variational integrator preserves the noncanonical
    # two-form ω(q) not exactly but to its own order: the solution satisfies the constraint
    # p = ϑ(q) only up to the truncation error. Asserting that order is a far sharper statement
    # about the wiring than any threshold — a form belonging to the wrong gauge, or an invariant
    # built over the wrong module, would still be small but would not converge cleanly.
    @testset "invariant error converges at the order of the method" begin
        function invariant_error(Δt)
            pinv = lv2si.poincare_invariant_1st(NLOOP)
            prob = lv2si.iodeproblem(; timespan = (0.0, 4Δt), timestep = Δt)
            sol = integrate(PIEnsembleProblem(prob, pinv, lv2si.f_loop), VPRKGauss(2);
                            f_abstol = 1E-14, f_reltol = 1E-14)
            I = compute!(pinv, sol, parameters(prob))
            maximum(abs, (I .- I[begin]) ./ I[begin])
        end

        hs = [0.02, 0.01, 0.005]
        εs = invariant_error.(hs)
        orders = [log2(εs[i] / εs[i+1]) for i in 1:(length(εs)-1)]

        @test all(o -> isapprox(o, 3.0; atol = 0.2), orders)
    end

    # The pre-0.4 names are superseded by the two constructors above but kept for compatibility.
    # Pin them to the error they throw, so that a repair — or an accidental revival against an
    # interface that happens to fit — has to come past this test.
    @testset "the pre-0.4 interface is still dead" begin
        for M in LV2D
            @test_throws UndefVarError M.ode_loop(4)
            @test_throws UndefVarError M.iode_loop(4)
            @test_throws UndefVarError M.ode_poincare_invariant_1st(0.01, 4, 10, 1)
            @test_throws UndefVarError M.iode_poincare_invariant_1st(0.01, 4, 10, 1)
        end
    end

end
