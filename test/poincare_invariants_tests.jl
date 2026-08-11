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
#
# What pins the *values* is the pair of geometric identities below: the second invariant of the
# singular gauge against its closed form, and `I₂` over the disc against `I₁` over the loop bounding
# it (Stokes' theorem), which is the one check that has to hold across both invariants, both forms
# and both plans at once.
#
# What does *not* pin them is any comparison between gauges. Within a family the gauges differ by a
# gauge transformation, so their one-forms differ by an exact form, whose integral over a closed loop
# vanishes: substituting another gauge's ϑ reproduces the invariant and its convergence rate to nine
# digits. And ω = -dϑ is gauge invariant, so the four Lotka-Volterra 2d two-forms are numerically the
# same function. Form provenance is therefore checked by identity (`getform(...) === ...`) rather
# than by value.

const LV2D = (lv2, lv2g, lv2si, lv2sy)
const MCP  = (mcp, mcpsi)
const MODULES = (LV2D..., MCP...)

# The loop takes any number of sample points; the surface's Chebyshev plan samples at Padua points
# and rounds up to the next Padua number, of which 45 = 9·10/2 and 91 = 13·14/2 are two.
const NLOOP = 32
const NSURFACE = 45
const NDISC = 91

# Four steps is enough to see an invariant move if it is going to; the point here is the wiring,
# not the long-time behaviour.
const TIMESTEP = 0.01
const TIMESPAN = (0.0, 4TIMESTEP)

# The disc that `f_loop` bounds, in polar form over the unit square:
# (s, t) ↦ centre + s ⋅ radii ⋅ (cos 2πt, sin 2πt). Built from the same constants `f_loop` reads,
# so it cannot describe a different loop. The map is degenerate at s = 0 and seamed at t ∈ {0, 1},
# which costs the Chebyshev plan its spectral accuracy but nothing else.
function f_disc(M)
    x0, y0 = M.loop_centre
    rx, ry = M.loop_radii

    (s, t) -> [x0 + rx * s * cos(2π*t), y0 + ry * s * sin(2π*t)]
end


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

            x0, y0 = M.loop_centre
            rx, ry = M.loop_radii

            # The loop closes — the Fourier plan of the first invariant needs it periodic — and it
            # is the curve `loop_centre`/`loop_radii` describe, which `f_disc` relies on.
            @test M.f_loop(0.0) ≈ M.f_loop(1.0)
            @test M.f_loop(0.0) ≈ [x0 + rx, y0]
            @test M.f_loop(0.25) ≈ [x0, y0 + ry]

            # The surface is concentric with the loop and contained in it. That is the whole of
            # their geometric relationship: it is *not* the region the loop bounds, so its invariant
            # is a different number — which is why `f_disc` exists for the Stokes check below.
            @test M.f_surface(0.5, 0.5) ≈ [x0, y0]
            for s in (0.0, 1.0), t in (0.0, 1.0)
                x, y = M.f_surface(s, t)
                @test ((x - x0) / rx)^2 + ((y - y0) / ry)^2 < 1
            end
        end
    end

    # The second invariant of a two-dimensional system is the ω-weighted area of the surface, which
    # for the singular Lotka-Volterra gauge (ω₁₂ = 1/q₁q₂ over a rectangle) integrates in closed
    # form. This pins the whole chain — parameterisation, two-form, quadrature — against a number
    # that owes nothing to the implementation.
    #
    # The sign is part of the assertion. This repo's convention is Ωᵢⱼ = ∂ϑᵢ/∂qⱼ - ∂ϑⱼ/∂qᵢ, i.e.
    # Ω = -dϑ, so the computed invariant is -∫∫ Ω₁₂ dA, which is what gives it the same sign as the
    # first invariant of the loop, as Stokes' theorem demands. Taking `abs` here would let an
    # orientation flip in either the plan or the two-form through.
    @testset "second invariant matches the analytic area integral" begin
        pinv = lv2si.poincare_invariant_2nd(NSURFACE)
        points = getpoints(lv2si.f_surface, pinv)
        I = compute!(pinv, points, 0.0, lv2si.default_parameters())

        # f_surface spans x ∈ 1 ± 0.1, y ∈ 1 ± 0.15, where ω₁₂ = 1/q₁q₂
        exact = log(1.1 / 0.9) * log(1.15 / 0.85)

        @test I ≈ -exact rtol = 1E-10
    end

    # Stokes' theorem: ∮ ϑ over a loop and ∫∫ ω over the surface it bounds are the same number. This
    # is the sharpest value-level statement available, because it is the only one that has to hold
    # across *both* invariants, both forms and both plans at once — an inconsistent ϑ/ω pair, or two
    # plans whose orientations disagree, fails it while every per-invariant test still passes.
    #
    # The polar map of `f_disc` is degenerate at its centre, which costs the Chebyshev plan its
    # spectral accuracy; at NDISC points the two sides agree to a few parts in 10⁶ for the
    # Lotka-Volterra gauges and a few in 10⁸ for the massless charged particles.
    @testset "$(nameof(M)) — the two invariants agree by Stokes' theorem" for M in MODULES
        par = M.default_parameters()

        pinv1 = M.poincare_invariant_1st(NLOOP)
        I₁ = compute!(pinv1, getpoints(M.f_loop, pinv1), 0.0, par)

        pinv2 = M.poincare_invariant_2nd(NDISC)
        I₂ = compute!(pinv2, getpoints(f_disc(M), pinv2), 0.0, par)

        @test !iszero(I₁)
        @test I₁ ≈ I₂ rtol = 1E-4
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
    # p = ϑ(q) only up to the truncation error. Asserting a convergence rate rather than a threshold
    # makes this a statement about the discretisation instead of about the size of the loop.
    #
    # The observed rate is 3 (2.97 and 2.99 at these step sizes). That is an empirical number and not
    # "the order of the method": `order(VPRKGauss(2))` reports 2, while the underlying Gauss
    # q-tableau is order 4. The assertion is therefore a one-sided bound — losing convergence, or
    # dropping to first or second order, fails it; a faster rate is not a regression.
    #
    # Note what this does *not* catch: substituting another gauge's ϑ reproduces these errors to nine
    # digits, since the gauges differ by an exact form. See the note at the top of the file.
    @testset "invariant error converges faster than second order" begin
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

        @test all(o -> o > 2.5, orders)
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
