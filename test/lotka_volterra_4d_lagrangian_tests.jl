using Test
using GeometricEquations
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricSolutions

import GeometricProblems.LotkaVolterra4d
import GeometricProblems.LotkaVolterra4dLagrangian
import GeometricProblems.LotkaVolterra4dLagrangian: reference_solution


@testset "$(rpad("Lotka-Volterra 4D (Lagrangian)",80))" begin
    lag_ode = LotkaVolterra4dLagrangian.odeproblem()
    lag_lode = LotkaVolterra4dLagrangian.lodeproblem()
    lag_ldae = LotkaVolterra4dLagrangian.ldaeproblem()

    ref_ode = LotkaVolterra4d.odeproblem()
    ref_lode = LotkaVolterra4d.lodeproblem()
    ref_ldae = LotkaVolterra4d.ldaeproblem()

    ref = integrate(ref_ode, Gauss(8))

    @test initialtime(lag_ode) == initialtime(ref_ode)
    @test initialtime(lag_lode) == initialtime(ref_lode)
    @test initialtime(lag_ldae) == initialtime(ref_ldae)

    @test lag_ode.ics == ref_ode.ics
    @test lag_lode.ics == ref_lode.ics
    @test lag_ldae.ics == ref_ldae.ics


    lag_equs = functions(lag_ode)
    lag_invs = invariants(lag_ode)
    ref_equs = functions(ref_ode)
    ref_invs = invariants(ref_ode)

    v1 = zero(lag_ode.ics.q)
    v2 = zero(ref_ode.ics.q)
    lag_equs.v(v1, initialtime(lag_ode), lag_ode.ics.q, parameters(lag_ode))
    ref_equs.v(v2, initialtime(ref_ode), ref_ode.ics.q, parameters(ref_ode))
    @test v1 ≈ v2 atol = eps()

    h1 = lag_invs.h(initialtime(lag_ode), lag_ode.ics.q, parameters(lag_ode))
    h2 = ref_invs.h(initialtime(ref_ode), ref_ode.ics.q, parameters(ref_ode))
    @test h1 ≈ h2 atol = eps()

    lag_sol = integrate(lag_ode, Gauss(2))
    @test relative_maximum_error(lag_sol.q, ref.q) < 1E-8

    ref_sol = integrate(ref_ode, Gauss(2))
    @test relative_maximum_error(lag_sol.q, ref_sol.q) < 1E-14


    lag_equs = functions(lag_lode)
    lag_invs = invariants(lag_lode)
    lag_igs = initialguess(lag_lode)
    ref_equs = functions(ref_lode)
    ref_invs = invariants(ref_lode)

    @test lag_lode.ics.q == ref_lode.ics.q
    @test lag_lode.ics.p == ref_lode.ics.p
    @test parameters(lag_lode) == parameters(ref_lode)

    ϑ1 = zero(lag_lode.ics.q)
    ϑ2 = zero(ref_lode.ics.q)
    lag_equs.ϑ(ϑ1, initialtime(lag_lode), lag_lode.ics.q, v1, parameters(lag_lode))
    ref_equs.ϑ(ϑ2, initialtime(ref_lode), ref_lode.ics.q, v2, parameters(ref_lode))
    @test ϑ1 ≈ ϑ2 atol = eps()

    f1 = zero(lag_lode.ics.q)
    f2 = zero(ref_lode.ics.q)
    lag_equs.f(f1, initialtime(lag_lode), lag_lode.ics.q, v1, parameters(lag_lode))
    ref_equs.f(f2, initialtime(ref_lode), ref_lode.ics.q, v2, parameters(ref_lode))
    @test f1 ≈ f2 atol = eps()

    g1 = zero(lag_lode.ics.q)
    g2 = zero(ref_lode.ics.q)
    lag_equs.g(g1, initialtime(lag_lode), lag_lode.ics.q, v1, zero(v1), parameters(lag_lode))
    ref_equs.g(g2, initialtime(ref_lode), ref_lode.ics.q, v2, zero(v2), parameters(ref_lode))
    @test g1 ≈ g2 atol = eps()

    v1 = zero(lag_lode.ics.q)
    v2 = zero(ref_lode.ics.q)
    lag_igs.v(v1, initialtime(lag_lode), lag_lode.ics.q, zero(f1), parameters(lag_lode))
    lag_igs.v(v2, initialtime(ref_lode), ref_lode.ics.q, zero(f2), parameters(ref_lode))
    @test v1 ≈ v2 atol = eps()

    f1 = zero(lag_lode.ics.q)
    f2 = zero(ref_lode.ics.q)
    lag_igs.f(f1, initialtime(lag_lode), lag_lode.ics.q, v1, parameters(lag_lode))
    lag_igs.f(f2, initialtime(ref_lode), ref_lode.ics.q, v2, parameters(ref_lode))
    @test f1 ≈ f2 atol = eps()

    h1 = lag_invs.h(initialtime(lag_lode), lag_lode.ics.q, lag_lode.ics.v, parameters(lag_lode))
    h2 = ref_invs.h(initialtime(ref_lode), ref_lode.ics.q, ref_lode.ics.v, parameters(ref_lode))
    @test h1 ≈ h2 atol = eps()

    lag_sol = integrate(lag_lode, Gauss(2))
    @test relative_maximum_error(lag_sol.q, ref.q) < 8E-4

    ref_sol = integrate(ref_lode, Gauss(2))
    @test relative_maximum_error(ref_sol.q, ref.q) < 8E-4
    @test relative_maximum_error(ref_sol.q, lag_sol.q) < 4E-13

    lag_sol = integrate(lag_lode, MidpointProjection(VPRKGauss(2)))
    ref_sol = integrate(ref_lode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(lag_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(ref_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(ref_sol.q, lag_sol.q) < 4E-14

    lag_sol = integrate(lag_lode, SymmetricProjection(VPRKGauss(2)))
    ref_sol = integrate(ref_lode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(lag_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(ref_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(ref_sol.q, lag_sol.q) < 2E-14


    lag_equs = functions(lag_ldae)
    lag_invs = invariants(lag_ldae)
    lag_igs = initialguess(lag_ldae)
    ref_equs = functions(ref_ldae)
    ref_invs = invariants(ref_ldae)

    @test lag_ldae.ics.q == ref_ldae.ics.q
    @test lag_ldae.ics.p == ref_ldae.ics.p
    @test parameters(lag_ldae) == parameters(ref_ldae)

    ϑ1 = zero(lag_ldae.ics.q)
    ϑ2 = zero(ref_ldae.ics.q)
    lag_equs.ϑ(ϑ1, initialtime(lag_ldae), lag_ldae.ics.q, v1, parameters(lag_ldae))
    ref_equs.ϑ(ϑ2, initialtime(ref_ldae), ref_ldae.ics.q, v2, parameters(ref_ldae))
    @test ϑ1 ≈ ϑ2 atol = eps()

    f1 = zero(lag_ldae.ics.q)
    f2 = zero(ref_ldae.ics.q)
    lag_equs.f(f1, initialtime(lag_ldae), lag_ldae.ics.q, v1, parameters(lag_ldae))
    ref_equs.f(f2, initialtime(ref_ldae), ref_ldae.ics.q, v2, parameters(ref_ldae))
    @test f1 ≈ f2 atol = eps()

    g1 = zero(lag_ldae.ics.q)
    g2 = zero(ref_ldae.ics.q)
    lag_equs.g(g1, initialtime(lag_ldae), lag_ldae.ics.q, v1, lag_ldae.ics.p, lag_ldae.ics.λ, parameters(lag_ldae))
    ref_equs.g(g2, initialtime(ref_ldae), ref_ldae.ics.q, v2, ref_ldae.ics.p, ref_ldae.ics.λ, parameters(ref_ldae))
    @test g1 ≈ g2 atol = eps()

    u1 = zero(lag_ldae.ics.q)
    u2 = zero(ref_ldae.ics.q)
    lag_equs.u(u1, initialtime(lag_ldae), lag_ldae.ics.q, v1, lag_ldae.ics.p, lag_ldae.ics.λ, parameters(lag_ldae))
    ref_equs.u(u2, initialtime(ref_ldae), ref_ldae.ics.q, v2, ref_ldae.ics.p, ref_ldae.ics.λ, parameters(ref_ldae))
    @test u1 ≈ u2 atol = eps()

    ϕ1 = zero(lag_ldae.ics.q)
    ϕ2 = zero(ref_ldae.ics.q)
    lag_equs.ϕ(ϕ1, initialtime(lag_ldae), lag_ldae.ics.q, v1, lag_ldae.ics.p, parameters(lag_ldae))
    ref_equs.ϕ(ϕ2, initialtime(ref_ldae), ref_ldae.ics.q, v2, ref_ldae.ics.p, parameters(ref_ldae))
    @test ϕ1 ≈ ϕ2 atol = eps()

    v1 = zero(lag_ldae.ics.q)
    v2 = zero(ref_ldae.ics.q)
    lag_igs.v(v1, initialtime(lag_ldae), lag_ldae.ics.q, zero(f1), parameters(lag_ldae))
    lag_igs.v(v2, initialtime(ref_ldae), ref_ldae.ics.q, zero(f2), parameters(ref_ldae))
    @test v1 ≈ v2 atol = eps()

    f1 = zero(lag_ldae.ics.q)
    f2 = zero(ref_ldae.ics.q)
    lag_igs.f(f1, initialtime(lag_ldae), lag_ldae.ics.q, v1, parameters(lag_ldae))
    lag_igs.f(f2, initialtime(ref_ldae), ref_ldae.ics.q, v2, parameters(ref_ldae))
    @test f1 ≈ f2 atol = eps()

    h1 = lag_invs.h(initialtime(lag_ldae), lag_ldae.ics.q, zero(lag_ldae.ics.q), parameters(lag_ldae))
    h2 = ref_invs.h(initialtime(ref_ldae), ref_ldae.ics.q, zero(ref_ldae.ics.q), parameters(ref_ldae))
    @test h1 ≈ h2 atol = eps()

    lag_sol = integrate(lag_ldae, TableauVSPARKGLRKpMidpoint(2))
    ref_sol = integrate(ref_ldae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(lag_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(ref_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(lag_sol.q, ref_sol.q) < 1E-14

    lag_sol = integrate(lag_ldae, TableauVSPARKGLRKpSymmetric(2))
    ref_sol = integrate(ref_ldae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(lag_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(ref_sol.q, ref.q) < 1E-8
    @test relative_maximum_error(lag_sol.q, ref_sol.q) < 8E-14

end
