using Test
using GeometricEquations
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricSolutions

import GeometricProblems.LotkaVolterra4d
import GeometricProblems.LotkaVolterra4dLagrangian
import GeometricProblems.LotkaVolterra4dLagrangian: reference_solution


@testset "$(rpad("Lotka-Volterra 4D (Lagrangian)",80))" begin
    lag_ode  = LotkaVolterra4dLagrangian.lotka_volterra_4d_ode()
    lag_iode = LotkaVolterra4dLagrangian.lotka_volterra_4d_iode()
    lag_idae = LotkaVolterra4dLagrangian.lotka_volterra_4d_idae()

    ref_ode  = LotkaVolterra4d.lotka_volterra_4d_ode()
    ref_iode = LotkaVolterra4d.lotka_volterra_4d_iode()
    ref_idae = LotkaVolterra4d.lotka_volterra_4d_idae()

    ref = integrate(ref_ode, Gauss(8))

    @assert tbegin(lag_ode)  == tbegin(ref_ode)
    @assert tbegin(lag_iode) == tbegin(ref_iode)
    @assert tbegin(lag_idae) == tbegin(ref_idae)

    @assert lag_ode.ics  == ref_ode.ics
    @assert lag_iode.ics == ref_iode.ics
    @assert lag_idae.ics == ref_idae.ics


    lag_equs = functions(lag_ode)
    lag_invs = invariants(lag_ode)
    ref_equs = functions(ref_ode)
    ref_invs = invariants(ref_ode)

    v1 = zero(lag_ode.ics.q)
    v2 = zero(ref_ode.ics.q)
    lag_equs[:v](v1, tbegin(lag_ode), lag_ode.ics.q)
    ref_equs[:v](v2, tbegin(ref_ode), ref_ode.ics.q)
    @test v1 ≈ v2  atol=eps()

    h1 = lag_invs[:h](tbegin(lag_ode), lag_ode.ics.q)
    h2 = ref_invs[:h](tbegin(ref_ode), ref_ode.ics.q)
    @test h1 ≈ h2  atol=eps()

    sol = integrate(lag_ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    ref_sol = integrate(ref_ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref_sol.q) < 2E-14


    lag_equs = functions(lag_iode)
    lag_invs = invariants(lag_iode)
    ref_equs = functions(ref_iode)
    ref_invs = invariants(ref_iode)

    @test lag_iode.ics.p == ref_iode.ics.p

    ϑ1 = zero(lag_iode.ics.q)
    ϑ2 = zero(ref_iode.ics.q)
    lag_equs[:ϑ](ϑ1, tbegin(lag_iode), lag_iode.ics.q, v1)
    ref_equs[:ϑ](ϑ2, tbegin(ref_iode), ref_iode.ics.q, v2)
    @test ϑ1 ≈ ϑ2  atol=eps()

    f1 = zero(lag_iode.ics.q)
    f2 = zero(ref_iode.ics.q)
    lag_equs[:f](f1, tbegin(lag_iode), lag_iode.ics.q, v1)
    ref_equs[:f](f2, tbegin(ref_iode), ref_iode.ics.q, v2)
    @test f1 ≈ f2  atol=eps()

    g1 = zero(lag_iode.ics.q)
    g2 = zero(ref_iode.ics.q)
    λ1 = zero(v1)
    λ2 = zero(v2)
    lag_equs[:g](g1, tbegin(lag_iode), lag_iode.ics.q, v1, λ1)
    ref_equs[:g](g2, tbegin(ref_iode), ref_iode.ics.q, v2, λ2)
    @test g1 ≈ g2  atol=eps()

    v1 = zero(lag_iode.ics.q)
    v2 = zero(ref_iode.ics.q)
    lag_equs[:v̄](v1, tbegin(lag_iode), lag_iode.ics.q, zero(f1))
    ref_equs[:v̄](v2, tbegin(ref_iode), ref_iode.ics.q, zero(f2))
    @test v1 ≈ v2  atol=eps()

    f1 = zero(lag_iode.ics.q)
    f2 = zero(ref_iode.ics.q)
    lag_equs[:f̄](f1, tbegin(lag_iode), lag_iode.ics.q, v1)
    ref_equs[:f̄](f2, tbegin(ref_iode), ref_iode.ics.q, v2)
    @test f1 ≈ f2  atol=eps()

    h1 = lag_invs[:h](tbegin(lag_iode), lag_iode.ics.q, zero(lag_iode.ics.q))
    h2 = ref_invs[:h](tbegin(ref_iode), ref_iode.ics.q, zero(ref_iode.ics.q))
    @test h1 ≈ h2  atol=eps()

    sol = integrate(lag_iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3

    ref_sol = integrate(ref_iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref_sol.q) < 8E-14

    sol = integrate(lag_iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    ref_sol = integrate(ref_iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref_sol.q) < 4E-14


    lag_equs = functions(lag_idae)
    lag_invs = invariants(lag_idae)
    ref_equs = functions(ref_idae)
    ref_invs = invariants(ref_idae)

    @test lag_idae.ics.p == ref_idae.ics.p

    v1 = zero(lag_idae.ics.q)
    v2 = zero(ref_idae.ics.q)
    lag_equs[:v̄](v1, tbegin(lag_idae), lag_idae.ics.q, zero(f1))
    ref_equs[:v̄](v2, tbegin(ref_idae), ref_idae.ics.q, zero(f2))
    @test v1 ≈ v2  atol=eps()

    f1 = zero(lag_idae.ics.q)
    f2 = zero(ref_idae.ics.q)
    lag_equs[:f̄](f1, tbegin(lag_idae), lag_idae.ics.q, v1)
    ref_equs[:f̄](f2, tbegin(ref_idae), ref_idae.ics.q, v2)
    @test f1 ≈ f2  atol=eps()

    h1 = lag_invs[:h](tbegin(lag_idae), lag_idae.ics.q, zero(lag_idae.ics.q))
    h2 = ref_invs[:h](tbegin(ref_idae), ref_idae.ics.q, zero(ref_idae.ics.q))
    @test h1 ≈ h2  atol=eps()

    sol = integrate(lag_idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3

    ref_sol = integrate(ref_idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref_sol.q) < 4E-14

    sol = integrate(lag_idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    ref_sol = integrate(ref_idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref_sol.q) < 8E-14

end
