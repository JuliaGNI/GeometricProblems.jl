
using GeometricEquations
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using SimpleSolvers
using Test

import GeometricProblems.LotkaVolterra4d
import GeometricProblems.LotkaVolterra4dLagrangian
import GeometricProblems.LotkaVolterra4dLagrangian: reference_solution


SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

SimpleSolvers.set_config(:nls_atol_break, Inf)
SimpleSolvers.set_config(:nls_rtol_break, Inf)
SimpleSolvers.set_config(:nls_stol_break, Inf)


@testset "$(rpad("Lotka-Volterra 4D",80))" begin
    ode  = LotkaVolterra4dLagrangian.lotka_volterra_4d_ode()
    iode = LotkaVolterra4dLagrangian.lotka_volterra_4d_iode()
    idae = LotkaVolterra4dLagrangian.lotka_volterra_4d_idae()

    ref_ode  = LotkaVolterra4d.lotka_volterra_4d_ode()
    ref_iode = LotkaVolterra4d.lotka_volterra_4d_iode()
    ref_idae = LotkaVolterra4d.lotka_volterra_4d_idae()

    ref_equs = functions(ref_ode)
    ref_invs = invariants(ref_ode)

    v1 = zero(ode.ics.q)
    v2 = zero(ref_ode.ics.q)
    functions(ode)[:v](tbegin(ode), ode.ics.q, v1)
    ref_equs[:v](tbegin(ref_ode), ref_ode.ics.q, v2)
    @test v1 ≈ v2  atol=1E-14

    h1 = invariants(ode)[:h](tbegin(ode), ode.ics.q)
    h2 = ref_invs[:h](tbegin(ref_ode), ref_ode.ics.q)
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(ode, TableauGauss(2))
    sol = integrate(ode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_ode, TableauGauss(2))
    ref_sol = integrate(ref_ode, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 2E-14


    ref_equs = functions(ref_iode)
    ref_invs = invariants(ref_iode)

    @test iode.p₀ == ref_iode.p₀

    ϑ1 = zero(iode.ics.q)
    ϑ2 = zero(ref_iode.ics.q)
    functions(iode)[:ϑ](tbegin(iode), iode.ics.q, v1, ϑ1)
    ref_equs[:ϑ](tbegin(ref_iode), ref_iode.ics.q, v2, ϑ2)
    @test ϑ1 ≈ ϑ2  atol=1E-14

    f1 = zero(iode.ics.q)
    f2 = zero(ref_iode.ics.q)
    functions(iode)[:f](tbegin(iode), iode.ics.q, v1, f1)
    ref_equs[:f](tbegin(ref_iode), ref_iode.ics.q, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    g1 = zero(iode.ics.q)
    g2 = zero(ref_iode.ics.q)
    functions(iode)[:g](tbegin(iode), iode.ics.q, v1, g1)
    ref_equs[:g](tbegin(ref_iode), ref_iode.ics.q, v2, g2)
    @test g1 ≈ g2  atol=1E-14

    v1 = zero(iode.ics.q)
    v2 = zero(ref_iode.ics.q)
    functions(iode)[:v̄](tbegin(iode), iode.ics.q, v1)
    ref_equs[:v̄](tbegin(ref_iode), ref_iode.ics.q, v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(iode.ics.q)
    f2 = zero(ref_iode.ics.q)
    functions(iode)[:f̄](tbegin(iode), iode.ics.q, v1, f1)
    ref_equs[:f̄](tbegin(ref_iode), ref_iode.ics.q, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = invariants(iode)[:h](tbegin(iode), iode.ics.q, zero(iode.ics.q))
    h2 = ref_invs[:h](tbegin(ref_iode), ref_iode.ics.q, zero(ref_iode.ics.q))
    @test h1 ≈ h2  atol=1E-14

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2))
    sol = integrate(iode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    ref_int = IntegratorVPRKpMidpoint(ref_iode, TableauVPGLRK(2))
    ref_sol = integrate(ref_iode, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 8E-14

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2))
    sol = integrate(iode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    ref_int = IntegratorVPRKpSymmetric(ref_iode, TableauVPGLRK(2))
    ref_sol = integrate(ref_iode, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 4E-14


    ref_equs = functions(ref_idae)
    ref_invs = invariants(ref_idae)

    @test idae.p₀ == ref_idae.p₀

    v1 = zero(idae.ics.q)
    v2 = zero(ref_idae.ics.q)
    functions(idae)[:v̄](tbegin(idae), idae.ics.q, v1)
    ref_equs[:v̄](tbegin(ref_idae), ref_idae.ics.q, v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(idae.ics.q)
    f2 = zero(ref_idae.ics.q)
    functions(idae)[:f̄](tbegin(idae), idae.ics.q, v1, f1)
    ref_equs[:f̄](tbegin(ref_idae), ref_idae.ics.q, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = invariants(idae)[:h](tbegin(idae), idae.ics.q, zero(idae.ics.q))
    h2 = ref_invs[:h](tbegin(ref_idae), ref_idae.ics.q, zero(ref_idae.ics.q))
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2))
    sol = integrate(idae, int)
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpMidpoint(2))
    ref_sol = integrate(ref_idae, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 4E-14

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2))
    sol = integrate(idae, int)
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpSymmetric(2))
    ref_sol = integrate(ref_idae, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 8E-14

end
