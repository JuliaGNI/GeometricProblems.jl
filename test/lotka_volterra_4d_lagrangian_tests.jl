
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


@testset "$(rpad("Lotka-Volterra 4D (Lagrangian)",80))" begin
    lag_ode  = LotkaVolterra4dLagrangian.lotka_volterra_4d_ode()
    lag_iode = LotkaVolterra4dLagrangian.lotka_volterra_4d_iode()
    lag_idae = LotkaVolterra4dLagrangian.lotka_volterra_4d_idae()

    ref_ode  = LotkaVolterra4d.lotka_volterra_4d_ode()
    ref_iode = LotkaVolterra4d.lotka_volterra_4d_iode()
    ref_idae = LotkaVolterra4d.lotka_volterra_4d_idae()

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
    lag_equs[:v](tbegin(lag_ode), lag_ode.ics.q, v1)
    ref_equs[:v](tbegin(ref_ode), ref_ode.ics.q, v2)
    @test v1 ≈ v2  atol=1E-14

    h1 = lag_invs[:h](tbegin(lag_ode), lag_ode.ics.q)
    h2 = ref_invs[:h](tbegin(ref_ode), ref_ode.ics.q)
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(lag_ode, TableauGauss(2))
    sol = integrate(lag_ode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_ode, TableauGauss(2))
    ref_sol = integrate(ref_ode, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 2E-14


    lag_equs = functions(lag_iode)
    lag_invs = invariants(lag_iode)
    ref_equs = functions(ref_iode)
    ref_invs = invariants(ref_iode)

    @test lag_iode.ics.p == ref_iode.ics.p

    ϑ1 = zero(lag_iode.ics.q)
    ϑ2 = zero(ref_iode.ics.q)
    lag_equs[:ϑ](tbegin(lag_iode), lag_iode.ics.q, v1, ϑ1)
    ref_equs[:ϑ](tbegin(ref_iode), ref_iode.ics.q, v2, ϑ2)
    @test ϑ1 ≈ ϑ2  atol=1E-14

    f1 = zero(lag_iode.ics.q)
    f2 = zero(ref_iode.ics.q)
    lag_equs[:f](tbegin(lag_iode), lag_iode.ics.q, v1, f1)
    ref_equs[:f](tbegin(ref_iode), ref_iode.ics.q, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    g1 = zero(lag_iode.ics.q)
    g2 = zero(ref_iode.ics.q)
    lag_equs[:g](tbegin(lag_iode), lag_iode.ics.q, v1, g1)
    ref_equs[:g](tbegin(ref_iode), ref_iode.ics.q, v2, g2)
    @test g1 ≈ g2  atol=1E-14

    v1 = zero(lag_iode.ics.q)
    v2 = zero(ref_iode.ics.q)
    lag_equs[:v̄](tbegin(lag_iode), lag_iode.ics.q, v1)
    ref_equs[:v̄](tbegin(ref_iode), ref_iode.ics.q, v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(lag_iode.ics.q)
    f2 = zero(ref_iode.ics.q)
    lag_equs[:f̄](tbegin(lag_iode), lag_iode.ics.q, v1, f1)
    ref_equs[:f̄](tbegin(ref_iode), ref_iode.ics.q, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = lag_invs[:h](tbegin(lag_iode), lag_iode.ics.q, zero(lag_iode.ics.q))
    h2 = ref_invs[:h](tbegin(ref_iode), ref_iode.ics.q, zero(ref_iode.ics.q))
    @test h1 ≈ h2  atol=1E-14

    int = IntegratorVPRKpMidpoint(lag_iode, TableauVPGLRK(2))
    sol = integrate(lag_iode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    ref_int = IntegratorVPRKpMidpoint(ref_iode, TableauVPGLRK(2))
    ref_sol = integrate(ref_iode, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 8E-14

    int = IntegratorVPRKpSymmetric(lag_iode, TableauVPGLRK(2))
    sol = integrate(lag_iode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    ref_int = IntegratorVPRKpSymmetric(ref_iode, TableauVPGLRK(2))
    ref_sol = integrate(ref_iode, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 4E-14


    lag_equs = functions(lag_idae)
    lag_invs = invariants(lag_idae)
    ref_equs = functions(ref_idae)
    ref_invs = invariants(ref_idae)

    @test lag_idae.ics.p == ref_idae.ics.p

    v1 = zero(lag_idae.ics.q)
    v2 = zero(ref_idae.ics.q)
    lag_equs[:v̄](tbegin(lag_idae), lag_idae.ics.q, v1)
    ref_equs[:v̄](tbegin(ref_idae), ref_idae.ics.q, v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(lag_idae.ics.q)
    f2 = zero(ref_idae.ics.q)
    lag_equs[:f̄](tbegin(lag_idae), lag_idae.ics.q, v1, f1)
    ref_equs[:f̄](tbegin(ref_idae), ref_idae.ics.q, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = lag_invs[:h](tbegin(lag_idae), lag_idae.ics.q, zero(lag_idae.ics.q))
    h2 = ref_invs[:h](tbegin(ref_idae), ref_idae.ics.q, zero(ref_idae.ics.q))
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(lag_idae, TableauVSPARKGLRKpMidpoint(2))
    sol = integrate(lag_idae, int)
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpMidpoint(2))
    ref_sol = integrate(ref_idae, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 4E-14

    int = Integrator(lag_idae, TableauVSPARKGLRKpSymmetric(2))
    sol = integrate(lag_idae, int)
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpSymmetric(2))
    ref_sol = integrate(ref_idae, ref_int)
    @test relative_maximum_error(sol.q, ref_sol.q[end]) < 8E-14

end
