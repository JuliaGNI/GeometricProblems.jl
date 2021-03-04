
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using SimpleSolvers
using Test

import GeometricProblems.LotkaVolterra4d
import GeometricProblems.LotkaVolterra4dLagrangian
import GeometricProblems.LotkaVolterra4dLagrangian: Δt, nt, reference_solution


SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

SimpleSolvers.set_config(:nls_atol_break, 1E3)
SimpleSolvers.set_config(:nls_rtol_break, 1E3)
SimpleSolvers.set_config(:nls_stol_break, 1E3)


@testset "$(rpad("Lotka-Volterra 4D",80))" begin
    ode  = LotkaVolterra4dLagrangian.lotka_volterra_4d_ode()
    iode = LotkaVolterra4dLagrangian.lotka_volterra_4d_iode()
    idae = LotkaVolterra4dLagrangian.lotka_volterra_4d_idae()

    ref_ode  = LotkaVolterra4d.lotka_volterra_4d_ode()
    ref_iode = LotkaVolterra4d.lotka_volterra_4d_iode()
    ref_idae = LotkaVolterra4d.lotka_volterra_4d_idae()

    ref_equs = get_functions(ref_ode)
    ref_invs = get_invariants(ref_ode)

    v1 = zero(ode.q₀[begin])
    v2 = zero(ref_ode.q₀[begin])
    ode.v(ode.t₀, ode.q₀[begin], v1)
    ref_equs[:v](ref_ode.t₀, ref_ode.q₀[begin], v2)
    @test v1 ≈ v2  atol=1E-14

    h1 = ode.invariants[:h](ode.t₀, ode.q₀[begin])
    h2 = ref_invs[:h](ref_ode.t₀, ref_ode.q₀[begin])
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(ode, TableauGauss(2), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_ode, TableauGauss(2), Δt)
    ref_sol = integrate(ref_ode, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[end]) < 1E-14


    ref_equs = get_functions(ref_iode)
    ref_invs = get_invariants(ref_iode)

    @test iode.p₀ == ref_iode.p₀

    ϑ1 = zero(iode.q₀[begin])
    ϑ2 = zero(ref_iode.q₀[begin])
    iode.ϑ(iode.t₀, iode.q₀[begin], v1, ϑ1)
    ref_equs[:ϑ](ref_iode.t₀, ref_iode.q₀[begin], v2, ϑ2)
    @test ϑ1 ≈ ϑ2  atol=1E-14

    f1 = zero(iode.q₀[begin])
    f2 = zero(ref_iode.q₀[begin])
    iode.f(iode.t₀, iode.q₀[begin], v1, f1)
    ref_equs[:f](ref_iode.t₀, ref_iode.q₀[begin], v2, f2)
    @test f1 ≈ f2  atol=1E-14

    g1 = zero(iode.q₀[begin])
    g2 = zero(ref_iode.q₀[begin])
    iode.g(iode.t₀, iode.q₀[begin], v1, g1)
    ref_equs[:g](ref_iode.t₀, ref_iode.q₀[begin], v2, g2)
    @test g1 ≈ g2  atol=1E-14

    v1 = zero(iode.q₀[begin])
    v2 = zero(ref_iode.q₀[begin])
    iode.v̄(iode.t₀, iode.q₀[begin], v1)
    ref_equs[:v̄](ref_iode.t₀, ref_iode.q₀[begin], v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(iode.q₀[begin])
    f2 = zero(ref_iode.q₀[begin])
    iode.f̄(iode.t₀, iode.q₀[begin], v1, f1)
    ref_equs[:f̄](ref_iode.t₀, ref_iode.q₀[begin], v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = iode.invariants[:h](iode.t₀, iode.q₀[begin], zero(iode.q₀[begin]))
    h2 = ref_invs[:h](ref_iode.t₀, ref_iode.q₀[begin], zero(ref_iode.q₀[begin]))
    @test h1 ≈ h2  atol=1E-14

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 2E-3

    ref_int = IntegratorVPRKpMidpoint(ref_iode, TableauVPGLRK(2), Δt)
    ref_sol = integrate(ref_iode, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[end]) < 1E-14

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    ref_int = IntegratorVPRKpSymmetric(ref_iode, TableauVPGLRK(2), Δt)
    ref_sol = integrate(ref_iode, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[end]) < 1E-14


    ref_equs = get_functions(ref_idae)
    ref_invs = get_invariants(ref_idae)

    @test idae.p₀ == ref_idae.p₀

    v1 = zero(idae.q₀[begin])
    v2 = zero(ref_idae.q₀[begin])
    idae.v̄(idae.t₀, idae.q₀[begin], v1)
    ref_equs[:v̄](ref_idae.t₀, ref_idae.q₀[begin], v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(idae.q₀[begin])
    f2 = zero(ref_idae.q₀[begin])
    idae.f̄(idae.t₀, idae.q₀[begin], v1, f1)
    ref_equs[:f̄](ref_idae.t₀, ref_idae.q₀[begin], v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = idae.invariants[:h](idae.t₀, idae.q₀[begin], zero(idae.q₀[begin]))
    h2 = ref_invs[:h](ref_idae.t₀, ref_idae.q₀[begin], zero(ref_idae.q₀[begin]))
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 2E-3

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpMidpoint(2), Δt)
    ref_sol = integrate(ref_idae, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[end]) < 2E-14

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpSymmetric(2), Δt)
    ref_sol = integrate(ref_idae, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[end]) < 6E-14

end
