
using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
import GeometricProblems.LotkaVolterra4dLagrangian
using GeometricProblems.LotkaVolterra4dLagrangian: Δt, nt, reference_solution

using GeometricProblems.LotkaVolterra4d

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

set_config(:nls_atol_break, 1E3)
set_config(:nls_rtol_break, 1E3)
set_config(:nls_stol_break, 1E3)


@testset "$(rpad("Lotka-Volterra 4D",80))" begin
    ode  = LotkaVolterra4dLagrangian.lotka_volterra_4d_ode()
    iode = LotkaVolterra4dLagrangian.lotka_volterra_4d_iode()
    idae = LotkaVolterra4dLagrangian.lotka_volterra_4d_idae()

    ref_ode  = LotkaVolterra4d.lotka_volterra_4d_ode()
    ref_iode = LotkaVolterra4d.lotka_volterra_4d_iode()
    ref_idae = LotkaVolterra4d.lotka_volterra_4d_idae()

    ref_eqs = get_function_tuple(ref_ode)

    v1 = zero(ode.q₀)
    v2 = zero(ref_ode.q₀)
    ode.v(ode.t₀, ode.q₀, v1)
    ref_eqs[:v](ref_ode.t₀, ref_ode.q₀, v2)
    @test v1 ≈ v2  atol=1E-14

    h1 = ode.h(ode.t₀, ode.q₀)
    h2 = ref_eqs[:h](ref_ode.t₀, ref_ode.q₀)
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(ode, getTableauGLRK(2), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_ode, getTableauGLRK(2), Δt)
    ref_sol = integrate(ref_ode, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[:,end]) < 1E-14


    ref_eqs = get_function_tuple(ref_iode)

    @test iode.p₀ == ref_iode.p₀

    v1 = zero(iode.q₀)
    v2 = zero(ref_iode.q₀)
    iode.v̄(iode.t₀, iode.q₀, v1)
    ref_eqs[:v̄](ref_iode.t₀, ref_iode.q₀, v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(iode.q₀)
    f2 = zero(ref_iode.q₀)
    iode.f̄(iode.t₀, iode.q₀, v1, f1)
    ref_eqs[:f̄](ref_iode.t₀, ref_iode.q₀, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = iode.h(iode.t₀, iode.q₀)
    h2 = ref_eqs[:h](ref_iode.t₀, ref_iode.q₀)
    @test h1 ≈ h2  atol=1E-14

    int = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 2E-3

    ref_int = IntegratorVPRKpMidpoint(ref_iode, getTableauVPGLRK(2), Δt)
    ref_sol = integrate(ref_iode, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[:,end]) < 1E-14

    int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    ref_int = IntegratorVPRKpSymmetric(ref_iode, getTableauVPGLRK(2), Δt)
    ref_sol = integrate(ref_iode, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[:,end]) < 1E-14


    ref_eqs = get_function_tuple(ref_idae)

    @test idae.p₀ == ref_idae.p₀

    v1 = zero(idae.q₀)
    v2 = zero(ref_idae.q₀)
    idae.v̄(idae.t₀, idae.q₀, v1)
    ref_eqs[:v̄](ref_idae.t₀, ref_idae.q₀, v2)
    @test v1 ≈ v2  atol=1E-14

    f1 = zero(idae.q₀)
    f2 = zero(ref_idae.q₀)
    idae.f̄(idae.t₀, idae.q₀, v1, f1)
    ref_eqs[:f̄](ref_idae.t₀, ref_idae.q₀, v2, f2)
    @test f1 ≈ f2  atol=1E-14

    h1 = idae.h(idae.t₀, idae.q₀)
    h2 = ref_eqs[:h](ref_idae.t₀, ref_idae.q₀)
    @test h1 ≈ h2  atol=1E-14

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 2E-3

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpMidpoint(2), Δt)
    ref_sol = integrate(ref_idae, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[:,end]) < 1E-14

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    ref_int = Integrator(ref_idae, TableauVSPARKGLRKpSymmetric(2), Δt)
    ref_sol = integrate(ref_idae, ref_int, nt)
    @test rel_err(sol.q, ref_sol.q[:,end]) < 4E-14

end
