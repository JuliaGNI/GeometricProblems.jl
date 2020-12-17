using SimpleSolvers
using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra4d
using GeometricProblems.LotkaVolterra4d: Δt, nt, reference_solution

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())


@testset "$(rpad("Lotka-Volterra 4D",80))" begin
    ode  = lotka_volterra_4d_ode()
    iode = lotka_volterra_4d_iode()
    idae = lotka_volterra_4d_idae()

    int = Integrator(ode, TableauGLRK(2), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 2E-3

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 2E-3

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 8E-4
end
