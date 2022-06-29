using SimpleSolvers
using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2dGauge
using GeometricProblems.LotkaVolterra2dGauge: Δt, nt, reference_solution

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())


@testset "$(rpad("Lotka-Volterra 2D with symmetric Lagrangian with gauge terms",80))" begin
    ode  = lotka_volterra_2d_ode()
    iode = lotka_volterra_2d_iode()
    idae = lotka_volterra_2d_idae()

    int = Integrator(ode, TableauGauss(2))
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 4E-4

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2))
    sol = integrate(iode, int, nt)
    # println(relative_maximum_error(sol.q, reference_solution))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2))
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2))
    sol = integrate(idae, int, nt)
    # println(relative_maximum_error(sol.q, reference_solution))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2))
    sol = integrate(idae, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4
end
