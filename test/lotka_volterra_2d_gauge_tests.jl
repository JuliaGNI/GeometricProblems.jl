using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2dGauge
using GeometricProblems.LotkaVolterra2dGauge: reference_solution


@testset "$(rpad("Lotka-Volterra 2D with symmetric Lagrangian with gauge terms",80))" begin
    ode  = lotka_volterra_2d_ode()
    iode = lotka_volterra_2d_iode()
    idae = lotka_volterra_2d_idae()

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    # println(relative_maximum_error(sol.q, reference_solution))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    # println(relative_maximum_error(sol.q, reference_solution))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4
end
