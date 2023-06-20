using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra4d
using GeometricProblems.LotkaVolterra4d: reference_solution


@testset "$(rpad("Lotka-Volterra 4D",80))" begin
    ode  = lotka_volterra_4d_ode()
    iode = lotka_volterra_4d_iode()
    idae = lotka_volterra_4d_idae()

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-3

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, reference_solution) < 8E-4
end
