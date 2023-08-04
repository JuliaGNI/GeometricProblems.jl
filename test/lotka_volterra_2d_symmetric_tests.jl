using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2dSymmetric
using GeometricSolutions


@testset "$(rpad("Lotka-Volterra 2D with symmetric Lagrangian",80))" begin
    ode  = lotka_volterra_2d_ode()
    iode = lotka_volterra_2d_iode()
    idae = lotka_volterra_2d_idae()
    ref  = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    # println(relative_maximum_error(sol.q, ref.q))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4
end
