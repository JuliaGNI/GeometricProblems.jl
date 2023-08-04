using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d
using GeometricSolutions


@testset "$(rpad("Lotka-Volterra 2d",80))" begin

    ode  = lotka_volterra_2d_ode()
    hode = lotka_volterra_2d_hode()
    iode = lotka_volterra_2d_iode()
    pode = lotka_volterra_2d_pode()
    lode = lotka_volterra_2d_lode()
    dae  = lotka_volterra_2d_dae()
    hdae = lotka_volterra_2d_hdae()
    idae = lotka_volterra_2d_idae()
    pdae = lotka_volterra_2d_pdae()
    ldae = lotka_volterra_2d_ldae()

    ref  = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4
end
