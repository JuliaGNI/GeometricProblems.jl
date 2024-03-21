using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d
using GeometricSolutions


@testset "$(rpad("Lotka-Volterra 2d",80))" begin

    ode  = odeproblem()
    hode = hodeproblem()
    iode = iodeproblem()
    pode = podeproblem()
    lode = lodeproblem()
    dae  = daeproblem()
    hdae = hdaeproblem()
    idae = idaeproblem()
    pdae = pdaeproblem()
    ldae = ldaeproblem()

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
