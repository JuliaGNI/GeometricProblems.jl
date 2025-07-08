using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra4d
using GeometricSolutions


@testset "$(rpad("Lotka-Volterra 4D",80))" begin

    ode = odeproblem()
    iode = iodeproblem()
    lode = lodeproblem()
    pode = podeproblem()
    idae = idaeproblem()
    ldae = ldaeproblem()
    pdae = pdaeproblem()

    ldae_secondary = ldaeproblem_secondary()

    ref = integrate(ode, Gauss(8))


    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(iode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(lode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(pode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(idae, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(ldae, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(pdae, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8


    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(lode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(lode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(ldae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(ldae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8


    sol = integrate(ldae_secondary, TableauVSPARKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-5

    # TODO: Reactivate!
    # sol = integrate(pdae, TableauHSPARKLobattoIIIE(2))
    # @test relative_maximum_error(sol.q, ref.q) < 4E-4

end
