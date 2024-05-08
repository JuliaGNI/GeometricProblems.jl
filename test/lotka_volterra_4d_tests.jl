using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra4d
using GeometricSolutions


@testset "$(rpad("Lotka-Volterra 4D",80))" begin
    ode  = odeproblem()
    iode = iodeproblem()
    idae = idaeproblem()
    ref  = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-3

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4
end
