using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.PointVortices
using GeometricSolutions


@testset "$(rpad("Point Vortices",80))" begin
    ode  = odeproblem()
    iode = iodeproblem()
    idae = idaeproblem()
    ref  = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-5

    sol = integrate(iode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 2E-5

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-5

    # the formal (phase-space) Lagrangian LODE reproduces the ODE dynamics
    lode = lodeproblem_formal_lagrangian()
    sol = integrate(lode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-4
end
