using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.PointVortices
using GeometricSolutions


@testset "$(rpad("Point Vortices",80))" begin
    ode  = point_vortices_ode()
    iode = point_vortices_iode()
    idae = point_vortices_idae()
    ref  = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-3

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-3

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-3
end
