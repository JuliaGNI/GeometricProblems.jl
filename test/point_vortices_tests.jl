using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricIntegrators.Utils
using GeometricProblems.PointVortices
using GeometricProblems.PointVortices: reference_solution


@testset "$(rpad("Point Vortices",80))" begin
    ode  = point_vortices_ode()
    iode = point_vortices_iode()
    idae = point_vortices_idae()

    sol = integrate(ode, Gauss(1))
    @test relative_maximum_error(sol.q, reference_solution) < 3E-2

    sol = integrate(iode, SymmetricProjection(VPRKGauss(1)))
    @test relative_maximum_error(sol.q, reference_solution) < 3E-2

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(1))
    @test relative_maximum_error(sol.q, reference_solution) < 3E-2
end
