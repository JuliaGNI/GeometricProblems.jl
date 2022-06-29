using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.PointVortices
using GeometricProblems.PointVortices: reference_solution, nt


@testset "$(rpad("Point Vortices",80))" begin
    ode  = point_vortices_ode()
    iode = point_vortices_iode()
    idae = point_vortices_idae()

    int = Integrator(ode, TableauGauss(1))
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 3E-2

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(1))
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 3E-2

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(1))
    sol = integrate(idae, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 3E-2
end
