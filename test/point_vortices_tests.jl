using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.PointVortices


const Δt = 0.01
const nt = 1000

const ref = [0.18722529318641928, 0.38967432450068706, 0.38125332930294187, 0.4258020604293123]


@testset "$(rpad("Point Vortices",80))" begin
    ode  = point_vortices_ode()
    iode = point_vortices_iode()
    idae = point_vortices_idae()

    int = Integrator(ode, TableauGauss(1), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, ref) < 3E-2

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, ref) < 3E-2

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, ref) < 3E-2
end
