
using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2d

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

const Δt = 0.01
const nt = 1000

const ref = [0.530592013081557, 1.1995663801610505]


@testset "$(rpad("Lotka-Volterra 2D",80))" begin
    ode  = lotka_volterra_2d_ode()
    iode = lotka_volterra_2d_iode()
    idae = lotka_volterra_2d_idae()

    int = Integrator(ode, getTableauGLRK(1), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = Integrator(idae, getTableauVSPARKGLRKpMidpoint(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = Integrator(idae, getTableauVSPARKGLRKpSymmetric(1), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, ref) < 1E-4
end
