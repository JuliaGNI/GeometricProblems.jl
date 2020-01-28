
using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2dsingular

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

const Δt = 0.01
const nt = 1000

const ref = [0.35005370250145684, 2.115866819186877]


@testset "$(rpad("Lotka-Volterra 2D with singular Lagrangian",80))" begin
    ode  = lotka_volterra_2d_ode()
    iode = lotka_volterra_2d_iode()
    idae = lotka_volterra_2d_idae()

    int = Integrator(ode, getTableauGLRK(2), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 4E-11

    int = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 5E-5

    int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 5E-5

    int = Integrator(idae, getTableauVSPARKGLRKpMidpoint(1), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 5E-5

    int = Integrator(idae, getTableauVSPARKGLRKpSymmetric(1), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 5E-5
end
