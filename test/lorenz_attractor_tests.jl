
using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.LorenzAttractor

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

const Δt = 0.01
const nt = 1000

const ref = [-4.902687541134471, -3.743872921802973, 24.690858102790042]


@testset "$(rpad("Lorenz Attractor",80))" begin
    ode  = lorenz_attractor_ode()

    int = Integrator(ode, getTableauGLRK(1), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 4E-2

    int = Integrator(ode, getTableauGLRK(2), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 2E-5

end
