using SimpleSolvers
using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.LorenzAttractor

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())

const Δt = 0.01
const nt = 1000

const ref = [-4.902687541134471, -3.743872921802973, 24.690858102790042]


@testset "$(rpad("Lorenz Attractor",80))" begin
    ode  = lorenz_attractor_ode()

    sol = integrate(ode, TableauGauss(1), Δt, nt)
    @test relative_maximum_error(sol.q, ref) < 4E-2

    sol = integrate(ode, TableauGauss(2), Δt, nt)
    @test relative_maximum_error(sol.q, ref) < 2E-5

end
