using SimpleSolvers
using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.LorenzAttractor
using GeometricProblems.LorenzAttractor: reference_solution, nt

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())


@testset "$(rpad("Lorenz Attractor",80))" begin
    ode  = lorenz_attractor_ode()

    sol = integrate(ode, TableauGauss(1), nt)
    @test relative_maximum_error(sol.q, reference_solution) < 4E-2

    sol = integrate(ode, TableauGauss(2), nt)
    @test relative_maximum_error(sol.q, reference_solution) < 2E-5

end
