using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.LorenzAttractor
using GeometricProblems.LorenzAttractor: reference_solution


@testset "$(rpad("Lorenz Attractor",80))" begin
    ode  = lorenz_attractor_ode()

    sol = integrate(ode, Gauss(1))
    @test relative_maximum_error(sol.q, reference_solution) < 4E-2

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-5

end
