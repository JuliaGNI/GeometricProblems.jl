using Test
using GeometricIntegrators
using GeometricProblems.LorenzAttractor
using GeometricProblems.LorenzAttractor: reference_solution
using GeometricSolutions


@testset "$(rpad("Lorenz Attractor",80))" begin
    ode  = lorenz_attractor_ode()

    sol = integrate(ode, Gauss(1))
    @test relative_maximum_error(sol.q[end], reference_solution) < 4E-2

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q[end], reference_solution) < 2E-5

end
