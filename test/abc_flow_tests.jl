using Test
using GeometricIntegrators
using GeometricProblems.ABCFlow
using GeometricSolutions


@testset "$(rpad("ABC Flow",80))" begin

    @test_nowarn odeproblem()

    ode = odeproblem([0.5, 0., 0.])

    ref_sol = integrate(ode, Gauss(8))
    ode_sol = integrate(ode, Gauss(2))

    @test relative_maximum_error(ode_sol.q, ref_sol.q) < 1E-5

end
