using Test
using GeometricIntegrators
using GeometricSolutions
using GeometricProblems.PerturbedPendulum

@testset "$(rpad("Perturbed Pendulum",80))" begin

    lode = lodeproblem([0.5],[0.6],timestep = 0.5,timespan = (0.0,100.0))
    hode = hodeproblem([0.5],[0.6],timestep = 0.5,timespan = (0.0,100.0))

    hode_sol = integrate(hode, ImplicitMidpoint())
    lode_sol = integrate(lode, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-12

end
