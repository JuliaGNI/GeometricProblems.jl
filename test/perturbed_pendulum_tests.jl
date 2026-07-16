using Test
using GeometricIntegrators
using GeometricProblems.PerturbedPendulum
using GeometricSolutions


@testset "$(rpad("Perturbed Pendulum",80))" begin

    hode = @test_nowarn hodeproblem()
    lode = @test_nowarn lodeproblem()

    hode_sol = integrate(hode, ImplicitMidpoint())
    lode_sol = integrate(lode, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1e-12
    @test isapprox(collect(lode_sol.p[:,1]), collect(hode_sol.p[:,1]); atol=1e-12) 
end