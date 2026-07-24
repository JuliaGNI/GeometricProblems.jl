using Test
using GeometricIntegrators
using GeometricProblems.OuterSolarSystem
using GeometricSolutions

@testset "$(rpad("Outer Solar System",80))" begin

    hode = @test_nowarn hodeproblem(timespan=(0.0,10.0))
    lode = @test_nowarn lodeproblem(timespan=(0.0,10.0))

    hode_sol = integrate(hode, ImplicitMidpoint())
    lode_sol = integrate(lode, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-11

end
