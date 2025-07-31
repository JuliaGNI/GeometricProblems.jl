using Test
using GeometricIntegrators
using GeometricSolutions
using GeometricProblems.OuterSolarSystem

@testset "$(rpad("Outer Solar System",80))" begin

    lode = OuterSolarSystem.lodeproblem(timespan = (0,3000),n=2)
    hode = OuterSolarSystem.hodeproblem(timespan = (0,3000),n=2)

    hode_sol = integrate(hode, Gauss(2))
    lode_sol = integrate(lode, Gauss(2))

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-6
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-6

end