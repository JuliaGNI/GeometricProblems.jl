using Test
using GeometricIntegrators
using GeometricProblems.DuffingOscillator
using GeometricSolutions


@testset "$(rpad("Duffing Oscillator",80))" begin

    @test_nowarn hodeproblem()
    @test_nowarn lodeproblem()

    hode = hodeproblem([1.2], [1.2]; timespan = (0, 100), timestep = 0.1)
    lode = lodeproblem()

    hode_sol = integrate(hode, Gauss(2))
    lode_sol = integrate(lode, Gauss(2))

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-9
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-9

end
