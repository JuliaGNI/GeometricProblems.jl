using Test
using GeometricIntegrators
using GeometricProblems.DoublePendulum
using GeometricSolutions


@testset "$(rpad("Double Pendulum",80))" begin

    @test_nowarn hodeproblem()
    @test_nowarn lodeproblem()

    hode = hodeproblem()
    lode = lodeproblem()

    hode_sol = integrate(hode, Gauss(2))
    lode_sol = integrate(lode, Gauss(2))

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-13
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-11

end
