using Test
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillator: reference_solution, reference_solution_q, reference_solution_p
using GeometricSolutions


@testset "$(rpad("Harmonic Oscillator",80))" begin

    @test_nowarn odeproblem()
    @test_nowarn hodeproblem()
    @test_nowarn iodeproblem()
    @test_nowarn lodeproblem()
    @test_nowarn podeproblem()
    @test_nowarn sodeproblem()

    @test_nowarn degenerate_iodeproblem()
    @test_nowarn degenerate_lodeproblem()

    @test_nowarn daeproblem()
    @test_nowarn hdaeproblem()
    @test_nowarn idaeproblem()
    @test_nowarn ldaeproblem()
    @test_nowarn pdaeproblem()


    ode = odeproblem()
    iode = degenerate_iodeproblem()
    pode = podeproblem()
    hode = hodeproblem()
    ref  = exact_solution(ode)
    
    
    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-4

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-4


    sol = exact_solution(ode)
    @test sol.q[end] == reference_solution

    sol = exact_solution(pode)
    @test sol.q[end] == [reference_solution_q]
    @test sol.p[end] == [reference_solution_p]

    sol = exact_solution(hode)
    @test sol.q[end] == [reference_solution_q]
    @test sol.p[end] == [reference_solution_p]

end
