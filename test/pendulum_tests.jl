using Test
using GeometricIntegrators
using GeometricProblems.Pendulum
using GeometricSolutions


@testset "$(rpad("Pendulum",80))" begin

    @test_nowarn odeproblem()
    @test_nowarn podeproblem()
    @test_nowarn iodeproblem()
    @test_nowarn idaeproblem()


    ode  = odeproblem()
    pode = podeproblem()
    iode = iodeproblem()

    ode_sol  = integrate(ode,  Gauss(2))
    pode_sol = integrate(pode, Gauss(2))
    iode_sol = integrate(iode, MidpointProjection(VPRKGauss(2)))

    # ODE state is [θ, θ̇]; PODE has q=[θ], p=[θ̇]: compare final values
    @test ode_sol.q[end][1] ≈ pode_sol.q[end][1]
    @test ode_sol.q[end][2] ≈ pode_sol.p[end][1]

    # IODE degenerate state is the same [θ, θ̇] as ODE
    @test relative_maximum_error(ode_sol.q, iode_sol.q) < 1E-4

    iode_sol2 = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(iode_sol.q, iode_sol2.q) < 1E-4

end
