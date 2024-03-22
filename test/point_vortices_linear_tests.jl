using Test
using GeometricIntegrators
using GeometricProblems.PointVorticesLinear
using GeometricSolutions


@testset "$(rpad("Point Vortices (linear)",80))" begin
    ode  = odeproblem()
    iode = iodeproblem()
    ref  = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-9

    sol = integrate(iode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-9
end
