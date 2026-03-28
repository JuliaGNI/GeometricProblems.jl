using Test
using GeometricIntegrators
using GeometricProblems.LotkaVolterra3d
using GeometricSolutions


@testset "$(rpad("Lotka-Volterra 3D",80))" begin
    ode = odeproblem()
    ref = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(1))
    H, ΔH = compute_energy_error(sol.t, sol.q, parameters(ode))
    C, ΔC = compute_casimir_error(sol.t, sol.q, parameters(ode))
    @test relative_maximum_error(sol.q, ref.q) < 5E-4
    @test ΔH[end] < 4E-6
    @test ΔC[end] < 8E-6

    sol = integrate(ode, Gauss(2))
    H, ΔH = compute_energy_error(sol.t, sol.q, parameters(ode))
    C, ΔC = compute_casimir_error(sol.t, sol.q, parameters(ode))
    @test relative_maximum_error(sol.q, ref.q) < 2E-9
    @test ΔH[end] < 5E-11
    @test ΔC[end] < 2E-11

end
