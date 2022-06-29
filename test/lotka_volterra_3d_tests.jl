
using SimpleSolvers
using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra3d
using GeometricProblems.LotkaVolterra3d: reference_solution, nt

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())


@testset "$(rpad("Lotka-Volterra 3D",80))" begin
    ode  = lotka_volterra_3d_ode()

    int = Integrator(ode, TableauGauss(1))
    sol = integrate(ode, int, nt)
    H, ΔH = compute_energy_error(sol.t, sol.q, parameters(ode))
    C, ΔC = compute_casimir_error(sol.t, sol.q, parameters(ode))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4
    @test ΔH[end] < 4E-6
    @test ΔC[end] < 8E-6

    int = Integrator(ode, TableauGauss(2))
    sol = integrate(ode, int, nt)
    H, ΔH = compute_energy_error(sol.t, sol.q, parameters(ode))
    C, ΔC = compute_casimir_error(sol.t, sol.q, parameters(ode))
    @test relative_maximum_error(sol.q, reference_solution) < 2E-9
    @test ΔH[end] < 5E-11
    @test ΔC[end] < 2E-11

end
