
using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra3d

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())

const Δt = 0.01
const nt = 1000

const ref = [0.39947308320241187, 1.9479527336244262, 2.570183075433086]


@testset "$(rpad("Lotka-Volterra 3D",80))" begin
    ode  = lotka_volterra_3d_ode()

    int = Integrator(ode, TableauGLRK(1), Δt)
    sol = integrate(ode, int, nt)
    H, ΔH = compute_energy_error(sol.t, sol.q)
    C, ΔC = compute_casimir_error(sol.t, sol.q)
    @test rel_err(sol.q, ref) < 5E-4
    @test ΔH[end] < 4E-6
    @test ΔC[end] < 8E-6

    int = Integrator(ode, TableauGLRK(2), Δt)
    sol = integrate(ode, int, nt)
    H, ΔH = compute_energy_error(sol.t, sol.q)
    C, ΔC = compute_casimir_error(sol.t, sol.q)
    @test rel_err(sol.q, ref) < 2E-9
    @test ΔH[end] < 5E-11
    @test ΔC[end] < 2E-11

end
