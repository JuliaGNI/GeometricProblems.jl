
using Test
using GeometricIntegrators
using GeometricIntegrators.Utils
using GeometricProblems.HarmonicOscillator

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())
set_config(:nls_stol_break, 1E3)

const Δt = 0.1
const nt = 1000

const ref = [-0.012420428712136283, -0.35344429048339926]

ode  = harmonic_oscillator_ode()
dae  = harmonic_oscillator_dae()
iode = harmonic_oscillator_iode()
idae = harmonic_oscillator_idae()
pode = harmonic_oscillator_pode()
pdae = harmonic_oscillator_pdae()

@testset "$(rpad("Harmonic Oscillator",80))" begin

    int = Integrator(ode, getTableauGLRK(2), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpMidpoint(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(2), Δt)
    sol = integrate(int, nt)
    @test rel_err(sol.q, ref) < 1E-4

end