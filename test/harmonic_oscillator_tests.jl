
using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.HarmonicOscillator

set_config(:nls_atol, 8eps())
set_config(:nls_rtol, 2eps())
set_config(:nls_stol_break, 1E3)

const Δt = 0.1
const nt = 1000

const ref = [-0.012420428712136283, -0.35344429048339926]

@testset "$(rpad("Harmonic Oscillator",80))" begin

    ode  = harmonic_oscillator_ode()
    pode = harmonic_oscillator_pode()
    iode = harmonic_oscillator_iode()
    dae  = harmonic_oscillator_dae()
    idae = harmonic_oscillator_idae()
    pdae = harmonic_oscillator_pdae()

    int = Integrator(ode, TableauGLRK(2), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

end
