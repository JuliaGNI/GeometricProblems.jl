using SimpleSolvers
using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.HarmonicOscillator

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())
SimpleSolvers.set_config(:nls_stol_break, 1E3)

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

    ode_equs = get_function_tuple(ode)
    @test_nowarn ode_equs[:v](ode.t₀, ode.q₀[begin], zero(ode.q₀[begin]))

    iode_equs = get_function_tuple(iode)
    @test_nowarn iode_equs[:ϑ](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.p₀[begin]))
    @test_nowarn iode_equs[:v̄](iode.t₀, iode.q₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f̄](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.q₀[begin]))

    pode_equs = get_function_tuple(pode)
    @test_nowarn pode_equs[:v](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.q₀[begin]))
    @test_nowarn pode_equs[:f](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.p₀[begin]))

    dae_equs = get_function_tuple(dae)
    @test_nowarn dae_equs[:v](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:u](dae.t₀, dae.q₀[begin], dae.λ₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:ϕ](dae.t₀, dae.q₀[begin], zero(dae.λ₀[begin]))
    @test_nowarn dae_equs[:v̄](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))

    idae_equs = get_function_tuple(idae)
    @test_nowarn idae_equs[:ϑ](idae.t₀, idae.q₀[begin], idae.λ₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f](idae.t₀, idae.q₀[begin], idae.λ₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:u](idae.t₀, idae.q₀[begin], pdae.p₀[begin], idae.λ₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:g](idae.t₀, idae.q₀[begin], pdae.p₀[begin], idae.λ₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:ϕ](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.λ₀[begin]))
    @test_nowarn idae_equs[:v̄](idae.t₀, idae.q₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f̄](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.q₀[begin]))

    pdae_equs = get_function_tuple(pdae)
    @test_nowarn pdae_equs[:v](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:u](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:g](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:ϕ](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.λ₀[begin]))
    @test_nowarn pdae_equs[:v̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))

    
    int = Integrator(ode, TableauGauss(2), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, ref) < 1E-4

end
