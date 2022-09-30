using SimpleSolvers
using Test
using GeometricBase
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillator: reference_solution

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())
SimpleSolvers.set_config(:nls_stol_break, Inf)


@testset "$(rpad("Harmonic Oscillator",80))" begin

    ode  = harmonic_oscillator_ode()
    pode = harmonic_oscillator_pode()
    iode = harmonic_oscillator_iode()
    dae  = harmonic_oscillator_dae()
    idae = harmonic_oscillator_idae()
    pdae = harmonic_oscillator_pdae()

    ode_equs = functions(ode)
    @test_nowarn ode_equs[:v](zero(ode.ics.q), tbegin(ode), ode.ics.q)

    iode_equs = functions(iode)
    @test_nowarn iode_equs[:ϑ](zero(iode.ics.p), tbegin(iode), iode.ics.q, iode.ics.p)
    @test_nowarn iode_equs[:f](zero(iode.ics.p), tbegin(iode), iode.ics.q, iode.ics.p)
    @test_nowarn iode_equs[:v̄](zero(iode.ics.q), tbegin(iode), iode.ics.q)
    @test_nowarn iode_equs[:f̄](zero(iode.ics.p), tbegin(iode), iode.ics.q, iode.ics.p)

    pode_equs = functions(pode)
    @test_nowarn pode_equs[:v](zero(pode.ics.q), tbegin(pode), pode.ics.q, pode.ics.p)
    @test_nowarn pode_equs[:f](zero(pode.ics.p), tbegin(pode), pode.ics.q, pode.ics.p)

    dae_equs = functions(dae)
    @test_nowarn dae_equs[:v](zero(dae.ics.q), tbegin(dae), dae.ics.q)
    @test_nowarn dae_equs[:u](zero(dae.ics.q), tbegin(dae), dae.ics.q, dae.ics.λ)
    @test_nowarn dae_equs[:ϕ](zero(dae.ics.λ), tbegin(dae), dae.ics.q)
    @test_nowarn dae_equs[:v̄](zero(dae.ics.q), tbegin(dae), dae.ics.q)

    idae_equs = functions(idae)
    @test_nowarn idae_equs[:ϑ](zero(idae.ics.q), tbegin(idae), idae.ics.q, idae.ics.λ)
    @test_nowarn idae_equs[:f](zero(idae.ics.p), tbegin(idae), idae.ics.q, idae.ics.λ)
    @test_nowarn idae_equs[:u](zero(idae.ics.q), tbegin(idae), idae.ics.q, pdae.ics.p, idae.ics.λ)
    @test_nowarn idae_equs[:g](zero(idae.ics.p), tbegin(idae), idae.ics.q, pdae.ics.p, idae.ics.λ)
    @test_nowarn idae_equs[:ϕ](zero(idae.ics.λ), tbegin(idae), idae.ics.q, idae.ics.p)
    @test_nowarn idae_equs[:v̄](zero(idae.ics.q), tbegin(idae), idae.ics.q)
    @test_nowarn idae_equs[:f̄](zero(idae.ics.q), tbegin(idae), idae.ics.q, idae.ics.p)

    pdae_equs = functions(pdae)
    @test_nowarn pdae_equs[:v](zero(pdae.ics.q), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    @test_nowarn pdae_equs[:f](zero(pdae.ics.p), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    @test_nowarn pdae_equs[:u](zero(pdae.ics.q), tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ)
    @test_nowarn pdae_equs[:g](zero(pdae.ics.p), tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ)
    @test_nowarn pdae_equs[:ϕ](zero(pdae.ics.λ), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    @test_nowarn pdae_equs[:v̄](zero(pdae.ics.q), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    @test_nowarn pdae_equs[:f̄](zero(pdae.ics.p), tbegin(pdae), pdae.ics.q, pdae.ics.p)

    
    int = Integrator(ode, TableauGauss(2))
    sol = integrate(ode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 1E-4

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2))
    sol = integrate(iode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 1E-4

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2))
    sol = integrate(iode, int)
    @test relative_maximum_error(sol.q, reference_solution) < 1E-4

end
