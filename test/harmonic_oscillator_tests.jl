using SimpleSolvers
using Test
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
    @test_nowarn ode_equs[:v](tbegin(ode), ode.ics.q, zero(ode.ics.q))

    iode_equs = functions(iode)
    @test_nowarn iode_equs[:ϑ](tbegin(iode), iode.ics.q, iode.ics.p, zero(iode.ics.q))
    @test_nowarn iode_equs[:f](tbegin(iode), iode.ics.q, iode.ics.p, zero(iode.ics.p))
    @test_nowarn iode_equs[:v̄](tbegin(iode), iode.ics.q, zero(iode.ics.q))
    @test_nowarn iode_equs[:f̄](tbegin(iode), iode.ics.q, iode.ics.p, zero(iode.ics.q))

    pode_equs = functions(pode)
    @test_nowarn pode_equs[:v](tbegin(pode), pode.ics.q, pode.ics.p, zero(pode.ics.q))
    @test_nowarn pode_equs[:f](tbegin(pode), pode.ics.q, pode.ics.p, zero(pode.ics.p))

    dae_equs = functions(dae)
    @test_nowarn dae_equs[:v](tbegin(dae), dae.ics.q, zero(dae.ics.q))
    @test_nowarn dae_equs[:u](tbegin(dae), dae.ics.q, dae.ics.λ, zero(dae.ics.q))
    @test_nowarn dae_equs[:ϕ](tbegin(dae), dae.ics.q, zero(dae.ics.λ))
    @test_nowarn dae_equs[:v̄](tbegin(dae), dae.ics.q, zero(dae.ics.q))

    idae_equs = functions(idae)
    @test_nowarn idae_equs[:ϑ](tbegin(idae), idae.ics.q, idae.ics.λ, zero(idae.ics.q))
    @test_nowarn idae_equs[:f](tbegin(idae), idae.ics.q, idae.ics.λ, zero(idae.ics.p))
    @test_nowarn idae_equs[:u](tbegin(idae), idae.ics.q, pdae.ics.p, idae.ics.λ, zero(idae.ics.q))
    @test_nowarn idae_equs[:g](tbegin(idae), idae.ics.q, pdae.ics.p, idae.ics.λ, zero(idae.ics.p))
    @test_nowarn idae_equs[:ϕ](tbegin(idae), idae.ics.q, idae.ics.p, zero(idae.ics.λ))
    @test_nowarn idae_equs[:v̄](tbegin(idae), idae.ics.q, zero(idae.ics.q))
    @test_nowarn idae_equs[:f̄](tbegin(idae), idae.ics.q, idae.ics.p, zero(idae.ics.q))

    pdae_equs = functions(pdae)
    @test_nowarn pdae_equs[:v](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.q))
    @test_nowarn pdae_equs[:f](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.p))
    @test_nowarn pdae_equs[:u](tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ, zero(pdae.ics.q))
    @test_nowarn pdae_equs[:g](tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ, zero(pdae.ics.p))
    @test_nowarn pdae_equs[:ϕ](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.λ))
    @test_nowarn pdae_equs[:v̄](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.q))
    @test_nowarn pdae_equs[:f̄](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.p))

    
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
