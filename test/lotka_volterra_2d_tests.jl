using SimpleSolvers
using Test
using GeometricEquations
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2d
using GeometricProblems.LotkaVolterra2d: nt, reference_solution

SimpleSolvers.set_config(:nls_atol, 8eps())
SimpleSolvers.set_config(:nls_rtol, 2eps())


@testset "$(rpad("Lotka-Volterra 2d",80))" begin

    ode  = lotka_volterra_2d_ode()
    hode = lotka_volterra_2d_hode()
    iode = lotka_volterra_2d_iode()
    pode = lotka_volterra_2d_pode()
    lode = lotka_volterra_2d_lode()
    dae  = lotka_volterra_2d_dae()
    hdae = lotka_volterra_2d_hdae()
    idae = lotka_volterra_2d_idae()
    pdae = lotka_volterra_2d_pdae()
    ldae = lotka_volterra_2d_ldae()

    ode_equs = functions(ode)
    ode_invs = invariants(ode)
    @test_nowarn ode_equs[:v](tbegin(ode), ode.ics.q, zero(ode.ics.q))
    @test_nowarn ode_invs[:h](tbegin(ode), ode.ics.q)

    hode_equs = functions(hode)
    @test_nowarn hode_equs[:v](tbegin(hode), hode.ics.q, hode.ics.p, zero(hode.ics.q))
    @test_nowarn hode_equs[:f](tbegin(hode), hode.ics.q, hode.ics.p, zero(hode.ics.p))
    @test_nowarn hode_equs[:h](tbegin(hode), hode.ics.q, hode.ics.p)

    iode_equs = functions(iode)
    iode_invs = invariants(iode)
    @test_nowarn iode_equs[:ϑ](tbegin(iode), iode.ics.q, iode.ics.p, zero(iode.ics.q))
    @test_nowarn iode_equs[:f](tbegin(iode), iode.ics.q, iode.ics.p, zero(iode.ics.p))
    @test_nowarn iode_equs[:v̄](tbegin(iode), iode.ics.q, zero(iode.ics.q))
    @test_nowarn iode_equs[:f̄](tbegin(iode), iode.ics.q, iode.ics.p, zero(iode.ics.p))
    @test_nowarn iode_invs[:h](tbegin(iode), iode.ics.q, zero(iode.ics.q))

    pode_equs = functions(pode)
    pode_invs = invariants(pode)
    @test_nowarn pode_equs[:v](tbegin(pode), pode.ics.q, pode.ics.p, zero(pode.ics.q))
    @test_nowarn pode_equs[:f](tbegin(pode), pode.ics.q, pode.ics.p, zero(pode.ics.p))
    @test_nowarn pode_invs[:h](tbegin(pode), pode.ics.q, pode.ics.p)

    lode_equs = functions(lode)
    lode_invs = invariants(lode)
    @test_nowarn lode_equs[:ϑ](tbegin(lode), lode.ics.q, lode.ics.p, zero(lode.ics.q))
    @test_nowarn lode_equs[:f](tbegin(lode), lode.ics.q, lode.ics.p, zero(lode.ics.p))
    @test_nowarn lode_equs[:v̄](tbegin(lode), lode.ics.q, zero(lode.ics.q))
    @test_nowarn lode_equs[:f̄](tbegin(lode), lode.ics.q, lode.ics.p, zero(lode.ics.p))
    @test_nowarn lode_invs[:h](tbegin(lode), lode.ics.q, zero(iode.ics.q))

    dae_equs = functions(dae)
    dae_invs = invariants(dae)
    @test_nowarn dae_equs[:v](tbegin(dae), dae.ics.q, zero(dae.ics.q))
    @test_nowarn dae_equs[:u](tbegin(dae), dae.ics.q, dae.ics.λ, zero(dae.ics.q))
    @test_nowarn dae_equs[:ϕ](tbegin(dae), dae.ics.q, zero(dae.ics.λ))
    @test_nowarn dae_equs[:v̄](tbegin(dae), dae.ics.q, zero(dae.ics.q))
    @test_nowarn dae_invs[:h](tbegin(dae), dae.ics.q)

    hdae_equs = functions(hdae)
    @test_nowarn hdae_equs[:v](tbegin(hdae), hdae.ics.q, hdae.ics.p, zero(hdae.ics.q))
    @test_nowarn hdae_equs[:f](tbegin(hdae), hdae.ics.q, hdae.ics.p, zero(hdae.ics.p))
    @test_nowarn hdae_equs[:u](tbegin(hdae), hdae.ics.q, hdae.ics.p, hdae.ics.λ, zero(hdae.ics.q))
    @test_nowarn hdae_equs[:g](tbegin(hdae), hdae.ics.q, hdae.ics.p, hdae.ics.λ, zero(hdae.ics.p))
    @test_nowarn hdae_equs[:ϕ](tbegin(hdae), hdae.ics.q, hdae.ics.p, zero(hdae.ics.λ))
    @test_nowarn hdae_equs[:v̄](tbegin(hdae), hdae.ics.q, hdae.ics.p, zero(hdae.ics.q))
    @test_nowarn hdae_equs[:f̄](tbegin(hdae), hdae.ics.q, hdae.ics.p, zero(hdae.ics.p))
    @test_nowarn hdae_equs[:h](tbegin(hdae), hdae.ics.q, hdae.ics.p)

    idae_equs = functions(idae)
    idae_invs = invariants(idae)
    @test_nowarn idae_equs[:ϑ](tbegin(idae), idae.ics.q, idae.ics.p, zero(idae.ics.q))
    @test_nowarn idae_equs[:f](tbegin(idae), idae.ics.q, idae.ics.p, zero(idae.ics.p))
    @test_nowarn idae_equs[:u](tbegin(idae), idae.ics.q, idae.ics.p, idae.ics.λ, zero(idae.ics.q))
    @test_nowarn idae_equs[:g](tbegin(idae), idae.ics.q, idae.ics.p, idae.ics.λ, zero(idae.ics.p))
    @test_nowarn idae_equs[:ϕ](tbegin(idae), idae.ics.q, idae.ics.p, zero(idae.ics.λ))
    @test_nowarn idae_equs[:v̄](tbegin(idae), idae.ics.q, zero(idae.ics.q))
    @test_nowarn idae_equs[:f̄](tbegin(idae), idae.ics.q, idae.ics.p, zero(idae.ics.p))
    @test_nowarn idae_invs[:h](tbegin(idae), idae.ics.q, zero(iode.ics.q))

    pdae_equs = functions(pdae)
    pdae_invs = invariants(pdae)
    @test_nowarn pdae_equs[:v](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.q))
    @test_nowarn pdae_equs[:f](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.p))
    @test_nowarn pdae_equs[:u](tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ, zero(pdae.ics.q))
    @test_nowarn pdae_equs[:g](tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ, zero(pdae.ics.p))
    @test_nowarn pdae_equs[:ϕ](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.λ))
    @test_nowarn pdae_equs[:v̄](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.q))
    @test_nowarn pdae_equs[:f̄](tbegin(pdae), pdae.ics.q, pdae.ics.p, zero(pdae.ics.p))
    @test_nowarn pdae_invs[:h](tbegin(pdae), pdae.ics.q, pdae.ics.p)

    ldae_equs = functions(ldae)
    ldae_invs = invariants(ldae)
    @test_nowarn ldae_equs[:ϑ](tbegin(ldae), ldae.ics.q, ldae.ics.λ, zero(ldae.ics.q))
    @test_nowarn ldae_equs[:f](tbegin(ldae), ldae.ics.q, ldae.ics.λ, zero(ldae.ics.p))
    @test_nowarn ldae_equs[:g](tbegin(ldae), ldae.ics.q, zero(ldae.ics.q), ldae.ics.λ, zero(ldae.ics.p))
    @test_nowarn ldae_equs[:ḡ](tbegin(ldae), ldae.ics.q, zero(ldae.ics.q), ldae.ics.λ, zero(ldae.ics.p))
    @test_nowarn ldae_equs[:ϕ](tbegin(ldae), ldae.ics.q, ldae.ics.p, zero(ldae.ics.λ))
    @test_nowarn ldae_equs[:ψ](tbegin(ldae), ldae.ics.q, ldae.ics.p, zero(ldae.ics.q), zero(ldae.ics.p), zero(ldae.ics.λ))
    @test_nowarn ldae_equs[:v̄](tbegin(ldae), ldae.ics.q, zero(ldae.ics.q))
    @test_nowarn ldae_equs[:f̄](tbegin(ldae), ldae.ics.q, ldae.ics.λ, zero(ldae.ics.p))
    @test_nowarn ldae_invs[:h](tbegin(ldae), ldae.ics.q, zero(iode.ics.q))


    int = Integrator(ode, TableauGauss(2))
    sol = integrate(ode, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2))
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2))
    sol = integrate(iode, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2))
    sol = integrate(idae, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2))
    sol = integrate(idae, int, nt)
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4
end
