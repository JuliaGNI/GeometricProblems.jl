using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2d
using GeometricProblems.LotkaVolterra2d: reference_solution


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

    # ode_equs = functions(ode)
    # ode_invs = invariants(ode)
    # @test_nowarn ode_equs[:v](zero(ode.ics.q), tbegin(ode), ode.ics.q)
    # @test_nowarn ode_invs[:h](tbegin(ode), ode.ics.q)

    # hode_equs = functions(hode)
    # @test_nowarn hode_equs[:v](zero(hode.ics.q), tbegin(hode), hode.ics.q, hode.ics.p)
    # @test_nowarn hode_equs[:f](zero(hode.ics.p), tbegin(hode), hode.ics.q, hode.ics.p)
    # @test_nowarn hode_equs[:h](tbegin(hode), hode.ics.q, hode.ics.p)

    # iode_equs = functions(iode)
    # iode_invs = invariants(iode)
    # @test_nowarn iode_equs[:ϑ](zero(iode.ics.q), tbegin(iode), iode.ics.q, iode.ics.p)
    # @test_nowarn iode_equs[:f](zero(iode.ics.p), tbegin(iode), iode.ics.q, iode.ics.p)
    # @test_nowarn iode_equs[:v̄](zero(iode.ics.q), tbegin(iode), iode.ics.q)
    # @test_nowarn iode_equs[:f̄](zero(iode.ics.p), tbegin(iode), iode.ics.q, iode.ics.p)
    # @test_nowarn iode_invs[:h](tbegin(iode), iode.ics.q, zero(iode.ics.q))

    # pode_equs = functions(pode)
    # pode_invs = invariants(pode)
    # @test_nowarn pode_equs[:v](zero(pode.ics.q), tbegin(pode), pode.ics.q, pode.ics.p)
    # @test_nowarn pode_equs[:f](zero(pode.ics.p), tbegin(pode), pode.ics.q, pode.ics.p)
    # @test_nowarn pode_invs[:h](tbegin(pode), pode.ics.q, pode.ics.p)

    # lode_equs = functions(lode)
    # lode_invs = invariants(lode)
    # @test_nowarn lode_equs[:ϑ](zero(lode.ics.q), tbegin(lode), lode.ics.q, lode.ics.p)
    # @test_nowarn lode_equs[:f](zero(lode.ics.p), tbegin(lode), lode.ics.q, lode.ics.p)
    # @test_nowarn lode_equs[:v̄](zero(lode.ics.q), tbegin(lode), lode.ics.q)
    # @test_nowarn lode_equs[:f̄](zero(lode.ics.p), tbegin(lode), lode.ics.q, lode.ics.p)
    # @test_nowarn lode_invs[:h](tbegin(lode), lode.ics.q, zero(iode.ics.q))

    # dae_equs = functions(dae)
    # dae_invs = invariants(dae)
    # @test_nowarn dae_equs[:v](zero(dae.ics.q), tbegin(dae), dae.ics.q)
    # @test_nowarn dae_equs[:u](zero(dae.ics.q), tbegin(dae), dae.ics.q, dae.ics.λ)
    # @test_nowarn dae_equs[:ϕ](zero(dae.ics.λ), tbegin(dae), dae.ics.q)
    # @test_nowarn dae_equs[:v̄](zero(dae.ics.q), tbegin(dae), dae.ics.q)
    # @test_nowarn dae_invs[:h](tbegin(dae), dae.ics.q)

    # hdae_equs = functions(hdae)
    # @test_nowarn hdae_equs[:v](zero(hdae.ics.q), tbegin(hdae), hdae.ics.q, hdae.ics.p)
    # @test_nowarn hdae_equs[:f](zero(hdae.ics.p), tbegin(hdae), hdae.ics.q, hdae.ics.p)
    # @test_nowarn hdae_equs[:u](zero(hdae.ics.q), tbegin(hdae), hdae.ics.q, hdae.ics.p, hdae.ics.λ)
    # @test_nowarn hdae_equs[:g](zero(hdae.ics.p), tbegin(hdae), hdae.ics.q, hdae.ics.p, hdae.ics.λ)
    # @test_nowarn hdae_equs[:ϕ](zero(hdae.ics.λ), tbegin(hdae), hdae.ics.q, hdae.ics.p)
    # @test_nowarn hdae_equs[:v̄](zero(hdae.ics.q), tbegin(hdae), hdae.ics.q, hdae.ics.p)
    # @test_nowarn hdae_equs[:f̄](zero(hdae.ics.p), tbegin(hdae), hdae.ics.q, hdae.ics.p)
    # @test_nowarn hdae_equs[:h](tbegin(hdae), hdae.ics.q, hdae.ics.p)

    # idae_equs = functions(idae)
    # idae_invs = invariants(idae)
    # @test_nowarn idae_equs[:ϑ](zero(idae.ics.q), tbegin(idae), idae.ics.q, idae.ics.p)
    # @test_nowarn idae_equs[:f](zero(idae.ics.p), tbegin(idae), idae.ics.q, idae.ics.p)
    # @test_nowarn idae_equs[:u](zero(idae.ics.q), tbegin(idae), idae.ics.q, idae.ics.p, idae.ics.λ)
    # @test_nowarn idae_equs[:g](zero(idae.ics.p), tbegin(idae), idae.ics.q, idae.ics.p, idae.ics.λ)
    # @test_nowarn idae_equs[:ϕ](zero(idae.ics.λ), tbegin(idae), idae.ics.q, idae.ics.p)
    # @test_nowarn idae_equs[:v̄](zero(idae.ics.q), tbegin(idae), idae.ics.q)
    # @test_nowarn idae_equs[:f̄](zero(idae.ics.p), tbegin(idae), idae.ics.q, idae.ics.p)
    # @test_nowarn idae_invs[:h](tbegin(idae), idae.ics.q, zero(iode.ics.q))

    # pdae_equs = functions(pdae)
    # pdae_invs = invariants(pdae)
    # @test_nowarn pdae_equs[:v](zero(pdae.ics.q), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    # @test_nowarn pdae_equs[:f](zero(pdae.ics.p), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    # @test_nowarn pdae_equs[:u](zero(pdae.ics.q), tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ)
    # @test_nowarn pdae_equs[:g](zero(pdae.ics.p), tbegin(pdae), pdae.ics.q, pdae.ics.p, pdae.ics.λ)
    # @test_nowarn pdae_equs[:ϕ](zero(pdae.ics.λ), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    # @test_nowarn pdae_equs[:v̄](zero(pdae.ics.q), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    # @test_nowarn pdae_equs[:f̄](zero(pdae.ics.p), tbegin(pdae), pdae.ics.q, pdae.ics.p)
    # @test_nowarn pdae_invs[:h](tbegin(pdae), pdae.ics.q, pdae.ics.p)

    # ldae_equs = functions(ldae)
    # ldae_invs = invariants(ldae)
    # @test_nowarn ldae_equs[:ϑ](zero(ldae.ics.q), tbegin(ldae), ldae.ics.q, ldae.ics.λ)
    # @test_nowarn ldae_equs[:f](zero(ldae.ics.p), tbegin(ldae), ldae.ics.q, ldae.ics.λ)
    # @test_nowarn ldae_equs[:g](zero(ldae.ics.p), tbegin(ldae), ldae.ics.q, zero(ldae.ics.q), ldae.ics.λ)
    # @test_nowarn ldae_equs[:ḡ](zero(ldae.ics.p), tbegin(ldae), ldae.ics.q, zero(ldae.ics.q), ldae.ics.λ)
    # @test_nowarn ldae_equs[:ϕ](zero(ldae.ics.λ), tbegin(ldae), ldae.ics.q, ldae.ics.p)
    # @test_nowarn ldae_equs[:ψ](zero(ldae.ics.λ), tbegin(ldae), ldae.ics.q, ldae.ics.p, zero(ldae.ics.q), zero(ldae.ics.p))
    # @test_nowarn ldae_equs[:v̄](zero(ldae.ics.q), tbegin(ldae), ldae.ics.q)
    # @test_nowarn ldae_equs[:f̄](zero(ldae.ics.p), tbegin(ldae), ldae.ics.q, ldae.ics.λ)
    # @test_nowarn ldae_invs[:h](tbegin(ldae), ldae.ics.q, zero(iode.ics.q))


    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, reference_solution) < 5E-4
end
