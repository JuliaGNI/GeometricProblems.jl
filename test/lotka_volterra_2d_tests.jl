using SimpleSolvers
using Test
using GeometricIntegrators
using GeometricIntegrators.Integrators.VPRK
using GeometricIntegrators.Utils
using GeometricProblems.LotkaVolterra2d
using GeometricProblems.LotkaVolterra2d: Δt, nt, reference_solution

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

    ode_equs = get_function_tuple(ode)
    @test_nowarn ode_equs[:v](ode.t₀, ode.q₀[begin], zero(ode.q₀[begin]))
    @test_nowarn ode_equs[:h](ode.t₀, ode.q₀[begin])

    hode_equs = get_function_tuple(hode)
    @test_nowarn hode_equs[:v](hode.t₀, hode.q₀[begin], hode.p₀[begin], zero(hode.q₀[begin]))
    @test_nowarn hode_equs[:f](hode.t₀, hode.q₀[begin], hode.p₀[begin], zero(hode.p₀[begin]))
    @test_nowarn hode_equs[:h](hode.t₀, hode.q₀[begin], hode.p₀[begin])

    iode_equs = get_function_tuple(iode)
    @test_nowarn iode_equs[:ϑ](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.p₀[begin]))
    @test_nowarn iode_equs[:v̄](iode.t₀, iode.q₀[begin], zero(iode.q₀[begin]))
    @test_nowarn iode_equs[:f̄](iode.t₀, iode.q₀[begin], iode.p₀[begin], zero(iode.p₀[begin]))
    @test_nowarn iode_equs[:h](iode.t₀, iode.q₀[begin])

    pode_equs = get_function_tuple(pode)
    @test_nowarn pode_equs[:v](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.q₀[begin]))
    @test_nowarn pode_equs[:f](pode.t₀, pode.q₀[begin], pode.p₀[begin], zero(pode.p₀[begin]))
    @test_nowarn pode_equs[:h](pode.t₀, pode.q₀[begin], pode.p₀[begin])

    lode_equs = get_function_tuple(lode)
    @test_nowarn lode_equs[:ϑ](lode.t₀, lode.q₀[begin], lode.p₀[begin], zero(lode.q₀[begin]))
    @test_nowarn lode_equs[:f](lode.t₀, lode.q₀[begin], lode.p₀[begin], zero(lode.p₀[begin]))
    @test_nowarn lode_equs[:v̄](lode.t₀, lode.q₀[begin], zero(lode.q₀[begin]))
    @test_nowarn lode_equs[:f̄](lode.t₀, lode.q₀[begin], lode.p₀[begin], zero(lode.p₀[begin]))
    @test_nowarn lode_equs[:h](lode.t₀, lode.q₀[begin])

    dae_equs = get_function_tuple(dae)
    @test_nowarn dae_equs[:v](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:u](dae.t₀, dae.q₀[begin], dae.λ₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:ϕ](dae.t₀, dae.q₀[begin], zero(dae.λ₀[begin]))
    @test_nowarn dae_equs[:v̄](dae.t₀, dae.q₀[begin], zero(dae.q₀[begin]))
    @test_nowarn dae_equs[:h](dae.t₀, dae.q₀[begin])

    hdae_equs = get_function_tuple(hdae)
    @test_nowarn hdae_equs[:v](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.q₀[begin]))
    @test_nowarn hdae_equs[:f](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.p₀[begin]))
    @test_nowarn hdae_equs[:u](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], hdae.λ₀[begin], zero(hdae.q₀[begin]))
    @test_nowarn hdae_equs[:g](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], hdae.λ₀[begin], zero(hdae.p₀[begin]))
    @test_nowarn hdae_equs[:ϕ](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.λ₀[begin]))
    @test_nowarn hdae_equs[:v̄](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.q₀[begin]))
    @test_nowarn hdae_equs[:f̄](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin], zero(hdae.p₀[begin]))
    @test_nowarn hdae_equs[:h](hdae.t₀, hdae.q₀[begin], hdae.p₀[begin])

    idae_equs = get_function_tuple(idae)
    @test_nowarn idae_equs[:ϑ](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:u](idae.t₀, idae.q₀[begin], idae.p₀[begin], idae.λ₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:g](idae.t₀, idae.q₀[begin], idae.p₀[begin], idae.λ₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:ϕ](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.λ₀[begin]))
    @test_nowarn idae_equs[:v̄](idae.t₀, idae.q₀[begin], zero(idae.q₀[begin]))
    @test_nowarn idae_equs[:f̄](idae.t₀, idae.q₀[begin], idae.p₀[begin], zero(idae.p₀[begin]))
    @test_nowarn idae_equs[:h](idae.t₀, idae.q₀[begin])

    pdae_equs = get_function_tuple(pdae)
    @test_nowarn pdae_equs[:v](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:u](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:g](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], pdae.λ₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:ϕ](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.λ₀[begin]))
    @test_nowarn pdae_equs[:v̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.q₀[begin]))
    @test_nowarn pdae_equs[:f̄](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin], zero(pdae.p₀[begin]))
    @test_nowarn pdae_equs[:h](pdae.t₀, pdae.q₀[begin], pdae.p₀[begin])

    ldae_equs = get_function_tuple(ldae)
    @test_nowarn ldae_equs[:ϑ](ldae.t₀, ldae.q₀[begin], ldae.λ₀[begin], zero(ldae.q₀[begin]))
    @test_nowarn ldae_equs[:f](ldae.t₀, ldae.q₀[begin], ldae.λ₀[begin], zero(ldae.p₀[begin]))
    @test_nowarn ldae_equs[:g](ldae.t₀, ldae.q₀[begin], ldae.λ₀[begin], zero(ldae.p₀[begin]))
    @test_nowarn ldae_equs[:ḡ](ldae.t₀, ldae.q₀[begin], ldae.λ₀[begin], zero(ldae.p₀[begin]))
    @test_nowarn ldae_equs[:ϕ](ldae.t₀, ldae.q₀[begin], ldae.p₀[begin], zero(ldae.λ₀[begin]))
    @test_nowarn ldae_equs[:ψ](ldae.t₀, ldae.q₀[begin], ldae.p₀[begin], zero(ldae.q₀[begin]), zero(ldae.p₀[begin]), zero(ldae.λ₀[begin]))
    @test_nowarn ldae_equs[:v̄](ldae.t₀, ldae.q₀[begin], zero(ldae.q₀[begin]))
    @test_nowarn ldae_equs[:f̄](ldae.t₀, ldae.q₀[begin], ldae.λ₀[begin], zero(ldae.p₀[begin]))
    @test_nowarn ldae_equs[:h](ldae.t₀, ldae.q₀[begin])


    int = Integrator(ode, TableauGauss(2), Δt)
    sol = integrate(ode, int, nt)
    @test rel_err(sol.q, reference_solution) < 5E-4

    int = IntegratorVPRKpMidpoint(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 5E-4

    int = IntegratorVPRKpSymmetric(iode, TableauVPGLRK(2), Δt)
    sol = integrate(iode, int, nt)
    @test rel_err(sol.q, reference_solution) < 5E-4

    int = Integrator(idae, TableauVSPARKGLRKpMidpoint(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 5E-4

    int = Integrator(idae, TableauVSPARKGLRKpSymmetric(2), Δt)
    sol = integrate(idae, int, nt)
    @test rel_err(sol.q, reference_solution) < 5E-4
end
