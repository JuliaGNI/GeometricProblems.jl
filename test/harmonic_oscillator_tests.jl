using Test
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillator: reference_solution, reference_solution_q, reference_solution_p
using GeometricSolutions


@testset "$(rpad("Harmonic Oscillator",80))" begin

    @test_nowarn odeproblem()
    @test_nowarn hodeproblem()
    @test_nowarn iodeproblem()
    @test_nowarn lodeproblem()
    @test_nowarn podeproblem()
    @test_nowarn sodeproblem()

    @test_nowarn degenerate_iodeproblem()
    @test_nowarn degenerate_lodeproblem()

    @test_nowarn daeproblem()
    @test_nowarn hdaeproblem()
    @test_nowarn idaeproblem()
    @test_nowarn ldaeproblem()
    @test_nowarn pdaeproblem()

    @test_nowarn odeensemble()
    @test_nowarn podeensemble()
    @test_nowarn hodeensemble()


    ode = odeproblem()
    iode = degenerate_iodeproblem()
    pode = podeproblem()
    hode = hodeproblem()
    ref  = exact_solution(ode)
    
    
    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-4

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-4

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-4


    sol = exact_solution(ode)
    @test sol.q[end] == reference_solution

    sol = exact_solution(pode)
    @test sol.q[end] == [reference_solution_q]
    @test sol.p[end] == [reference_solution_p]

    sol = exact_solution(hode)
    @test sol.q[end] == [reference_solution_q]
    @test sol.p[end] == [reference_solution_p]

end


# The mass `m` has to be threaded consistently through every formulation. These tests use a
# non-unit mass and a non-zero initial velocity, the combination under which the earlier
# `m`-agnostic code silently disagreed with the Hamiltonian and the exact solution.

@testset "$(rpad("Harmonic Oscillator (non-unit mass)",80))" begin
    m = 2.0
    params = (m = m, k = 0.5, ω = sqrt(0.5 / m))

    # The ODE state is (q, v); the exact solution converts to and from the momentum convention.
    x₀ = [0.7, 0.3]
    ode = odeproblem(x₀; timespan = (0.0, 5.0), timestep = 0.001, parameters = params)
    sol = integrate(ode, Gauss(4))
    @test sol.q[end] ≈ HarmonicOscillator.exact_solution(5.0, x₀, 0.0, params) atol = 1e-12

    # The Hamiltonian formulation uses momenta directly.
    q₀, p₀ = [0.7], [0.3]
    hsol = integrate(hodeproblem(q₀, p₀; timespan = (0.0, 5.0), timestep = 0.001, parameters = params), Gauss(4))
    @test hsol.q[end][1] ≈ HarmonicOscillator.exact_solution_q(5.0, q₀, p₀, 0.0, params) atol = 1e-12
    @test hsol.p[end][1] ≈ HarmonicOscillator.exact_solution_p(5.0, q₀, p₀, 0.0, params) atol = 1e-12

    # The degenerate one-form is ϑ = (m·q₂, 0), so the default initial momentum of the degenerate
    # problems must carry the mass. It used to be built from a mass-free `ϑ`, which only agreed
    # with the problem's own `ϑ` at m = 1 (or for the default, zero-velocity initial condition).
    ϑ₀ = zeros(2)
    dprob = degenerate_iodeproblem(x₀; parameters = params)
    functions(dprob).ϑ(ϑ₀, 0.0, x₀, zeros(2), params)
    @test dprob.ics.p ≈ ϑ₀
    @test dprob.ics.p ≈ [m * x₀[2], 0.0]

    # Ω = dϑ of the degenerate one-form, in the package-wide convention Ωᵢⱼ = ∂ϑᵢ/∂qⱼ - ∂ϑⱼ/∂qᵢ.
    Ω = zeros(2, 2)
    HarmonicOscillator.degenerate_ω!(Ω, 0.0, x₀, params)
    @test Ω ≈ [0.0 m; -m 0.0]

    # Ω v = -∇H must reproduce the degenerate velocity field.
    v = zeros(2)
    HarmonicOscillator.degenerate_oscillator_iode_v(v, 0.0, x₀, zeros(2), params)
    ∇H = [params.k * x₀[1], m * x₀[2]]
    @test Ω * v ≈ -∇H atol = 1e-12
end


# The energy constraint of the differential-algebraic forms is ϕ = H - H₀. Its reference energy H₀
# has to come from the problem's own initial condition: it used to be read from the module-level
# default, so ϕ did not vanish at t₀ for any other initial data.

@testset "$(rpad("Harmonic Oscillator (DAE energy constraint)",80))" begin
    params = default_parameters()
    q₀, p₀ = [1.5], [0.7]          # deliberately *not* the module default (q₀ = [0.5], p₀ = [0.0])
    x₀ = [1.5, 0.7]
    ϕ = zeros(1)
    t₀ = 0.0

    functions(daeproblem(x₀, [0.0])).ϕ(ϕ, t₀, x₀, params)
    @test ϕ[1] ≈ 0 atol = 1e-14

    for prob in (pdaeproblem(q₀, p₀), hdaeproblem(q₀, p₀))
        functions(prob).ϕ(ϕ, t₀, q₀, p₀, params)
        @test ϕ[1] ≈ 0 atol = 1e-14
    end

    # The implicit/variational forms evaluate ϕ with an extra velocity slot.
    for prob in (idaeproblem(q₀, p₀), ldaeproblem(q₀, p₀))
        functions(prob).ϕ(ϕ, t₀, q₀, [0.2], p₀, params)
        @test ϕ[1] ≈ 0 atol = 1e-14
    end
end
