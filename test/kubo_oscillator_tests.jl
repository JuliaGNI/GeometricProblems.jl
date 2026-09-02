using GeometricEquations
using GeometricProblems.KuboOscillator
using Test

# The Kubo oscillator is a unit-frequency harmonic oscillator driven by multiplicative
# (Stratonovich) noise. Here we check that the SDE / PSDE / SPSDE problems build with the current
# GeometricEquations noise-object API and that their equation functions evaluate to the expected
# Kubo drift and diffusion terms. `functions(prob)` returns the raw functions, so we pass the
# parameters explicitly.

@testset "$(rpad("Kubo Oscillator",80))" begin
    params = KuboOscillator.default_parameters()
    ν = params.ν

    # --- SDE ---
    sde = sdeproblem()
    @test sde isa SDEProblem
    f = functions(sde)

    q = [0.5, 0.3]
    v = zero(q)
    f.v(v, 0.0, q, params)
    @test v ≈ [q[2], -q[1]]                          # drift: unit harmonic oscillator
    B = zeros(eltype(q), 2, 1)
    f.B(B, 0.0, q, params)
    @test B ≈ reshape([ν * q[2], -ν * q[1]], 2, 1)   # multiplicative diffusion

    # single-IC builders are problems; the multi-IC builders create ensembles
    @test sdeensemble() isa GeometricEquations.EnsembleProblem
    @test psdeensemble() isa GeometricEquations.EnsembleProblem
    @test spsdeensemble() isa GeometricEquations.EnsembleProblem

    # --- PSDE ---
    psde = psdeproblem()
    @test psde isa PSDEProblem
    fp = functions(psde)

    qq = [0.5]
    pp = [0.3]
    v = zero(qq)
    fp.v(v, 0.0, qq, pp, params)
    @test v ≈ [pp[1]]
    ff = zero(pp)
    fp.f(ff, 0.0, qq, pp, params)
    @test ff ≈ [-qq[1]]
    B = zeros(1, 1)
    fp.B(B, 0.0, qq, pp, params)
    @test B ≈ reshape([ν * pp[1]], 1, 1)
    G = zeros(1, 1)
    fp.G(G, 0.0, qq, pp, params)
    @test G ≈ reshape([-ν * qq[1]], 1, 1)

    # --- SPSDE ---
    spsde = spsdeproblem()
    @test spsde isa SPSDEProblem
    fs = functions(spsde)

    v = zero(qq)
    fs.v(v, 0.0, qq, pp, params)
    @test v ≈ [pp[1]]
    f1 = zero(pp)
    fs.f1(f1, 0.0, qq, pp, params)
    @test f1 ≈ [-qq[1]]
    f2 = zero(pp)
    fs.f2(f2, 0.0, qq, pp, params)
    @test f2 ≈ [0.0]
    B = zeros(1, 1)
    fs.B(B, 0.0, qq, pp, params)
    @test B ≈ reshape([ν * pp[1]], 1, 1)
    G1 = zeros(1, 1)
    fs.G1(G1, 0.0, qq, pp, params)
    @test G1 ≈ reshape([-ν * qq[1]], 1, 1)
    G2 = zeros(1, 1)
    fs.G2(G2, 0.0, qq, pp, params)
    @test G2 ≈ reshape([0.0], 1, 1)

    # ODE (deterministic drift only)
    @test odeproblem() isa ODEProblem
end

# The damped variant is the only one here with non-zero f2 and G2, so it is what exercises the
# split half of an SPSDE. The checks below pin both the individual terms and the agreement
# between the split and unsplit formulations, which is the property that makes the pair useful.

@testset "$(rpad("Kubo Oscillator (damped)",80))" begin
    params = KuboOscillator.damped_parameters()
    ν = params.ν
    γ = params.γ

    qq = [2.0]
    pp = [0.3]

    # --- PSDE: the damping folded into f and G ---
    psde = damped_psdeproblem()
    @test psde isa PSDEProblem
    fp = functions(psde)

    v = zero(qq)
    fp.v(v, 0.0, qq, pp, params)
    @test v ≈ [pp[1]]
    f = zero(pp)
    fp.f(f, 0.0, qq, pp, params)
    @test f ≈ [-qq[1] - γ * pp[1]]
    B = zeros(1, 1)
    fp.B(B, 0.0, qq, pp, params)
    @test B ≈ reshape([ν * pp[1]], 1, 1)
    G = zeros(1, 1)
    fp.G(G, 0.0, qq, pp, params)
    @test G ≈ reshape([-ν * qq[1] - ν * γ * pp[1]], 1, 1)

    # --- SPSDE: the same dynamics, Hamiltonian part in f1/G1 and damping in f2/G2 ---
    spsde = damped_spsdeproblem()
    @test spsde isa SPSDEProblem
    fs = functions(spsde)

    f1 = zero(pp)
    fs.f1(f1, 0.0, qq, pp, params)
    @test f1 ≈ [-qq[1]]
    f2 = zero(pp)
    fs.f2(f2, 0.0, qq, pp, params)
    @test f2 ≈ [-γ * pp[1]]
    @test f2 ≉ [0.0]                                 # unlike the undamped problem
    G1 = zeros(1, 1)
    fs.G1(G1, 0.0, qq, pp, params)
    @test G1 ≈ reshape([-ν * qq[1]], 1, 1)
    G2 = zeros(1, 1)
    fs.G2(G2, 0.0, qq, pp, params)
    @test G2 ≈ reshape([-ν * γ * pp[1]], 1, 1)
    @test G2 ≉ zeros(1, 1)

    # the split and unsplit formulations describe the same dynamics
    @test f1 + f2 ≈ f
    @test G1 + G2 ≈ G

    # The damped builders carry their own default timespan: over the undamped problems' 0.1 the
    # energy decays by 5e-06 relative and the damped problem is indistinguishable from the
    # undamped one.
    @test timespan(psde) == KuboOscillator.DEFAULT_DAMPED_TIMESPAN
    @test timespan(spsde) == KuboOscillator.DEFAULT_DAMPED_TIMESPAN
    @test timestep(psde) == KuboOscillator.DEFAULT_DAMPED_TIMESTEP
    @test timespan(psde)[end] > KuboOscillator.DEFAULT_TIMESPAN[end]
end

# `exact_solution` is a pathwise reference: given the Wiener increment it returns the solution
# exactly, which is what a convergence-order measurement needs. `exact_mean_energy` is the
# corresponding average over the noise. The full derivation checks are in
# `scripts/verify_kubo_exact_solution.jl`; what follows pins the interface and the invariants.

@testset "$(rpad("Kubo Oscillator (exact solution)",80))" begin
    undamped = KuboOscillator.default_parameters()
    damped = KuboOscillator.damped_parameters()

    q₀, p₀, t₀ = 2.0, 0.3, 1.5
    H(q, p) = (q^2 + p^2) / 2

    for params in (undamped, damped)
        # at t = t₀ along the path that has not moved, the solution is the initial condition
        @test all(exact_solution(t₀, 0.0, q₀, p₀, t₀, params) .≈ (q₀, p₀))

        # only the elapsed time enters, so shifting t and t₀ together changes nothing ...
        @test all(exact_solution(t₀ + 0.7, 0.4, q₀, p₀, t₀, params) .≈
                  exact_solution(0.7, 0.4, q₀, p₀, 0.0, params))

        # ... and t₀ is not simply ignored
        @test !all(exact_solution(t₀ + 0.7, 0.4, q₀, p₀, t₀, params) .≈
                   exact_solution(t₀ + 0.7, 0.4, q₀, p₀, 0.0, params))

        # the component accessors and the vector forms agree with the scalar form
        q, p = exact_solution(1.0, 0.4, q₀, p₀, t₀, params)
        @test exact_solution_q(1.0, 0.4, q₀, p₀, t₀, params) ≈ q
        @test exact_solution_p(1.0, 0.4, q₀, p₀, t₀, params) ≈ p
        @test all(exact_solution(1.0, 0.4, [q₀], [p₀], t₀, params) .≈ ([q], [p]))
        @test exact_solution(1.0, 0.4, [q₀, p₀], t₀, params) ≈ [q, p]
    end

    # without damping the energy is conserved exactly along every sample path
    for (t, W) in ((0.0, 0.0), (1.3, -0.8), (7.0, 2.4))
        @test H(exact_solution(t, W, q₀, p₀, 0.0, undamped)...) ≈ H(q₀, p₀)
    end

    # with damping it decays
    @test H(exact_solution(100.0, 0.0, q₀, p₀, 0.0, damped)...) < H(q₀, p₀)

    # the closed form is the underdamped solution: |γ| = 2 is the case that would otherwise
    # return NaN silently, |γ| > 2 the one with an imaginary frequency
    @test_throws DomainError exact_solution(1.0, 0.0, q₀, p₀, 0.0, (ν = 0.5, γ = 2.0))
    @test_throws DomainError exact_solution(1.0, 0.0, q₀, p₀, 0.0, (ν = 0.5, γ = 2.5))
    @test_throws DomainError exact_mean_energy(1.0, q₀, p₀, 0.0, (ν = 0.5, γ = 2.0))
end

@testset "$(rpad("Kubo Oscillator (mean energy)",80))" begin
    undamped = KuboOscillator.default_parameters()
    damped = KuboOscillator.damped_parameters()

    q₀, p₀, t₀ = 2.0, 0.3, 1.5
    H(q, p) = (q^2 + p^2) / 2
    H₀ = H(q₀, p₀)

    # without damping E(H) is the constant initial energy
    for t in (0.0, 1.0, 25.0)
        @test exact_mean_energy(t, q₀, p₀, 0.0, undamped) ≈ H₀
    end

    # with damping it starts there and decays; only the elapsed time enters
    @test exact_mean_energy(0.0, q₀, p₀, 0.0, damped) ≈ H₀
    @test exact_mean_energy(1000.0, q₀, p₀, 0.0, damped) < H₀
    @test exact_mean_energy(t₀ + 3.0, q₀, p₀, t₀, damped) ≈
          exact_mean_energy(3.0, q₀, p₀, 0.0, damped)
    @test exact_mean_energy(3.0, [q₀], [p₀], 0.0, damped) ≈
          exact_mean_energy(3.0, q₀, p₀, 0.0, damped)

    # E(H) is an average of the pathwise energy over the Wiener increment W ~ N(0, t - t₀), so
    # the closed form has to agree with a quadrature of `exact_solution` against that density.
    # The integrand is smooth and Gaussian-weighted, so the trapezoidal rule converges
    # spectrally and this is far more accurate than the tolerance asserted below.
    function quadrature_mean_energy(t, q₀, p₀, params; nodes = 2001, width = 8)
        σ = sqrt(t)
        z = range(-width * σ, width * σ; length = nodes)
        density = exp.(-(z .^ 2) ./ (2 * σ^2)) ./ (σ * sqrt(2π))
        energy = [H(exact_solution(t, W, q₀, p₀, 0.0, params)...) for W in z]
        step(z) * sum(density .* energy)
    end

    @testset "quadrature (ν = $(params.ν), γ = $(get(params, :γ, 0.0)), t = $t)" for params in (
            undamped, damped, (ν = 0.2, γ = 1.0)),
        t in (0.5, 2.0, 10.0)

        @test quadrature_mean_energy(t, q₀, p₀, params)≈exact_mean_energy(
            t, q₀, p₀, 0.0, params) rtol=1e-8
    end
end
