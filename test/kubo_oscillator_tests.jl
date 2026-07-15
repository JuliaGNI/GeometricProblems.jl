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
    v = zero(q); f.v(v, 0.0, q, params)
    @test v ≈ [q[2], -q[1]]                          # drift: unit harmonic oscillator
    B = zeros(eltype(q), 2, 1); f.B(B, 0.0, q, params)
    @test B ≈ reshape([ν * q[2], -ν * q[1]], 2, 1)   # multiplicative diffusion

    # single-IC builders are problems; the multi-IC builders create ensembles
    @test sdeensemble() isa GeometricEquations.EnsembleProblem
    @test psdeensemble() isa GeometricEquations.EnsembleProblem
    @test spsdeensemble() isa GeometricEquations.EnsembleProblem

    # --- PSDE ---
    psde = psdeproblem()
    @test psde isa PSDEProblem
    fp = functions(psde)

    qq = [0.5]; pp = [0.3]
    v = zero(qq); fp.v(v, 0.0, qq, pp, params); @test v ≈ [pp[1]]
    ff = zero(pp); fp.f(ff, 0.0, qq, pp, params); @test ff ≈ [-qq[1]]
    B = zeros(1, 1); fp.B(B, 0.0, qq, pp, params); @test B ≈ reshape([ν * pp[1]], 1, 1)
    G = zeros(1, 1); fp.G(G, 0.0, qq, pp, params); @test G ≈ reshape([-ν * qq[1]], 1, 1)

    # --- SPSDE ---
    spsde = spsdeproblem()
    @test spsde isa SPSDEProblem
    fs = functions(spsde)

    v = zero(qq); fs.v(v, 0.0, qq, pp, params); @test v ≈ [pp[1]]
    f1 = zero(pp); fs.f1(f1, 0.0, qq, pp, params); @test f1 ≈ [-qq[1]]
    f2 = zero(pp); fs.f2(f2, 0.0, qq, pp, params); @test f2 ≈ [0.0]
    B = zeros(1, 1); fs.B(B, 0.0, qq, pp, params); @test B ≈ reshape([ν * pp[1]], 1, 1)
    G1 = zeros(1, 1); fs.G1(G1, 0.0, qq, pp, params); @test G1 ≈ reshape([-ν * qq[1]], 1, 1)
    G2 = zeros(1, 1); fs.G2(G2, 0.0, qq, pp, params); @test G2 ≈ reshape([0.0], 1, 1)

    # ODE (deterministic drift only)
    @test odeproblem() isa ODEProblem
end
