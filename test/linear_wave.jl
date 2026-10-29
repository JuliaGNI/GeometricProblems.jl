using GeometricEquations: HODEEnsemble, HODEProblem, LODEEnsemble, LODEProblem, functions,
                          initialguess
using GeometricIntegrators
using GeometricProblems.LinearWave
using GeometricSolutions
using LinearAlgebra
using Test

include("helpers/integrate_quietly.jl")

const LW = LinearWave

# Regression test for the linear wave equation. Previously `lodeproblem` called a lowercase
# `lodeproblem(::LagrangianSystem, …)` (itself) and threw a `MethodError`; it must build a proper
# `LODEProblem`, mirroring `hodeproblem`.
#
# The two default-size constructions below used to generate the equations of motion symbolically for
# Ñ = 256, which cost 155 s for `lodeproblem` alone. They are now instantaneous, and the identity
# assertions are what keep them that way: if the default ever routes through EulerLagrange again,
# these fail rather than merely getting slow.

@testset "$(rpad("Linear Wave",80))" begin
    hode = @test_nowarn hodeproblem()
    lode = @test_nowarn lodeproblem()

    @test hode isa HODEProblem
    @test lode isa LODEProblem

    @test functions(hode).v === LW.linear_wave_v
    @test functions(hode).f === LW.linear_wave_f
    @test functions(hode).h === LW.hamiltonian

    @test functions(lode).ϑ === LW.linear_wave_ϑ
    @test functions(lode).f === LW.linear_wave_f
    @test functions(lode).g === LW.linear_wave_g
    @test functions(lode).ω === LW.ω!
    @test functions(lode).l === LW.lagrangian

    # was `GeometricEquations._lode_default_v̄`, which returns `nothing` without writing anything
    @test initialguess(lode).v === LW.v̄
end

@testset "$(rpad("Linear Wave: hand-written matches generated",80))" begin
    # Small on purpose: `lagrangian_system` grows as n^2.5 and the generated `ω` as n^2.
    N = 5
    n = N
    params = LW.default_parameters()

    q = LW.initial_condition(LW.μ̃, n).q
    p = LW.initial_condition(LW.μ̃, n).p
    v = copy(p)                 # the mass matrix is the identity, so q̇ = p

    hamiltonian_functions = functions(LW.hamiltonian_system(N, params))
    lagrangian_functions = functions(LW.lagrangian_system(N, params))

    generated = zeros(n)
    handwritten = zeros(n)

    # ∂H/∂p = p
    hamiltonian_functions.v(generated, 0.0, q, p, params)
    LW.linear_wave_v(handwritten, 0.0, q, p, params)
    @test generated ≈ handwritten
    @test handwritten == p

    # -∂H/∂q, and the identical ∂L/∂q
    hamiltonian_functions.f(generated, 0.0, q, p, params)
    LW.linear_wave_f(handwritten, 0.0, q, p, params)
    @test generated ≈ handwritten
    @test any(!iszero, handwritten)

    lagrangian_functions.f(generated, 0.0, q, v, params)
    LW.linear_wave_f(handwritten, 0.0, q, v, params)
    @test generated ≈ handwritten

    # ϑ = ∂L/∂q̇ = q̇
    lagrangian_functions.ϑ(generated, 0.0, q, v, params)
    LW.linear_wave_ϑ(handwritten, 0.0, q, v, params)
    @test generated ≈ handwritten
    @test any(!iszero, handwritten)

    # g = λ, in both arities
    λ = collect(range(0.1, 0.7, length = n))
    lagrangian_functions.g(generated, 0.0, q, v, λ, params)
    LW.linear_wave_g(handwritten, 0.0, q, v, λ, params)
    @test generated ≈ handwritten ≈ λ
    LW.linear_wave_g(handwritten, 0.0, q, λ, params)
    @test handwritten ≈ λ

    @test lagrangian_functions.L(0.0, q, v, params) ≈ lagrangian(0.0, q, v, params, N)
    @test hamiltonian_functions.H(0.0, q, p, params) ≈ hamiltonian(0.0, q, p, params, N)

    # the four-argument methods GeometricEquations requires recover N = length(q)
    @test lagrangian(0.0, q, v, params) == lagrangian(0.0, q, v, params, N)
    @test hamiltonian(0.0, q, p, params) == hamiltonian(0.0, q, p, params, N)

    # `∇V!` is where the hand-written formulation could go wrong: the Dirichlet values q₀ = q_{N+1} = 0
    # enter only through the two end rows of the stencil, which nothing else checks. Differentiate
    # `potential` itself, which shares no code with `∇V!`. N = 1 is included because there the two
    # end rows coincide.
    function reference_∇V!(dV, q, params, N)
        h = 1e-6
        for i in eachindex(q)
            qp = copy(q)
            qp[i] += h
            qm = copy(q)
            qm[i] -= h
            dV[i] = (LW.potential(qp, params, N) - LW.potential(qm, params, N)) / 2h
        end
        nothing
    end

    for M in (1, 2, 3, 5, 12)
        # perturbed off the initial condition, so that no stencil entry is accidentally zero
        qM = LW.initial_condition(LW.μ̃, M).q .+ 0.1 .* collect(1:M)
        kernel = zeros(M)
        reference = zeros(M)
        LW.∇V!(kernel, qM, params, M)
        reference_∇V!(reference, qM, params, M)
        @test kernel ≈ reference rtol = 1e-6
        @test any(!iszero, kernel)

        # α scales the gradient; α = -1 is what the force function uses
        negated = zeros(M)
        LW.∇V!(negated, qM, params, M, -1)
        @test negated ≈ -kernel
    end

    # The Lagrangian is regular (ϑ = q̇, M = I), so the two-form is the 2n×2n canonical [0 -I; I 0].
    Ω = ones(2n, 2n)
    LW.ω!(Ω, 0.0, q, v, params)
    @test Ω ≈ [zeros(n, n) -I(n); I(n) zeros(n, n)]
    @test Ω ≈ -transpose(Ω)
    @test count(!iszero, Ω) == 2n

    # ...and it agrees with what EulerLagrange generates.
    generated_ω = zeros(2n, 2n)
    lagrangian_functions.ω(generated_ω, 0.0, q, v, params)
    @test generated_ω ≈ Ω
end

@testset "$(rpad("Linear Wave: integration agreement",80))" begin
    N = 5
    timespan = (0.0, 0.5)
    timestep = 0.005

    hode = hodeproblem(N; timespan = timespan, timestep = timestep)
    lode = lodeproblem(N; timespan = timespan, timestep = timestep)

    hode_sol = integrate_quietly(hode, ImplicitMidpoint())
    lode_sol = integrate_quietly(lode, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-11

    # The hand-written vector fields exist precisely to replace the generated ones, so they have to
    # keep agreeing with them.
    hode_sym = hodeproblem(N; timespan = timespan, timestep = timestep, symbolic = true)
    lode_sym = lodeproblem(N; timespan = timespan, timestep = timestep, symbolic = true)

    hode_sym_sol = integrate_quietly(hode_sym, ImplicitMidpoint())
    lode_sym_sol = integrate_quietly(lode_sym, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, hode_sym_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, hode_sym_sol.p) < 1E-11
    @test relative_maximum_error(lode_sol.q, lode_sym_sol.q) < 1E-12
    @test relative_maximum_error(lode_sol.p, lode_sym_sol.p) < 1E-11

    # The two-argument form infers N from the state, in both branches.
    q₀ = LW.initial_condition(LW.μ̃, N).q
    p₀ = LW.initial_condition(LW.μ̃, N).p
    @test hodeproblem(q₀, p₀) isa HODEProblem
    @test lodeproblem(q₀, p₀) isa LODEProblem
    @test lodeproblem(q₀, p₀; symbolic = true) isa LODEProblem

    # The hand-written fields read the size off the state while the symbolic ones bake it in, so an
    # `N` that disagrees with the initial condition has to be rejected rather than silently meaning
    # two different systems depending on `symbolic`.
    @test_throws AssertionError hodeproblem(N + 1, q₀, p₀)
    @test_throws AssertionError lodeproblem(N + 1, q₀, p₀)

    # ...and every *sample* of an ensemble is checked, not just the first. A ragged ensemble is
    # accepted by neither branch on its own merits: the hand-written fields would integrate sample 2
    # as a lattice of its own size, and the generated ones — which have N baked in — would write only
    # the first N components of its force and leave the rest at whatever the buffer held. Both run
    # to completion without complaint, which is what makes the check the only thing standing here.
    q₁ = [q₀; 0.1]      # one component too many
    p₁ = [p₀; 0.0]

    @test_throws AssertionError hodeensemble(N, [q₀, q₁], [p₀, p₀])
    @test_throws AssertionError lodeensemble(N, [q₀, q₁], [p₀, p₀])
    @test_throws AssertionError hodeensemble(N, [q₀, q₀], [p₀, p₁])
    @test_throws AssertionError lodeensemble(N, [q₀, q₀], [p₀, p₁])

    # The two-argument forms infer N from the first sample, so a ragged tail has to be caught there
    # too rather than defining the lattice size out from under the rest.
    @test_throws AssertionError hodeensemble([q₀, q₁], [p₀, p₀])
    @test_throws AssertionError lodeensemble([q₀, q₁], [p₀, p₀])

    # ...while a well-formed ensemble is of course still accepted, in both branches.
    @test hodeensemble(N, [q₀, q₀ .+ 0.1], [p₀, p₀]) isa HODEEnsemble
    @test lodeensemble(N, [q₀, q₀ .+ 0.1], [p₀, p₀]; symbolic = true) isa LODEEnsemble
end

@testset "$(rpad("Linear Wave: Dirichlet full-order model",80))" begin
    # The FOM of Buchfink, Glas and Haasdonk: N interior points, q₀ = q_{N+1} = 0 and
    # H = ½ pᵀp - ½ μ² qᵀ D q with the central-difference matrix D = tridiag(1, -2, 1) / Δξ².
    N = 7
    params = LW.default_parameters()
    μ = params.μ
    Δξ = 1 / (N + 1)
    D = SymTridiagonal(fill(-2.0, N), fill(1.0, N - 1)) / Δξ^2

    q = collect(range(-0.3, 0.8, length = N)) .^ 2
    p = collect(range(0.4, -0.2, length = N))
    @test hamiltonian(0.0, q, p, params) ≈ p' * p / 2 - μ^2 * q' * D * q / 2
    f = zeros(N)
    LW.linear_wave_f(f, 0.0, q, p, params)
    @test f ≈ μ^2 * D * q

    # the grid is interior: N points at spacing 1/(N+1), none of them on the boundary
    Ω = LW.interior_domain(N)
    @test length(Ω) == N
    @test Ω[begin] ≈ -1 / 2 + Δξ
    @test Ω[end] ≈ 1 / 2 - Δξ
    @test diff(Ω) ≈ fill(Δξ, N - 1)
    @test LW.initial_condition(μ, N).q ≈ LW.u₀(collect(Ω), μ)

    # Integrated from `initial_condition`, the discrete solution approximates the travelling wave
    # u₀(ξ - μt) of the continuous problem while the pulse is away from the boundary, and the
    # implicit midpoint rule conserves the quadratic Hamiltonian to round-off.
    M = 256
    timespan = (0.0, 0.25)
    timestep = 1 / 400
    sol = integrate_quietly(hodeproblem(M; timespan = timespan, timestep = timestep),
        ImplicitMidpoint())
    Ω = collect(LW.interior_domain(M))
    exact = LW.u₀(Ω .- μ * timespan[2], μ)
    @test maximum(abs, sol.q[end] - exact) < 2e-2
    H = [hamiltonian(0.0, sol.q[n], sol.p[n], params) for n in eachindex(sol.t)]
    @test maximum(abs, H .- H[begin]) / H[begin] < 1e-10

    # An ensemble with one μ per sample and an initial condition of its own for each μ, which is how
    # a training set over a parameter range is built.
    μs = [5 / 12, 1 / 2, 2 / 3]
    ics = [LW.initial_condition(μᵢ, N) for μᵢ in μs]
    ensemble = hodeensemble(N, [ic.q for ic in ics], [ic.p for ic in ics];
        parameters = [(μ = μᵢ,) for μᵢ in μs], timespan = (0.0, 0.1), timestep = 0.01)
    esol = integrate_quietly(ensemble, ImplicitMidpoint())
    @test length(esol.s) == length(μs)
    for (i, ic) in pairs(ics)
        @test esol.s[i].q[0] ≈ ic.q
        @test esol.s[i].p[0] ≈ ic.p
    end
end
