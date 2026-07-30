using GeometricEquations: HODEEnsemble, HODEProblem, LODEEnsemble, LODEProblem, functions, initialguess
using GeometricIntegrators: Gauss, ImplicitMidpoint, integrate, relative_maximum_error
using GeometricProblems.TodaLattice
using LinearAlgebra
using Test

const TL = TodaLattice


# Parameters and initial conditions

N = 20
μ = TodaLattice.μ
q₀ = TodaLattice.compute_initial_q(μ, N)
p₀ = zero(q₀)
params = TodaLattice.default_parameters()


# Ensemble initial conditions and parameters

q₀_vec = [q₀ .+ α for α in -0.2 : 0.2 : +0.2]
p₀_vec = [p₀ .+ α for α in -0.2 : 0.2 : +0.2]

param_vec = [NamedTuple{keys(params)}(values(params) .+ β) for β in -0.1 : 0.1 : +0.1]


# HODE and LODE problems

hode_prb = hodeproblem(q₀, p₀)
lode_prb = lodeproblem(q₀, p₀)

href_sol = integrate(hode_prb, Gauss(1))
lref_sol = integrate(lode_prb, Gauss(1))

@test relative_maximum_error(href_sol.q, lref_sol.q) < 2E-14


# Ensemble problems with different initial conditions
hode_ens = hodeensemble(N, q₀_vec, p₀_vec)
lode_ens = lodeensemble(N, q₀_vec, p₀_vec)

hode_sol = integrate(hode_ens, Gauss(1))
lode_sol = integrate(lode_ens, Gauss(1))

@test relative_maximum_error(hode_sol[2].q, href_sol.q) < 2E-14
@test relative_maximum_error(lode_sol[2].q, lref_sol.q) < 2E-14

for (hsol,lsol) in zip(hode_sol,lode_sol)
    @test relative_maximum_error(hsol.q, lsol.q) < 2E-14
end


# Ensemble problems with different parameters
hode_ens = hodeensemble(N; parameters = param_vec)
lode_ens = lodeensemble(N; parameters = param_vec)

hode_sol = integrate(hode_ens, Gauss(1))
lode_sol = integrate(lode_ens, Gauss(1))

@test relative_maximum_error(hode_sol[2].q, href_sol.q) < 2E-14
@test relative_maximum_error(lode_sol[2].q, lref_sol.q) < 2E-14

for (hsol,lsol) in zip(hode_sol,lode_sol)
    @test relative_maximum_error(hsol.q, lsol.q) < 2E-14
end


# Ensemble problems with different initial conditions and parameters
hode_ens = hodeensemble(q₀_vec, p₀_vec; parameters = param_vec)
lode_ens = lodeensemble(q₀_vec, p₀_vec; parameters = param_vec)

hode_sol = integrate(hode_ens, Gauss(1))
lode_sol = integrate(lode_ens, Gauss(1))

@test relative_maximum_error(hode_sol[2].q, href_sol.q) < 2E-14
@test relative_maximum_error(lode_sol[2].q, lref_sol.q) < 2E-14

for (hsol,lsol) in zip(hode_sol,lode_sol)
    @test relative_maximum_error(hsol.q, lsol.q) < 2E-14
end


# The two default-size constructions below used to generate the equations of motion symbolically for
# Ñ = 200; `lodeproblem` alone cost minutes there, most of it spent on a dense 400×400 two-form that
# no integrator ever evaluates. They are now instantaneous, and the identity assertions are what keep
# them that way: if the default ever routes through EulerLagrange again, these fail rather than
# merely getting slow.

@testset "$(rpad("Toda Lattice: default constructors are hand-written",80))" begin
    hode = @test_nowarn hodeproblem()
    lode = @test_nowarn lodeproblem()

    @test hode isa HODEProblem
    @test lode isa LODEProblem

    @test functions(hode).v === TL.toda_lattice_v
    @test functions(hode).f === TL.toda_lattice_f
    @test functions(hode).h === TL.hamiltonian

    @test functions(lode).ϑ === TL.toda_lattice_ϑ
    @test functions(lode).f === TL.toda_lattice_f
    @test functions(lode).g === TL.toda_lattice_g
    @test functions(lode).ω === TL.ω!
    @test functions(lode).l === TL.lagrangian

    # was `functions(hamiltonian_system(N, …)).v`, for which an entire second `HamiltonianSystem`
    # had to be built although the function it yields is just `v .= p`
    @test initialguess(lode).v === TL.v̄
end

@testset "$(rpad("Toda Lattice: hand-written matches generated",80))" begin
    # Small on purpose: `lagrangian_system` grows as N^2.4 and the generated `ω` as N^2.
    M = 5

    q = TL.compute_initial_q(μ, M) .+ 0.05 .* collect(1 : M)
    p = collect(range(0.1, 0.4, length = M))
    v = copy(p)                 # the mass matrix is the identity, so q̇ = p

    hamiltonian_functions = functions(TL.hamiltonian_system(M, params))
    lagrangian_functions = functions(TL.lagrangian_system(M, params))

    generated = zeros(M)
    handwritten = zeros(M)

    # ∂H/∂p = p
    hamiltonian_functions.v(generated, 0.0, q, p, params)
    TL.toda_lattice_v(handwritten, 0.0, q, p, params)
    @test generated ≈ handwritten
    @test handwritten == p

    # -∂H/∂q, and the identical ∂L/∂q
    hamiltonian_functions.f(generated, 0.0, q, p, params)
    TL.toda_lattice_f(handwritten, 0.0, q, p, params)
    @test generated ≈ handwritten
    @test any(!iszero, handwritten)

    lagrangian_functions.f(generated, 0.0, q, v, params)
    TL.toda_lattice_f(handwritten, 0.0, q, v, params)
    @test generated ≈ handwritten

    # ϑ = ∂L/∂q̇ = q̇
    lagrangian_functions.ϑ(generated, 0.0, q, v, params)
    TL.toda_lattice_ϑ(handwritten, 0.0, q, v, params)
    @test generated ≈ handwritten
    @test any(!iszero, handwritten)

    # g = λ, in both arities
    λ = collect(range(0.1, 0.7, length = M))
    lagrangian_functions.g(generated, 0.0, q, v, λ, params)
    TL.toda_lattice_g(handwritten, 0.0, q, v, λ, params)
    @test generated ≈ handwritten ≈ λ
    TL.toda_lattice_g(handwritten, 0.0, q, λ, params)
    @test handwritten ≈ λ

    @test lagrangian_functions.L(0.0, q, v, params) ≈ lagrangian(0.0, q, v, params, M)
    @test hamiltonian_functions.H(0.0, q, p, params) ≈ hamiltonian(0.0, q, p, params, M)

    # the four-argument methods GeometricEquations requires recover N = length(q)
    @test lagrangian(0.0, q, v, params) == lagrangian(0.0, q, v, params, M)
    @test hamiltonian(0.0, q, p, params) == hamiltonian(0.0, q, p, params, M)

    # `∇V!` is where the hand-written formulation could go wrong: the lattice is periodic, so the
    # gradient α (Eⱼ - Eⱼ₋₁) wraps around at j = 1, and that wrap-around is checked by nothing else.
    # Differentiate `potential` itself, which shares no code with `∇V!`. N = 1 and N = 2 are included
    # because they are the cases where the wrap-around collapses: at N = 1 the neighbour of the only
    # point is itself, and at N = 2 the forward and backward neighbours coincide.
    function reference_∇V!(dV, q, params, N)
        h = 1e-6
        for i in eachindex(q)
            qp = copy(q); qp[i] += h
            qm = copy(q); qm[i] -= h
            dV[i] = (TL.potential(qp, params, N) - TL.potential(qm, params, N)) / 2h
        end
        nothing
    end

    for K in (1, 2, 3, 5, 12, 64)
        # perturbed off the initial condition, so that no entry is accidentally zero
        qK = TL.compute_initial_q(μ, K) .+ 0.1 .* collect(1:K)
        kernel = zeros(K)
        reference = zeros(K)
        TL.∇V!(kernel, qK, params, K)
        reference_∇V!(reference, qK, params, K)
        @test kernel ≈ reference rtol = 1e-6 atol = 1e-8
        # At K = 1 the only neighbour of the only point is itself, so V is constant and ∇V vanishes
        # identically — the one size at which an all-zero gradient is the right answer.
        K == 1 ? (@test all(iszero, kernel)) : (@test any(!iszero, kernel))

        # `scale` scales the gradient; scale = -1 is what the force function uses
        negated = zeros(K)
        TL.∇V!(negated, qK, params, K, -1)
        @test negated ≈ -kernel

        # Every entry of the output is written once and none is read back, and the wrap-around term
        # is evaluated before the first write, so the kernel is correct even when the output aliases
        # the state. Nothing in the package calls it that way; the assertion is here because the
        # single-pass formulation is what makes it true, and a return to a scratch-space version
        # would silently break it.
        aliased = copy(qK)
        TL.∇V!(aliased, aliased, params, K)
        @test aliased ≈ kernel
    end

    # The Lagrangian is regular (ϑ = q̇, M = I), so the two-form is the 2N×2N canonical [0 -I; I 0].
    Ω = ones(2M, 2M)
    TL.ω!(Ω, 0.0, q, v, params)
    @test Ω ≈ [zeros(M, M) -I(M); I(M) zeros(M, M)]
    @test Ω ≈ -transpose(Ω)
    @test count(!iszero, Ω) == 2M

    # ...and it agrees with what EulerLagrange generates.
    generated_ω = zeros(2M, 2M)
    lagrangian_functions.ω(generated_ω, 0.0, q, v, params)
    @test generated_ω ≈ Ω
end

@testset "$(rpad("Toda Lattice: integration agreement",80))" begin
    M = 5
    timespan = (0.0, 0.5)
    timestep = 0.005

    qM = TL.compute_initial_q(μ, M)
    pM = zero(qM)

    hode = hodeproblem(M, qM, pM; timespan = timespan, timestep = timestep)
    lode = lodeproblem(M, qM, pM; timespan = timespan, timestep = timestep)

    hode_sol = integrate(hode, ImplicitMidpoint())
    lode_sol = integrate(lode, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-11

    # The hand-written vector fields exist precisely to replace the generated ones, so they have to
    # keep agreeing with them.
    hode_sym = hodeproblem(M, qM, pM; timespan = timespan, timestep = timestep, symbolic = true)
    lode_sym = lodeproblem(M, qM, pM; timespan = timespan, timestep = timestep, symbolic = true)

    hode_sym_sol = integrate(hode_sym, ImplicitMidpoint())
    lode_sym_sol = integrate(lode_sym, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, hode_sym_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, hode_sym_sol.p) < 1E-11
    @test relative_maximum_error(lode_sol.q, lode_sym_sol.q) < 1E-12
    @test relative_maximum_error(lode_sol.p, lode_sym_sol.p) < 1E-11

    # The `symbolic = true` LODE keeps the hand-written `v̄` too — that is the whole point of no
    # longer building a second `HamiltonianSystem` for it.
    @test initialguess(lode_sym).v === TL.v̄

    # The two-argument form infers N from the state, in both branches.
    @test hodeproblem(qM, pM) isa HODEProblem
    @test lodeproblem(qM, pM) isa LODEProblem
    @test lodeproblem(qM, pM; symbolic = true) isa LODEProblem

    # ...and it does so through `_length`, so a vector of samples gives the size of one sample rather
    # than the number of samples.
    @test hodeensemble([qM, qM .+ 0.1], [pM, pM]) isa HODEEnsemble
    @test lodeensemble([qM, qM .+ 0.1], [pM, pM]) isa LODEEnsemble

    # The hand-written fields read the size off the state while the symbolic ones bake it in, so an
    # `N` that disagrees with the initial condition has to be rejected rather than silently meaning
    # two different systems depending on `symbolic`.
    @test_throws AssertionError hodeproblem(M + 1, qM, pM)
    @test_throws AssertionError lodeproblem(M + 1, qM, pM)
end

# The ensembles above are hand-written now, so the symbolic branch of `hodeensemble`/`lodeensemble` —
# which destructures `functions(...)`, a different code path from the problem constructors — is only
# covered if it is asked for explicitly. `test/eulerlagrange_ensembles_tests.jl` exists to exercise
# the EulerLagrange ensemble plumbing and names this file as its `TodaLattice` counterpart, so the
# coverage has to stay here. Small lattice, as there: construction is all this checks, and
# `lagrangian_system` grows as N^2.4.
@testset "$(rpad("Toda Lattice: symbolic ensembles (construction)",80))" begin
    M = 5
    qM = TL.compute_initial_q(μ, M)
    pM = zero(qM)
    pv = [params, NamedTuple{keys(params)}(values(params) .+ 0.1)]

    @test hodeensemble(M; parameters = pv, symbolic = true) isa HODEEnsemble
    @test lodeensemble(M; parameters = pv, symbolic = true) isa LODEEnsemble
    @test hodeensemble([qM, qM .+ 0.01], [pM, pM]; symbolic = true) isa HODEEnsemble
    @test lodeensemble([qM, qM .+ 0.01], [pM, pM]; symbolic = true) isa LODEEnsemble
end
