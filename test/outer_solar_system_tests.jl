using Test
using GeometricIntegrators
using GeometricProblems.OuterSolarSystem
using GeometricSolutions

include("integrate_quietly.jl")

const OSS = OuterSolarSystem

@testset "$(rpad("Outer Solar System",80))" begin

    # The hand-written vector fields are the default; `symbolic = true` routes through
    # EulerLagrange. Both are exercised here, since the hand-written functions exist precisely to
    # replace the generated ones and so have to keep agreeing with them.

    hode = @test_nowarn hodeproblem(timespan=(0.0, 10.0))
    lode = @test_nowarn lodeproblem(timespan=(0.0, 10.0))

    hode_sol = integrate_quietly(hode, ImplicitMidpoint())
    lode_sol = integrate_quietly(lode, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-11

    hode_sym = @test_nowarn hodeproblem(timespan=(0.0, 10.0), symbolic=true)
    lode_sym = @test_nowarn lodeproblem(timespan=(0.0, 10.0), symbolic=true)

    hode_sym_sol = integrate_quietly(hode_sym, ImplicitMidpoint())
    lode_sym_sol = integrate_quietly(lode_sym, ImplicitMidpoint())

    @test relative_maximum_error(hode_sol.q, hode_sym_sol.q) < 1E-12
    @test relative_maximum_error(hode_sol.p, hode_sym_sol.p) < 1E-11
    @test relative_maximum_error(lode_sol.q, lode_sym_sol.q) < 1E-12
    @test relative_maximum_error(lode_sol.p, lode_sym_sol.p) < 1E-11

end

@testset "$(rpad("Outer Solar System: hand-written matches generated",80))" begin

    params = default_parameters()
    D = OSS.N_BODIES * OSS.N_DIM

    q = OSS.q₀
    p = OSS.p₀
    v = [p[OSS.N_DIM*i-2+k] / params.m[i] for i in 1:OSS.N_BODIES for k in 0:2]

    hamiltonian_functions = functions(OSS.hamiltonian_system(params))
    lagrangian_functions = functions(OSS.lagrangian_system(params))

    generated = zeros(D)
    handwritten = zeros(D)

    # ∂H/∂p
    hamiltonian_functions.v(generated, 0.0, q, p, params)
    OSS.outer_solar_system_v(handwritten, 0.0, q, p, params)
    @test generated ≈ handwritten

    # -∂H/∂q, and the identical ∂L/∂q
    hamiltonian_functions.f(generated, 0.0, q, p, params)
    OSS.outer_solar_system_f(handwritten, 0.0, q, p, params)
    @test generated ≈ handwritten
    @test any(!iszero, handwritten)

    lagrangian_functions.f(generated, 0.0, q, v, params)
    OSS.outer_solar_system_f(handwritten, 0.0, q, v, params)
    @test generated ≈ handwritten

    # ϑ = ∂L/∂q̇
    lagrangian_functions.ϑ(generated, 0.0, q, v, params)
    OSS.outer_solar_system_ϑ(handwritten, 0.0, q, v, params)
    @test generated ≈ handwritten
    @test any(!iszero, handwritten)

    @test lagrangian_functions.L(0.0, q, v, params) ≈ lagrangian(0.0, q, v, params)
    @test hamiltonian_functions.H(0.0, q, p, params) ≈ hamiltonian(0.0, q, p, params)

    # `∇V!` accumulates each pair into both bodies at once, so check it against a direct
    # body-by-body double sum that shares no code with it.
    function reference_∇V!(dV, q, params)
        fill!(dV, 0)
        for i in 1:OSS.N_BODIES, j in 1:OSS.N_BODIES
            i == j && continue
            Δ = [q[OSS.N_DIM*i-2+k] - q[OSS.N_DIM*j-2+k] for k in 0:2]
            r = sqrt(sum(abs2, Δ))
            for k in 0:2
                dV[OSS.N_DIM*i-2+k] += params.G * params.m[i] * params.m[j] * Δ[k+1] / r^3
            end
        end
        nothing
    end

    reference = zeros(D)
    OSS.∇V!(handwritten, q, params)
    reference_∇V!(reference, q, params)
    @test handwritten ≈ reference

    # The Lagrangian is regular, so the two-form is the 2n×2n canonical [0 -M; M 0].
    M = zeros(D, D)
    for i in 1:OSS.N_BODIES, k in 0:OSS.N_DIM-1
        M[OSS.N_DIM*i-2+k, OSS.N_DIM*i-2+k] = params.m[i]
    end

    Ω = ones(2D, 2D)
    OSS.ω!(Ω, 0.0, q, v, params)
    @test Ω ≈ [zeros(D, D) -M; M zeros(D, D)]
    @test Ω ≈ -transpose(Ω)
    @test any(!iszero, Ω)

    # ...and it agrees with what EulerLagrange generates.
    generated_ω = zeros(2D, 2D)
    lagrangian_functions.ω(generated_ω, 0.0, q, v, params)
    @test generated_ω ≈ Ω

end
