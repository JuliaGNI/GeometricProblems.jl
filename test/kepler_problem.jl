using Test
using GeometricIntegrators: Gauss, integrate, relative_maximum_error

import GeometricProblems.KeplerProblem as kepler

const params = kepler.default_parameters()

@testset "$(rpad("Kepler Problem",80))" begin
    @test_nowarn kepler.hodeproblem()
    @test_nowarn kepler.lodeproblem()

    q₀, p₀ = kepler.initial_condition(0.9)
    hsol = integrate(kepler.hodeproblem(q₀, p₀; timespan = (0.0, 2π), timestep = 2π / 1000), Gauss(8))
    lsol = integrate(kepler.lodeproblem(q₀, p₀; timespan = (0.0, 2π), timestep = 2π / 1000), Gauss(8))
    @test relative_maximum_error(hsol.q, lsol.q) < 1E-12
    @test relative_maximum_error(hsol.p, lsol.p) < 1E-12

    # one period of the orbit with e = 0.9 agrees with the closed-form solution
    for n in eachindex(hsol.t)
        q, p = kepler.exact_solution(hsol.t[n], 0.0, q₀, p₀, params)
        @test hsol.q[n] ≈ q atol = 1E-10
        @test hsol.p[n] ≈ p atol = 1E-10 * p₀[2]
    end
end

@testset "$(rpad("Kepler Problem (exact solution)",80))" begin
    # Kepler's equation is solved to round-off, also for e → 1.
    for e in (0.01, 0.5, 0.9, 0.99, 0.999), M in range(-π, π; length = 1001)
        E = kepler.eccentric_anomaly(M, e)
        @test E - e * sin(E) - M ≈ 0 atol = 1E-14
    end

    # perihelion data, and a prograde and a retrograde orbit through a generic point
    for (q₀, p₀) in (kepler.initial_condition(0.99), ([0.3, -0.4], [0.9, 1.1]),
        ([0.3, -0.4], [-0.9, -1.1]))
        H₀ = kepler.hamiltonian(0.0, q₀, p₀, params)
        ℓ₀ = kepler.angular_momentum(0.0, q₀, p₀, params)
        R₀ = kepler.runge_lenz_vector(0.0, q₀, p₀, params)
        period = 2π * sqrt((-1 / (2H₀))^3)

        # relative tolerances: for e → 1 the semi-major axis a = -μ/2H is computed from a
        # cancellation of the kinetic and potential energy at the perihelion
        q, p = kepler.exact_solution(0.0, 0.0, q₀, p₀, params)
        @test [q; p] ≈ [q₀; p₀] rtol = 1E-13
        q, p = kepler.exact_solution(period, 0.0, q₀, p₀, params)
        @test [q; p] ≈ [q₀; p₀] rtol = 1E-12

        for t in range(0, 7; length = 50)
            q, p = kepler.exact_solution(t, 0.0, q₀, p₀, params)
            @test kepler.hamiltonian(t, q, p, params) ≈ H₀ atol = 1E-13
            @test kepler.angular_momentum(t, q, p, params) ≈ ℓ₀ atol = 1E-13
            @test kepler.runge_lenz_vector(t, q, p, params) ≈ R₀ atol = 1E-13
        end

        # q̇ = p and ṗ = -μ q / |q|³, checked by central differences
        δ = 1E-5
        q₊, p₊ = kepler.exact_solution(1.234 + δ, 0.0, q₀, p₀, params)
        q₋, p₋ = kepler.exact_solution(1.234 - δ, 0.0, q₀, p₀, params)
        q, p = kepler.exact_solution(1.234, 0.0, q₀, p₀, params)
        @test (q₊ - q₋) / 2δ ≈ p atol = 1E-8
        @test (p₊ - p₋) / 2δ ≈ -params.μ * q / sqrt(q[1]^2 + q[2]^2)^3 atol = 1E-8
    end
end
