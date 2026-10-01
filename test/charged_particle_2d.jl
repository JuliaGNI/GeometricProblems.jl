using Test
using GeometricIntegrators: Gauss, integrate, relative_maximum_error

import GeometricProblems.ChargedParticle2d as cp
import GeometricProblems.MasslessChargedParticle as mcp

include("helpers/integrate_quietly.jl")

@testset "$(rpad("Charged Particle 2D",80))" begin
    @test_nowarn cp.hodeproblem()
    @test_nowarn cp.lodeproblem()

    params = cp.default_parameters()
    hsol = integrate(cp.hodeproblem(; timespan = (0.0, 1.0)), Gauss(8))
    lsol = integrate(cp.lodeproblem(; timespan = (0.0, 1.0)), Gauss(8))
    @test relative_maximum_error(hsol.q, lsol.q) < 1E-12
    @test relative_maximum_error(hsol.p, lsol.p) < 1E-12

    # the static fields conserve the energy
    H = [cp.hamiltonian(hsol.t[n], hsol.q[n], hsol.p[n], params) for n in eachindex(hsol.t)]
    @test maximum(abs.(H .- H[begin])) < 1E-13 * abs(H[begin])

    # μ = m |v - v_E|² / 2B with the velocity v = (p - A) / m and the E × B drift v_E
    v = zeros(2)
    cp.v̄(v, 0.0, hsol.q[end], hsol.p[end], params)
    w = v .- cp.drift_velocity(hsol.q[end], params)
    @test cp.magnetic_moment(0.0, hsol.q[end], hsol.p[end], params) ≈
          params.m * (w[1]^2 + w[2]^2) / (2 * mcp.B(hsol.q[end], params))
end

# For small masses the residual of the implicit stages has a round-off floor of about eps |A| / m,
# from the cancellation in v = (p - A) / m, which lies above the solver's default absolute
# tolerance; `f_abstol` is set just above that floor.
const f_abstol = 1E-12

@testset "$(rpad("Charged Particle 2D (guiding-centre limit)",80))" begin
    # Started with the E × B drift velocity, the particle follows the massless guiding-centre
    # trajectory up to O(m): reducing m fivefold reduces the maximal deviation fivefold. The time
    # step resolves the gyration, h B(q₀) / m = 1/2.
    T = 1.0
    params = cp.default_parameters()
    ref = integrate(mcp.odeproblem(cp.q₀; timespan = (0.0, T), timestep = 1E-3,
            parameters = (A₀ = params.A₀, E₀ = params.E₀)), Gauss(4))

    Δq = map((1E-2, 2E-3)) do m
        params = merge(cp.default_parameters(), (m = m,))
        sol = integrate_quietly(
            cp.hodeproblem(; timespan = (0.0, T), timestep = m / 10, parameters = params),
            Gauss(2); f_abstol = f_abstol)
        k = round(Int, 1E-3 / (m / 10))
        maximum(n -> maximum(abs.(sol.q[k * n] .- ref.q[n])), eachindex(ref.t))
    end
    @test 4.5 < Δq[1] / Δq[2] < 5.5
    @test Δq[1] < 3E-2 * 1E-2
end

@testset "$(rpad("Charged Particle 2D (magnetic moment)",80))" begin
    # Started with a gyration velocity on top of the drift, the magnetic moment deviates from its
    # initial value by O(m), while B(q) varies by about 25 % along the orbit.
    Δμ = map((1E-2, 2E-3)) do m
        params = merge(cp.default_parameters(), (m = m,))
        p₀ = cp.initial_momentum(
            cp.q₀, cp.drift_velocity(cp.q₀, params) .+ [0.5, 0.0], params)
        sol = integrate_quietly(
            cp.hodeproblem(cp.q₀, p₀; timespan = (0.0, 1.0), timestep = m / 20,
                parameters = params),
            Gauss(2); f_abstol = f_abstol)
        μ = [cp.magnetic_moment(sol.t[n], sol.q[n], sol.p[n], params) for n in eachindex(sol.t)]
        maximum(abs.(μ .- μ[begin])) / μ[begin]
    end
    @test 4.5 < Δμ[1] / Δμ[2] < 5.5
    @test Δμ[1] < 1E-2
end
