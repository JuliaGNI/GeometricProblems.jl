using GeometricIntegrators: Gauss, integrate, relative_maximum_error
using GeometricProblems.ThreeBody
using Test

# Regression tests for the three-body problem:
#  * the potential must contain all three pairwise gravitational interactions (previously the
#    body-1↔body-3 term was missing);
#  * `lodeproblem` must construct and integrate (it used to reference undefined names and threw).
# The default initial condition is the figure-eight choreography over one of its periods, which is
# what makes the assertions below possible: it is a closed, collision-free orbit, so the trajectory
# has to return to where it started, and the Hamiltonian and Lagrangian forms have to agree to
# roundoff. The `initial_conditions` grid cannot support either check — every member of it ends in a
# collision (see `ThreeBody.sympnets_initial_condition`).

@testset "$(rpad("Three-Body Problem",80))" begin
    params = ThreeBody.default_parameters()

    G  = params.G
    m₁ = params.m₁
    m₂ = params.m₂
    m₃ = params.m₃

    # positions: body 1 = (0,0), body 2 = (3,0), body 3 = (0,4)
    q = [0.0, 0.0, 3.0, 0.0, 0.0, 4.0]
    r₁₂ = 3.0
    r₂₃ = 5.0
    r₁₃ = 4.0
    V_expected = -G * m₁ * m₂ / r₁₂ - G * m₂ * m₃ / r₂₃ - G * m₁ * m₃ / r₁₃

    # all three pairwise interactions present (B3 regression)
    @test ThreeBody.V(q, params) ≈ V_expected

    # both problem types construct and integrate (B4 regression: `lodeproblem` previously threw an
    # UndefVarError), and agree to roundoff over the default window.
    hsol = integrate(hodeproblem(), Gauss(2))
    lsol = integrate(lodeproblem(), Gauss(2))

    @test relative_maximum_error(hsol.q, lsol.q) < 1E-14
    @test relative_maximum_error(hsol.p, lsol.p) < 1E-14

    # The choreography is periodic, so after one period both formulations must return to the initial
    # condition. The residual is set by the eight-digit precision of the tabulated constants, not by
    # the integrator: Gauss(2) closes to 1.3E-7 in q and 2.0E-7 in p, and conserves energy to 3E-15.
    # A sign or factor error anywhere in the vector fields would break the closure outright.
    @test maximum(abs, hsol.q[end] .- ThreeBody.figure_eight.q) < 1E-6
    @test maximum(abs, hsol.p[end] .- ThreeBody.figure_eight.p) < 1E-6
    @test maximum(abs, lsol.q[end] .- ThreeBody.figure_eight.q) < 1E-6

    H₀ = ThreeBody.hamiltonian(hsol.t[0], hsol.q[0], hsol.p[0], params)
    H₁ = ThreeBody.hamiltonian(hsol.t[end], hsol.q[end], hsol.p[end], params)
    @test abs(H₁ - H₀) < 1E-12

    # Hamiltonian and Lagrangian formulations must agree for non-unit masses too, so that a
    # mass-scaling error in the Lagrangian kinetic term (∝ q̇²/m instead of m·q̇²) shows up as an O(1)
    # discrepancy. Unequal masses destroy the choreography — the orbit then decays into a close
    # encounter within a period — so this check keeps an explicit short window instead of the default.
    mparams = (m₁ = 1.0, m₂ = 2.0, m₃ = 3.0, G = 1.0)
    h = integrate(hodeproblem(; timespan = (0.0, 0.05), timestep = 0.01, parameters = mparams), Gauss(2))
    l = integrate(lodeproblem(; timespan = (0.0, 0.05), timestep = 0.01, parameters = mparams), Gauss(2))
    @test relative_maximum_error(h.q, l.q) < 1E-6
end
