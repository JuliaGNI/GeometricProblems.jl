using GeometricIntegrators: Gauss, integrate, relative_maximum_error
using GeometricProblems.ThreeBody
using Test

# Regression tests for the three-body problem:
#  * the potential must contain all three pairwise gravitational interactions (previously the
#    body-1↔body-3 term was missing);
#  * `lodeproblem` must construct and integrate (it used to reference undefined names and threw).
# A direct HODE↔LODE trajectory comparison is not used here: the dynamics from the default
# initial conditions is sensitive/chaotic, so two distinct order-4 integrators diverge quickly.

@testset "$(rpad("Three-Body Problem",80))" begin
    params = ThreeBody.default_parameters

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

    # both problem types construct and integrate, returning finite solutions (B4 regression:
    # `lodeproblem` previously threw an UndefVarError). A short window keeps the stiff implicit
    # solves well-behaved; the coarse default step can trigger line-search warnings.
    hsol = integrate(hodeproblem(; timespan = (0.0, 1.0), timestep = 0.05), Gauss(2))
    lsol = integrate(lodeproblem(; timespan = (0.0, 1.0), timestep = 0.05), Gauss(2))

    @test all(isfinite, hsol.q[end])
    @test all(isfinite, lsol.q[end])

    # Hamiltonian and Lagrangian formulations must agree for non-unit masses too. Over a very short
    # window the two order-4 integrators differ only negligibly, so a mass-scaling error in the
    # Lagrangian kinetic term (∝ q̇²/m instead of m·q̇²) would show up as an O(1) discrepancy.
    mparams = (m₁ = 1.0, m₂ = 2.0, m₃ = 3.0, G = 1.0)
    h = integrate(hodeproblem(; timespan = (0.0, 0.05), timestep = 0.01, parameters = mparams), Gauss(2))
    l = integrate(lodeproblem(; timespan = (0.0, 0.05), timestep = 0.01, parameters = mparams), Gauss(2))
    @test relative_maximum_error(h.q, l.q) < 1E-6
end
