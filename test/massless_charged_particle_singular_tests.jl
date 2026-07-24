using Test
using ForwardDiff
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricSolutions

import GeometricProblems.MasslessChargedParticle as mcp
import GeometricProblems.MasslessChargedParticleSingular as mcps

# The singular problem uses a gauge-transformed, one-component vector potential (A₂ = 0). It
# describes the same physics as the standard massless charged particle: same magnetic field B,
# same electric field, hence the same velocity field. Only the one-form ϑ = A differs, which
# yields a different (degenerate) variational integrator. As for the standard problem, the ODE
# vector field must satisfy the Euler–Lagrange relation Ω v = -∇H, where
# Ω_ij = ∂ϑ_i/∂q_j - ∂ϑ_j/∂q_i is the (anti-symmetric) symplectic matrix.

@testset "$(rpad("Massless Charged Particle (singular)",80))" begin
    params = mcps.default_parameters()

    for q in ([1.0, 1.0], [0.5, -0.7], [1.3, 0.2], [-0.9, 0.4])
        v = zeros(2)
        mcps.massless_charged_particle_v(v, 0.0, q, params)

        ∇H = [mcps.dHd₁(0.0, q, params), mcps.dHd₂(0.0, q, params)]

        # Euler–Lagrange / IODE consistency of the ODE vector field: Ω v = -∇H
        @test mcps.ω(0.0, q, params) * v ≈ -∇H atol = 1e-12

        # the drift is orthogonal to ∇H, i.e. the Hamiltonian (energy) is conserved
        @test abs(sum(v .* ∇H)) < 1e-12

        # the singular gauge reproduces the same magnetic field and hence the same dynamics
        @test mcps.B(q, params) ≈ mcp.B(q, params) atol = 1e-12

        v_std = zeros(2)
        mcp.massless_charged_particle_v(v_std, 0.0, q, params)
        @test v ≈ v_std atol = 1e-12
    end
end

# The symplectic two-form Ω is assembled from the hand-coded derivatives dϑᵢdxⱼ. Comparing it to
# the automatic derivative of ϑ therefore validates both the assembly and those derivatives. Since
# the gauge only shifts ϑ by an exact one-form, Ω itself must be identical to the standard gauge.
@testset "$(rpad("Massless Charged Particle (singular, two-form and Lagrangian)",80))" begin
    params = mcps.default_parameters()
    t = 0.0

    for q in ([1.0, 1.0], [0.5, -0.7], [1.3, 0.2], [-0.9, 0.4])
        v = zeros(2)
        mcps.massless_charged_particle_v(v, t, q, params)

        # (∇ϑ)ᵢⱼ = ∂ϑᵢ/∂qⱼ, so the convention Ωᵢⱼ = ∂ϑᵢ/∂qⱼ - ∂ϑⱼ/∂qᵢ means Ω = ∇ϑ - (∇ϑ)ᵀ
        J = ForwardDiff.jacobian(x -> mcps.ϑ(t, x, params), q)
        Ω = mcps.ω(t, q, params)
        @test Ω ≈ J - transpose(J) atol = 1e-12
        @test Ω ≈ -transpose(Ω) atol = 1e-12
        @test Ω[1, 2] ≈ -mcps.B(q, params) atol = 1e-12

        # the gauge transformation leaves the two-form invariant
        @test Ω ≈ mcp.ω(t, q, params) atol = 1e-12

        # the in-place callback agrees with the allocating version
        Ω̃ = zeros(2, 2)
        mcps.massless_charged_particle_ω(Ω̃, t, q, v, params)
        @test Ω̃ ≈ Ω

        # L = ϑ⋅v - H, with ∂L/∂v = ϑ (degenerate, linear in v)
        @test mcps.lagrangian(t, q, v, params) ≈ sum(mcps.ϑ(t, q, params) .* v) - mcps.hamiltonian(t, q, params)
        @test ForwardDiff.gradient(u -> mcps.lagrangian(t, q, u, params), v) ≈ mcps.ϑ(t, q, params) atol = 1e-12
    end

    # No integrator in GeometricIntegrators evaluates ω, and `check_methods` skips it, so the
    # positional wiring of LODEProblem(ϑ, f, g, ω, l, …) / LDAEProblem(…, ω, l, …) has to be
    # checked explicitly — an ω/l swap would otherwise go unnoticed.
    q = [1.3, 0.2]
    v = [0.7, -0.4]
    Ω = mcps.ω(t, q, params)

    for prob in (mcps.lodeproblem(), mcps.ldaeproblem())
        equs = functions(prob)
        pars = parameters(prob)

        Ω̃ = zeros(2, 2)
        equs.ω(Ω̃, t, q, v, pars)
        @test Ω̃ ≈ Ω
        @test equs.l(t, q, v, pars) ≈ mcps.lagrangian(t, q, v, params)
        @test invariants(prob).h(t, q, v, pars) ≈ mcps.hamiltonian(t, q, params)
    end
end

# The variational (IDAE/LODE/LDAE) formulations must reproduce the ODE dynamics. `idaeproblem` is a
# variational integrator problem (VSPARK), `idaeproblem_spark` uses the force-split SPARK form, and
# `lodeproblem`/`ldaeproblem` are the Lagrangian forms built from ω and L = ϑ⋅v - H.
@testset "$(rpad("Massless Charged Particle (singular, variational)",80))" begin
    tspan = (0.0, 1.0)
    tstep = 0.05
    ref = integrate(mcps.odeproblem(; timespan = tspan, timestep = tstep), Gauss(8))

    sol = integrate(mcps.idaeproblem(; timespan = tspan, timestep = tstep), TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(mcps.idaeproblem_spark(; timespan = tspan, timestep = tstep), SPARKGLRK(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(mcps.lodeproblem(; timespan = tspan, timestep = tstep), MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    # SLRK is the second-order-accurate integrator for the degenerate LDAE; note that
    # SLRKLobattoIIIAB is singular for this gauge, hence LobattoIIIE.
    sol = integrate(mcps.ldaeproblem(; timespan = tspan, timestep = tstep), SLRKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-5
end
