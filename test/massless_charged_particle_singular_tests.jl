using Test
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

        Ω₁₂ = mcps.dϑ₁dx₂(0.0, q, params) - mcps.dϑ₂dx₁(0.0, q, params)
        ∇H  = [mcps.dHd₁(0.0, q, params), mcps.dHd₂(0.0, q, params)]
        Ωv  = [Ω₁₂ * v[2], -Ω₁₂ * v[1]]   # (anti-symmetric) Ω acting on v

        # Euler–Lagrange / IODE consistency of the ODE vector field: Ω v = -∇H
        @test Ωv ≈ -∇H atol = 1e-12

        # the drift is orthogonal to ∇H, i.e. the Hamiltonian (energy) is conserved
        @test abs(sum(v .* ∇H)) < 1e-12

        # the singular gauge reproduces the same magnetic field and hence the same dynamics
        @test mcps.B(q, params) ≈ mcp.B(q, params) atol = 1e-12

        v_std = zeros(2)
        mcp.massless_charged_particle_v(v_std, 0.0, q, params)
        @test v ≈ v_std atol = 1e-12
    end
end

# The variational (IDAE) formulations must reproduce the ODE dynamics. `idaeproblem` is a
# variational integrator problem (VSPARK), `idaeproblem_spark` uses the force-split SPARK form.
@testset "$(rpad("Massless Charged Particle (singular, variational)",80))" begin
    tspan = (0.0, 1.0)
    tstep = 0.05
    ref = integrate(mcps.odeproblem(; timespan = tspan, timestep = tstep), Gauss(8))

    sol = integrate(mcps.idaeproblem(; timespan = tspan, timestep = tstep), TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(mcps.idaeproblem_spark(; timespan = tspan, timestep = tstep), SPARKGLRK(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8
end
