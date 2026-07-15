using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricSolutions

import GeometricProblems.MasslessChargedParticle as mcp

# The ODE vector field must be consistent with the variational (IODE/IDAE) formulation built from
# the one-form ϑ = A. Along a solution the Euler–Lagrange relation reads Ω v = -∇H, where
# Ω_ij = ∂ϑ_i/∂q_j - ∂ϑ_j/∂q_i is the (anti-symmetric) symplectic matrix. This test would fail for
# the previous, sign-reversed ODE vector field (it gave Ω v = +∇H, i.e. a time-reversed flow).

@testset "$(rpad("Massless Charged Particle",80))" begin
    params = mcp.default_parameters()

    for q in ([1.0, 1.0], [0.5, -0.7], [1.3, 0.2], [-0.9, 0.4])
        v = zeros(2)
        mcp.massless_charged_particle_v(v, 0.0, q, params)

        Ω₁₂ = mcp.dϑ₁dx₂(0.0, q, params) - mcp.dϑ₂dx₁(0.0, q, params)
        ∇H  = [mcp.dHd₁(0.0, q, params), mcp.dHd₂(0.0, q, params)]
        Ωv  = [Ω₁₂ * v[2], -Ω₁₂ * v[1]]   # (anti-symmetric) Ω acting on v

        # Euler–Lagrange / IODE consistency of the ODE vector field: Ω v = -∇H
        @test Ωv ≈ -∇H atol = 1e-12

        # the drift is orthogonal to ∇H, i.e. the Hamiltonian (energy) is conserved
        @test abs(sum(v .* ∇H)) < 1e-12
    end
end

# The variational (IDAE) formulations must reproduce the ODE dynamics. `idaeproblem` is a
# variational integrator problem (VSPARK), `idaeproblem_spark` uses the force-split SPARK form.
@testset "$(rpad("Massless Charged Particle (variational)",80))" begin
    tspan = (0.0, 1.0)
    tstep = 0.05
    ref = integrate(mcp.odeproblem(; timespan = tspan, timestep = tstep), Gauss(8))

    sol = integrate(mcp.idaeproblem(; timespan = tspan, timestep = tstep), TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(mcp.idaeproblem_spark(; timespan = tspan, timestep = tstep), SPARKGLRK(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8
end
