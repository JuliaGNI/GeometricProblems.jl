using Test
using CairoMakie
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricSolutions
using GeometricEquations: invariants, parameters

import GeometricProblems.Diagnostics as diag
import GeometricProblems.HarmonicOscillator as ho
import GeometricProblems.Pendulum as pend
import GeometricProblems.LotkaVolterra2d as lv2
import GeometricProblems.LotkaVolterra3d as lv3
import GeometricProblems.LotkaVolterra4d as lv4
import GeometricProblems.MasslessChargedParticle as mcp

# Loading CairoMakie activates the `Makie` weakdep, and with it all `*Plots`
# extensions. These are smoke tests: each plot function must build and return a
# Makie `Figure` without error.

@testset "$(rpad("Plotting extensions",80))" begin

    @testset "Diagnostics + Lotka-Volterra 2d" begin
        ode  = lv2.odeproblem()
        idae = lv2.idaeproblem()

        sol  = integrate(ode, Gauss(2))
        dsol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))

        # problem-specific plots
        @test lv2.plot_phase_portrait(sol)  isa Figure
        @test lv2.plot_solution(sol, ode)   isa Figure
        @test lv2.plot_traces(sol, ode)     isa Figure

        # generic diagnostics
        _, ΔH = compute_invariant_error(sol.t, sol.q, parameters(ode), invariants(ode)[:h])
        @test diag.plot_energy_error(sol)         isa Figure
        @test diag.plot_energy_error(sol.t, ΔH)   isa Figure
        # `compute_drift` is not available in the pinned GeometricSolutions, so the
        # pointwise error series ΔH stands in as a (time, values) pair here — this
        # only smoke-tests that a Figure is built, not the drift semantics.
        @test diag.plot_energy_drift(sol.t, ΔH)   isa Figure
        @test diag.plot_constraint_error(dsol)    isa Figure
        @test diag.plot_lagrange_multiplier(dsol) isa Figure
    end

    @testset "Lotka-Volterra 3d" begin
        ode = lv3.odeproblem()
        sol = integrate(ode, Gauss(2))
        @test lv3.plot_phase_portrait(sol)          isa Figure
        @test lv3.plot_traces(sol, parameters(ode)) isa Figure
    end

    @testset "Lotka-Volterra 4d" begin
        ode = lv4.odeproblem()
        sol = integrate(ode, Gauss(2))
        @test lv4.plot_phase_portrait(sol)          isa Figure
        @test lv4.plot_traces(sol, parameters(ode)) isa Figure
    end

    @testset "Massless charged particle" begin
        ode = mcp.odeproblem()
        sol = integrate(ode, Gauss(2))
        @test mcp.plot_phase_portrait(sol)  isa Figure
        @test mcp.plot_solution(sol, ode)   isa Figure
        @test mcp.plot_traces(sol, ode)     isa Figure
    end

    @testset "Harmonic oscillator" begin
        sol = integrate(ho.hodeproblem(; timespan = (0.0, 10.0), timestep = 0.1), Gauss(2))
        @test ho.plot_phase_portrait(sol) isa Figure
        @test ho.plot_traces(sol)         isa Figure
        @test ho.plot_hamiltonian()       isa Figure
    end

    @testset "Pendulum" begin
        sol = integrate(pend.hodeproblem([1.0], [0.0]; timespan = (0.0, 10.0), timestep = 0.1), ImplicitMidpoint())
        @test pend.plot_phase_portrait(sol) isa Figure
        @test pend.plot_traces(sol)         isa Figure
        @test pend.plot_hamiltonian()       isa Figure
    end

end
