# Lotka-Volterra 2d: a survey of integrators
#
# Integrates the same system with a representative selection of the methods GeometricIntegrators
# offers for it, and writes one PDF per method next to this script. Each figure shows the time
# traces x₁(t), x₂(t) together with the relative energy error, which is the quantity these methods
# are built to preserve.
#
# The problem is available in several formulations — ODE, implicit ODE, and implicit DAE — and not
# every method applies to every one of them. The pairings below are the ones exercised in
# `test/lotka_volterra_2d_tests.jl`.

using CairoMakie
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2d

# The module defaults are t ∈ [0, 10] at Δt = 0.01, which is a thousand steps and runs in about a
# second per method. Widen the timespan to see the energy error develop.
const timespan = (0.0, 10.0)
const timestep = 0.01

const runs = (
    ("gauss2",                    odeproblem(;  timespan, timestep), Gauss(2)),
    ("gauss8",                    odeproblem(;  timespan, timestep), Gauss(8)),
    ("vprk-gauss2-midpoint",      iodeproblem(; timespan, timestep), MidpointProjection(VPRKGauss(2))),
    ("vprk-gauss2-symmetric",     iodeproblem(; timespan, timestep), SymmetricProjection(VPRKGauss(2))),
    ("vspark-glrk2-midpoint",     idaeproblem(; timespan, timestep), TableauVSPARKGLRKpMidpoint(2)),
    ("vspark-glrk2-symmetric",    idaeproblem(; timespan, timestep), TableauVSPARKGLRKpSymmetric(2)),
)

for (name, problem, method) in runs
    println("Integrating lotka_volterra_2d with $(name) ...")
    solution = integrate(problem, method)
    save(joinpath(@__DIR__, "lotka_volterra_2d-$(name).pdf"), plot_traces(solution, problem))
end
