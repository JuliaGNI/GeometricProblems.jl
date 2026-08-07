# Lotka-Volterra 4d: a survey of integrators
#
# The 4d system carries a degenerate Lagrangian, so it comes in a full set of formulations — ODE,
# implicit ODE, Lagrangian ODE, and the corresponding DAEs — and the point of this example is that
# the same trajectory can be reached through any of them. Writes one PDF per run next to this
# script, showing the time traces together with the relative energy error.
#
# The pairings below are the ones exercised in `test/lotka_volterra_4d_tests.jl`.

using CairoMakie
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra4d

const timespan = (0.0, 10.0)
const timestep = 0.01

const runs = (
    ("ode-gauss2",                 odeproblem(;  timespan, timestep), Gauss(2)),
    ("iode-gauss2",                iodeproblem(; timespan, timestep), Gauss(2)),
    ("lode-gauss2",                lodeproblem(; timespan, timestep), Gauss(2)),
    ("iode-vprk-gauss2-midpoint",  iodeproblem(; timespan, timestep), MidpointProjection(VPRKGauss(2))),
    ("iode-vprk-gauss2-symmetric", iodeproblem(; timespan, timestep), SymmetricProjection(VPRKGauss(2))),
    ("idae-vspark-glrk2-midpoint", idaeproblem(; timespan, timestep), TableauVSPARKGLRKpMidpoint(2)),
    ("ldae-vspark-glrk2-midpoint", ldaeproblem(; timespan, timestep), TableauVSPARKGLRKpMidpoint(2)),
)

for (name, problem, method) in runs
    println("Integrating lotka_volterra_4d with $(name) ...")
    solution = integrate(problem, method)
    save(joinpath(@__DIR__, "lotka_volterra_4d-$(name).pdf"), plot_traces(solution, problem))
end
