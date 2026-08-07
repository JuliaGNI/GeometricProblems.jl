# Lotka-Volterra 2d with a singular Lagrangian: a survey of integrators
#
# Same system and same methods as `lotka_volterra_2d.jl`, but built from a Lagrangian whose
# Hessian is singular, so the implicit formulations are genuinely degenerate. This is the case
# variational and SPARK integrators exist for, and the one where an unprojected method drifts.
#
# `LotkaVolterra2dSingular` has no plot extension of its own — only `LotkaVolterra2d` does — but
# those plot functions are generic in the solution and the problem, so they are reused here.

using CairoMakie
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra2dSingular

import GeometricProblems.LotkaVolterra2d: plot_traces

const timespan = (0.0, 10.0)
const timestep = 0.01

const runs = (
    ("gauss2",                 odeproblem(;  timespan, timestep), Gauss(2)),
    ("vprk-gauss2-midpoint",   iodeproblem(; timespan, timestep), MidpointProjection(VPRKGauss(2))),
    ("vprk-gauss2-symmetric",  iodeproblem(; timespan, timestep), SymmetricProjection(VPRKGauss(2))),
    ("vspark-glrk2-midpoint",  idaeproblem(; timespan, timestep), TableauVSPARKGLRKpMidpoint(2)),
    ("vspark-glrk2-symmetric", idaeproblem(; timespan, timestep), TableauVSPARKGLRKpSymmetric(2)),
)

for (name, problem, method) in runs
    println("Integrating lotka_volterra_2d_singular with $(name) ...")
    solution = integrate(problem, method)
    save(joinpath(@__DIR__, "lotka_volterra_2d_singular-$(name).pdf"), plot_traces(solution, problem))
end
