# Point vortices: a survey of integrators
#
# Four interacting point vortices in the plane. The system is noncanonical and its Lagrangian is
# degenerate, so as with the massless charged particle the interesting comparison is between plain
# Runge-Kutta on the ODE form and the projected/SPARK methods on the implicit forms. Writes one
# PDF per method next to this script, showing the time traces and the relative energy error.
#
# The pairings below are the ones exercised in `test/point_vortices_tests.jl`.

using CairoMakie
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.PointVortices

const timespan = (0.0, 10.0)
const timestep = 0.01

const runs = (
    ("gauss2",                 odeproblem(;  timespan, timestep), Gauss(2)),
    ("gauss8",                 odeproblem(;  timespan, timestep), Gauss(8)),
    ("iode-gauss2",            iodeproblem(; timespan, timestep), Gauss(2)),
    ("vspark-glrk2-symmetric", idaeproblem(; timespan, timestep), TableauVSPARKGLRKpSymmetric(2)),
)

for (name, problem, method) in runs
    println("Integrating point_vortices with $(name) ...")
    solution = integrate(problem, method)
    save("point_vortices-$(name).pdf", plot_traces(solution, problem))
end
