# Massless charged particle in 2d: a survey of integrators
#
# The massless limit makes the Lagrangian degenerate by construction — the momentum is the magnetic
# vector potential, not the velocity — so this is the problem the SPARK and variational families
# were built for, and the one where the choice of projection is visible in the energy error.
# Writes one PDF per method next to this script.
#
# The pairings below are the ones exercised in `test/massless_charged_particle_tests.jl`; note that
# the SPARK method takes `idaeproblem_spark` rather than `idaeproblem`.

using CairoMakie
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.MasslessChargedParticle

# The module default is t ∈ [0, 1000] at Δt = 0.2, which is five thousand steps.
const timespan = (0.0, 1000.0)
const timestep = 0.2

const runs = (
    ("gauss8",                 odeproblem(;       timespan, timestep), Gauss(8)),
    ("vspark-glrk2-midpoint",  idaeproblem(;      timespan, timestep), TableauVSPARKGLRKpMidpoint(2)),
    ("spark-glrk2",            idaeproblem_spark(; timespan, timestep), SPARKGLRK(2)),
    ("vprk-gauss2-midpoint",   lodeproblem(;      timespan, timestep), MidpointProjection(VPRKGauss(2))),
    ("slrk-lobatto-iiie2",     ldaeproblem(;      timespan, timestep), SLRKLobattoIIIE(2)),
)

for (name, problem, method) in runs
    println("Integrating massless_charged_particle_2d with $(name) ...")
    solution = integrate(problem, method)
    save("massless_charged_particle_2d-$(name).pdf", plot_traces(solution, problem))
end
