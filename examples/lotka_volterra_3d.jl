# Lotka-Volterra 3d: a survey of integrators
#
# The 3d system is a Poisson system with a Casimir as well as a Hamiltonian, and it is provided
# only in its ODE form — there is no variational or DAE formulation here, so the survey is over
# Runge-Kutta methods of increasing order. Writes one PDF per method next to this script, showing
# the time traces together with the relative energy error.

using CairoMakie
using GeometricIntegrators
using GeometricProblems.LotkaVolterra3d

const timespan = (0.0, 10.0)
const timestep = 0.01

const runs = (
    ("implicit-midpoint", ImplicitMidpoint()),
    ("gauss1",            Gauss(1)),
    ("gauss2",            Gauss(2)),
    ("gauss8",            Gauss(8)),
)

const problem = odeproblem(; timespan, timestep)

for (name, method) in runs
    println("Integrating lotka_volterra_3d with $(name) ...")
    solution = integrate(problem, method)
    save(joinpath(@__DIR__, "lotka_volterra_3d-$(name).pdf"), plot_traces(solution, problem))
end
