using Test
using GeometricIntegrators
using GeometricSolutions
using GeometricProblems

lode = GeometricProblems.PerturbedPendulum.lodeproblem(tspan = (0,50))
hode = hodeproblem()

hode_sol = integrate(hode, Gauss(2))
lode_sol = integrate(lode, Gauss(2))

@test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-12

using Plots
plot(lode_sol.q[:,1])

pen_lode = GeometricProblems.Pendulum.lodeproblem([0.5],[0.0],tspan = (0,50))
pen_sol = integrate(pen_lode, Gauss(2))
plot!(pen_sol.q[:,1])