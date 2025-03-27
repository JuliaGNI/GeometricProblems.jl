using Test
using GeometricIntegrators
using GeometricProblems.Pendulum
using GeometricSolutions

lode = lodeproblem()
hode = hodeproblem()

hode_sol = integrate(hode, Gauss(2))
lode_sol = integrate(lode, Gauss(2))

@test relative_maximum_error(hode_sol.q, lode_sol.q) < 1E-12
@test relative_maximum_error(hode_sol.p, lode_sol.p) < 1E-11
 