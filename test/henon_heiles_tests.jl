using Test
using GeometricIntegrators
using GeometricProblems.HenonHeilesPotential
using GeometricSolutions

lode = GeometricProblems.HenonHeilesPotential.lodeproblem(tspan =(0,300))
hode = GeometricProblems.HenonHeilesPotential.hodeproblem(tspan =(0,300))

hode_sol = integrate(hode, ImplicitMidpoint())
lode_sol = integrate(lode, ImplicitMidpoint())

@test relative_maximum_error(hode_sol[1].q, lode_sol[1].q) < 1E-3
@test relative_maximum_error(hode_sol[1].p, lode_sol[1].p) < 1E-3
 

using Plots
plot(lode_sol[1].q[:,1], label = "lode_q₁", size = (800,200))
plot!(hode_sol[1].q[:,1], label = "hode_q₁")

plot(lode_sol[1].q[:,2], label = "lode_q₂", size = (800,200))
plot!(hode_sol[1].q[:,2], label = "hode_q₂")