using Test
using GeometricIntegrators

t_end =5.
t_step = 0.5
lode_problem = OuterSolarSystem.lodeproblem(tspan = (0,t_end), tstep = t_step,n=3)
ImMi_sol = integrate(lode_problem, ImplicitMidpoint())
ImMi_sol.q


hode_problem = OuterSolarSystem.hodeproblem(tspan = (0,t_end), tstep = t_step,n=3)
sol = integrate(hode_problem, RK4())

relative_maximum_error(sol.q,ImMi_sol.q)
