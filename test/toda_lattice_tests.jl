using GeometricIntegrators: Gauss, integrate, relative_maximum_error
using GeometricProblems.TodaLattice
using Test


# Parameters and initial conditions

N = 20
μ = TodaLattice.μ
q₀ = TodaLattice.compute_initial_q(μ, N)
p₀ = zero(q₀)
params = TodaLattice.default_parameters


# Ensemble initial conditions and parameters

q₀_vec = [q₀ .+ α for α in -0.2 : 0.2 : +0.2]
p₀_vec = [p₀ .+ α for α in -0.2 : 0.2 : +0.2]

param_vec = [NamedTuple{keys(params)}(values(params) .+ β) for β in -0.1 : 0.1 : +0.1]


# HODE and LODE problems

hode_prb = hodeproblem(q₀, p₀)
lode_prb = lodeproblem(q₀, p₀)

href_sol = integrate(hode_prb, Gauss(1))
lref_sol = integrate(lode_prb, Gauss(1))

@test relative_maximum_error(href_sol.q, lref_sol.q) < 2E-14


# Ensemble problems with different initial conditions
hode_ens = hodeensemble(N, q₀_vec, p₀_vec)
lode_ens = lodeensemble(N, q₀_vec, p₀_vec)

hode_sol = integrate(hode_ens, Gauss(1))
lode_sol = integrate(lode_ens, Gauss(1))

@test relative_maximum_error(hode_sol[2].q, href_sol.q) < 2E-14
@test relative_maximum_error(lode_sol[2].q, lref_sol.q) < 2E-14

for (hsol,lsol) in zip(hode_sol,lode_sol)
    @test relative_maximum_error(hsol.q, lsol.q) < 2E-14
end


# Ensemble problems with different parameters
hode_ens = hodeensemble(N; parameters = param_vec)
lode_ens = lodeensemble(N; parameters = param_vec)

hode_sol = integrate(hode_ens, Gauss(1))
lode_sol = integrate(lode_ens, Gauss(1))

@test relative_maximum_error(hode_sol[2].q, href_sol.q) < 2E-14
@test relative_maximum_error(lode_sol[2].q, lref_sol.q) < 2E-14

for (hsol,lsol) in zip(hode_sol,lode_sol)
    @test relative_maximum_error(hsol.q, lsol.q) < 2E-14
end


# Ensemble problems with different initial conditions and parameters
hode_ens = hodeensemble(q₀_vec, p₀_vec; parameters = param_vec)
lode_ens = lodeensemble(q₀_vec, p₀_vec; parameters = param_vec)

hode_sol = integrate(hode_ens, Gauss(1))
lode_sol = integrate(lode_ens, Gauss(1))

@test relative_maximum_error(hode_sol[2].q, href_sol.q) < 2E-14
@test relative_maximum_error(lode_sol[2].q, lref_sol.q) < 2E-14

for (hsol,lsol) in zip(hode_sol,lode_sol)
    @test relative_maximum_error(hsol.q, lsol.q) < 2E-14
end
