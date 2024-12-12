using GeometricIntegrators: ImplicitMidpoint, integrate
using GeometricProblems.TodaLattice
using Test

# Initial conditions

N = 20
μ = TodaLattice.μ
q₀ = TodaLattice.compute_initial_q(μ, N)
p₀ = zero(q₀)
params = TodaLattice.default_parameters

# Ensemble initial conditions and parameters

q₀_vec = [q₀ .+ α for α in 0. : .4 : .4]
p₀_vec = [p₀ .+ α for α in 0. : .4 : .4]

param_vec = [NamedTuple{keys(params)}(values(params) .+ i) for i in 0:1]


# Ensemble problem with different initial conditions
epr = hodeensemble(N, q₀_vec, p₀_vec)

# ensemble solution
esol = integrate(epr, ImplicitMidpoint())

@test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent


# Ensemble problem with different parameters
epr = hodeensemble(N; parameters = param_vec)

# ensemble solution
esol = integrate(epr, ImplicitMidpoint())

@test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent


# Ensemble problem with different initial conditions and parameters
epr = hodeensemble(q₀_vec, p₀_vec; parameters = param_vec)

# ensemble solution
esol = integrate(epr, ImplicitMidpoint())

@test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent
