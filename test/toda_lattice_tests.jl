using GeometricIntegrators: ImplicitMidpoint, integrate
import GeometricProblems.TodaLattice as tl
using Test

N = 20
q₀ = tl.compute_initial_q(tl.μ, N)
p₀ = zero(q₀)

q₀_vec = [q₀ .+ α for α in 0. : .4 : .4]
p₀_vec = [p₀ .+ α for α in 0. : .4 : .4]

# ensemble problem
epr = tl.hodeensemble(N, q₀_vec, p₀_vec)

# ensemble solution
esol = integrate(epr, ImplicitMidpoint())

@test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent


function _params(i)
    NamedTuple{keys(tl.default_parameters)}(values(tl.default_parameters) .+ i)
end

# ensemble problem
epr = tl.hodeensemble(N; parameters = [_params(i) for i in 0:1])

# ensemble solution
esol = integrate(epr, ImplicitMidpoint())

@test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent
