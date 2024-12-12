using GeometricIntegrators: ImplicitMidpoint, integrate
import GeometricProblems.CoupledHarmonicOscillator as cho
using Test


q₀_vec = [cho.q₀ .+ α for α in 0. : .4 : .4]
p₀_vec = [cho.p₀ .+ α for α in 0. : .4 : .4]

# ensemble problem
epr = cho.hodeensemble(q₀_vec, p₀_vec)

# ensemble solution
esol = integrate(epr, ImplicitMidpoint())

@test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent


function _params(i)
    NamedTuple{keys(cho.default_parameters)}(values(cho.default_parameters) .+ i)
end

# ensemble problem
epr = cho.hodeensemble(; parameters = [_params(i) for i in 0:1])

# ensemble solution
esol = integrate(epr, ImplicitMidpoint())

@test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent
