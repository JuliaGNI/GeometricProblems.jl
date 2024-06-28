using GeometricIntegrators: ImplicitMidpoint, integrate
import GeometricProblems.CoupledHarmonicOscillator as cho
import GeometricProblems.TodaLattice as tl
using Test

function test_multiple_initial_conditions(cho::Module)
    q₀_vec = [cho.q₀ .+ α for α in 0. : .4 : .8]
    p₀_vec = [cho.p₀ .+ α for α in 0. : .4 : .8]

    # ensemble problem
    epr = cho.hodeensemble(q₀_vec, p₀_vec)
    # ensemble solution
    esol = integrate(epr, ImplicitMidpoint())

    sol = integrate(cho.hodeproblem(), ImplicitMidpoint())
    @test esol.s[1].q.d.parent ≈ sol.q.d.parent
    @test esol.s[2].q.d.parent ≉ sol.q.d.parent
end

function test_multiple_parameters(cho::Module)
    params_vec = Vector{NamedTuple}()

    for i in 0:10
        param_vals = ()
        for key in keys(cho.default_parameters)
            param_vals = (param_vals..., cho.default_parameters[key] .+ 1. * i)
        end
        params_vec = push!(params_vec, NamedTuple{keys(cho.default_parameters)}(param_vals))
    end

    # ensemble problem
    epr = cho.hodeensemble(; parameters = params_vec)
    # ensemble solution
    esol = integrate(epr, ImplicitMidpoint())

    sol = integrate(cho.hodeproblem(), ImplicitMidpoint())
    @test esol.s[1].q.d.parent ≈ sol.q.d.parent
    @test esol.s[2].q.d.parent ≉ sol.q.d.parent
end

test_multiple_initial_conditions(cho)
test_multiple_initial_conditions(tl)

test_multiple_parameters(cho)
test_multiple_parameters(tl)