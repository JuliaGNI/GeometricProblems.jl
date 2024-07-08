using GeometricIntegrators: ImplicitMidpoint, integrate
import GeometricProblems.CoupledHarmonicOscillator as cho
import GeometricProblems.TodaLattice as tl
using Test

function test_multiple_initial_conditions(cho::Module)
    q₀_vec = [cho.q₀ .+ α for α in 0. : .4 : .4]
    p₀_vec = [cho.p₀ .+ α for α in 0. : .4 : .4]

    # ensemble problem
    epr = cho.hodeensemble(q₀_vec, p₀_vec)
    # ensemble solution
    esol = integrate(epr, ImplicitMidpoint())

    @test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent
end

function test_multiple_parameters(cho::Module)
    params_vec = Vector{NamedTuple}()

    for i in 0:1
        param_vals = ()

        # take the first parameter (compute the ensemble by changing this parameter)
        key = keys(cho.default_parameters)[1]
    
        params_vec = push!(params_vec, NamedTuple{keys(cho.default_parameters)}(param_vals))
    end

    # ensemble problem
    epr = cho.hodeensemble(; parameters = params_vec)
    # ensemble solution
    esol = integrate(epr, ImplicitMidpoint())

    @test esol.s[2].q.d.parent ≉ esol.s[1].q.d.parent
end

test_multiple_initial_conditions(cho)
test_multiple_initial_conditions(tl)

test_multiple_parameters(cho)
test_multiple_parameters(tl)