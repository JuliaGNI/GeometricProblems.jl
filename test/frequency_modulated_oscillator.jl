using Test
using GeometricIntegrators: Gauss, integrate, relative_maximum_error

import GeometricProblems.FrequencyModulatedOscillator as fmo

@testset "$(rpad("Frequency-Modulated Oscillator",80))" begin
    @test_nowarn fmo.hodeproblem()
    @test_nowarn fmo.lodeproblem()

    # The Hamiltonian and Lagrangian formulations describe the same non-autonomous system.
    hsol = integrate(fmo.hodeproblem(; timespan = (0.0, 10.0)), Gauss(8))
    lsol = integrate(fmo.lodeproblem(; timespan = (0.0, 10.0)), Gauss(8))
    @test relative_maximum_error(hsol.q, lsol.q) < 1E-12
    @test relative_maximum_error(hsol.p, lsol.p) < 1E-12
end

@testset "$(rpad("Frequency-Modulated Oscillator (adiabatic invariant)",80))" begin
    # Over one modulation period 2π/ε the action J = H/ω deviates from its initial value by O(ε):
    # reducing ε tenfold reduces the maximal deviation tenfold.
    ΔJ = map((1E-1, 1E-2)) do ε
        params = merge(fmo.default_parameters(), (ε = ε,))
        sol = integrate(
            fmo.hodeproblem(; timespan = (0.0, 2π / ε), parameters = params), Gauss(4))
        J = [fmo.adiabatic_invariant(sol.t[n], sol.q[n], sol.p[n], params)
             for n in eachindex(sol.t)]
        maximum(abs.(J .- J[begin])) / J[begin]
    end
    @test 8 < ΔJ[1] / ΔJ[2] < 12
    @test ΔJ[2] < 1E-2
end
