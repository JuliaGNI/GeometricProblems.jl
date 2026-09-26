using GeometricIntegrators: Gauss, integrate, relative_maximum_error
using Test

import GeometricProblems.DuffingOscillator as duffing
import GeometricProblems.LennardJonesOscillator as lj
import GeometricProblems.MorseOscillator as morse
import GeometricProblems.MathewsLakshmananOscillator as ml

# For each oscillator the Hamiltonian (HODE) and Lagrangian (LODE) formulations must describe the
# same system, so integrating both from the same initial condition with the same method has to
# yield the same trajectory (up to solver tolerance). This exercises the Legendre transform /
# consistency of the newly implemented models.

@testset "$(rpad("Duffing oscillator",80))" begin
    href = integrate(duffing.hodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    lref = integrate(duffing.lodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    @test relative_maximum_error(href.q, lref.q) < 1E-10
end

@testset "$(rpad("Lennard-Jones oscillator",80))" begin
    href = integrate(lj.hodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    lref = integrate(lj.lodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    @test relative_maximum_error(href.q, lref.q) < 1E-10
end

@testset "$(rpad("Morse oscillator",80))" begin
    href = integrate(morse.hodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    lref = integrate(morse.lodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    @test relative_maximum_error(href.q, lref.q) < 1E-10
end

@testset "$(rpad("Mathews-Lakshmanan oscillator",80))" begin
    href = integrate(ml.hodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    lref = integrate(ml.lodeproblem(; timespan = (0.0, 1.0)), Gauss(2))
    @test relative_maximum_error(href.q, lref.q) < 1E-10
end
