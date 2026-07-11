using GeometricEquations: HODEProblem, LODEProblem
using GeometricProblems.LinearWave
using Test

# Regression test for the linear wave equation. Previously `lodeproblem` called a lowercase
# `lodeproblem(::LagrangianSystem, …)` (itself) and threw a `MethodError`; it must now build a
# proper `LODEProblem` via EulerLagrange, mirroring `hodeproblem`.

@testset "$(rpad("Linear Wave",80))" begin
    hode = hodeproblem()
    lode = lodeproblem()

    @test hode isa HODEProblem
    @test lode isa LODEProblem
end
