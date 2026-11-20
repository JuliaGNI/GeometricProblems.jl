using Test
using GeometricIntegrators
using GeometricIntegrators.SPARK
using GeometricProblems.LotkaVolterra4d
using GeometricSolutions
using StaticArrays: MVector, SVector

import GeometricProblems.LotkaVolterra4d as lv4

# The stage residual of an IODE in a kernel or an allocation-free solve writes `ϑ` into an
# `MVector`, which stays on the stack only while every call that receives it is inlined.
function stage_ϑ(t, q, v, params)
    Θ = MVector{4, eltype(q)}(undef)
    @inline lv4.lotka_volterra_4d_ϑ(Θ, t, q, v, params)
    return SVector(Θ)
end

function stage_ϑ_allocated(t, q, v, params)
    stage_ϑ(t, q, v, params)
    return @allocated stage_ϑ(t, q, v, params)
end

@testset "$(rpad("Lotka-Volterra 4D",80))" begin
    ode = odeproblem()
    iode = iodeproblem()
    lode = lodeproblem()
    pode = podeproblem()
    idae = idaeproblem()
    ldae = ldaeproblem()
    pdae = pdaeproblem()

    ldae_secondary = ldaeproblem_secondary()

    ref = integrate(ode, Gauss(8))

    sol = integrate(ode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(iode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(lode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(pode, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(idae, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(ldae, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 8E-4

    sol = integrate(pdae, Gauss(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(iode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(lode, MidpointProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(lode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(idae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(idae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(ldae, TableauVSPARKGLRKpMidpoint(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(ldae, TableauVSPARKGLRKpSymmetric(2))
    @test relative_maximum_error(sol.q, ref.q) < 1E-8

    sol = integrate(ldae_secondary, TableauVSPARKLobattoIIIE(2))
    @test relative_maximum_error(sol.q, ref.q) < 4E-5

    # TODO: Reactivate!
    # sol = integrate(pdae, TableauHSPARKLobattoIIIE(2))
    # @test relative_maximum_error(sol.q, ref.q) < 4E-4

end

@testset "Lotka-Volterra 4D: ϑ into an MVector stays on the stack" begin
    for T in (Float32, Float64)
        t = zero(T)
        q = SVector{4, T}(2.0, 1.5, 0.75, 1.25)
        v = SVector{4, T}(0.1, -0.2, 0.3, -0.4)
        params = lv4.default_parameters(T)
        @test stage_ϑ(t, q, v, params) == lv4.ϑ(t, Vector(q))
        @test stage_ϑ_allocated(t, q, v, params) == 0
    end
end
