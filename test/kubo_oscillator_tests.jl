using GeometricEquations
using GeometricIntegrators
using GeometricProblems.KuboOscillator
using Test

@testset "$(rpad("Kubo Oscillator",80))" begin

    sde   = kubo_oscillator_sde_1()
    psde  = kubo_oscillator_psde_1()
    spsde = kubo_oscillator_spsde_1()

    sde_equs = functions(sde)
    @test_nowarn sde_equs[:v](tbegin(sde), sde.ics.q, zero(sde.ics.q))
    @test_nowarn sde_equs[:B](tbegin(sde), sde.ics.q, zeros(eltype(sde.ics.q), sde.d, sde.m))
    @test_nowarn sde_equs[:B](tbegin(sde), sde.ics.q, zeros(eltype(sde.ics.q), sde.d, sde.m), 1)

    psde_equs = functions(psde)
    @test_nowarn psde_equs[:v](tbegin(psde), psde.ics.q, psde.ics.p, zero(psde.ics.q))
    @test_nowarn psde_equs[:f](tbegin(psde), psde.ics.q, psde.ics.p, zero(psde.ics.p))
    @test_nowarn psde_equs[:B](tbegin(psde), psde.ics.q, psde.ics.p, zero(psde.ics.q))
    @test_nowarn psde_equs[:G](tbegin(psde), psde.ics.q, psde.ics.p, zero(psde.ics.p))
    @test_nowarn psde_equs[:B](tbegin(psde), psde.ics.q, psde.ics.p, zeros(eltype(psde.ics.q), psde.d, psde.m))
    @test_nowarn psde_equs[:G](tbegin(psde), psde.ics.q, psde.ics.p, zeros(eltype(psde.ics.p), psde.d, psde.m))

    spsde_equs = functions(spsde)
    @test_nowarn spsde_equs[:v ](tbegin(spsde), spsde.ics.q, spsde.ics.p, zero(spsde.ics.q))
    @test_nowarn spsde_equs[:f1](tbegin(spsde), spsde.ics.q, spsde.ics.p, zero(spsde.ics.p))
    @test_nowarn spsde_equs[:f2](tbegin(spsde), spsde.ics.q, spsde.ics.p, zero(spsde.ics.p))
    @test_nowarn spsde_equs[:B ](tbegin(spsde), spsde.ics.q, spsde.ics.p, zero(spsde.ics.q))
    @test_nowarn spsde_equs[:G1](tbegin(spsde), spsde.ics.q, spsde.ics.p, zero(spsde.ics.p))
    @test_nowarn spsde_equs[:G2](tbegin(spsde), spsde.ics.q, spsde.ics.p, zero(spsde.ics.p))
    @test_nowarn spsde_equs[:B ](tbegin(spsde), spsde.ics.q, spsde.ics.p, zeros(eltype(spsde.ics.q), spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G1](tbegin(spsde), spsde.ics.q, spsde.ics.p, zeros(eltype(spsde.ics.p), spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G2](tbegin(spsde), spsde.ics.q, spsde.ics.p, zeros(eltype(spsde.ics.p), spsde.d, spsde.m))

end
