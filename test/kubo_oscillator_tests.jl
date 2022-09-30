using GeometricIntegrators
using GeometricProblems.KuboOscillator
using Test

@testset "$(rpad("Kubo Oscillator",80))" begin

    sde   = kubo_oscillator_sde_1()
    psde  = kubo_oscillator_psde_1()
    spsde = kubo_oscillator_spsde_1()

    sde_equs = functions(sde)
    @test_nowarn sde_equs[:v](zero(sde.ics.q), tbegin(sde), sde.ics.q)
    @test_nowarn sde_equs[:B](zeros(eltype(sde.ics.q), tbegin(sde), sde.ics.q, sde.d, sde.m))
    @test_nowarn sde_equs[:B](zeros(eltype(sde.ics.q), tbegin(sde), sde.ics.q, sde.d, sde.m), 1)

    psde_equs = functions(psde)
    @test_nowarn psde_equs[:v](zero(psde.ics.q), tbegin(psde), psde.ics.q, psde.ics.p)
    @test_nowarn psde_equs[:f](zero(psde.ics.p), tbegin(psde), psde.ics.q, psde.ics.p)
    @test_nowarn psde_equs[:B](zero(psde.ics.q), tbegin(psde), psde.ics.q, psde.ics.p)
    @test_nowarn psde_equs[:G](zero(psde.ics.p), tbegin(psde), psde.ics.q, psde.ics.p)
    @test_nowarn psde_equs[:B](zeros(eltype(psde.ics.q), tbegin(psde), psde.ics.q, psde.ics.p, psde.d, psde.m))
    @test_nowarn psde_equs[:G](zeros(eltype(psde.ics.p), tbegin(psde), psde.ics.q, psde.ics.p, psde.d, psde.m))

    spsde_equs = functions(spsde)
    @test_nowarn spsde_equs[:v ](zero(spsde.ics.q), tbegin(spsde), spsde.ics.q, spsde.ics.p)
    @test_nowarn spsde_equs[:f1](zero(spsde.ics.p), tbegin(spsde), spsde.ics.q, spsde.ics.p)
    @test_nowarn spsde_equs[:f2](zero(spsde.ics.p), tbegin(spsde), spsde.ics.q, spsde.ics.p)
    @test_nowarn spsde_equs[:B ](zero(spsde.ics.q), tbegin(spsde), spsde.ics.q, spsde.ics.p)
    @test_nowarn spsde_equs[:G1](zero(spsde.ics.p), tbegin(spsde), spsde.ics.q, spsde.ics.p)
    @test_nowarn spsde_equs[:G2](zero(spsde.ics.p), tbegin(spsde), spsde.ics.q, spsde.ics.p)
    @test_nowarn spsde_equs[:B ](zeros(eltype(spsde.ics.q), tbegin(spsde), spsde.ics.q, spsde.ics.p, spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G1](zeros(eltype(spsde.ics.p), tbegin(spsde), spsde.ics.q, spsde.ics.p, spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G2](zeros(eltype(spsde.ics.p), tbegin(spsde), spsde.ics.q, spsde.ics.p, spsde.d, spsde.m))

end
