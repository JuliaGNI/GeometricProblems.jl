using GeometricIntegrators
using GeometricProblems.KuboOscillator
using Test

@testset "$(rpad("Kubo Oscillator",80))" begin

    sde   = kubo_oscillator_sde_1()
    psde  = kubo_oscillator_psde_1()
    spsde = kubo_oscillator_spsde_1()

    sde_equs = get_function_tuple(sde)
    @test_nowarn sde_equs[:v](sde.t₀, sde.q₀, zero(sde.q₀))
    @test_nowarn sde_equs[:B](sde.t₀, sde.q₀, zeros(eltype(sde.q₀), sde.d, sde.m))
    @test_nowarn sde_equs[:B](sde.t₀, sde.q₀, zeros(eltype(sde.q₀), sde.d, sde.m), 1)

    psde_equs = get_function_tuple(psde)
    @test_nowarn psde_equs[:v](psde.t₀, psde.q₀, psde.p₀, zero(psde.q₀))
    @test_nowarn psde_equs[:f](psde.t₀, psde.q₀, psde.p₀, zero(psde.p₀))
    @test_nowarn psde_equs[:B](psde.t₀, psde.q₀, psde.p₀, zero(psde.q₀))
    @test_nowarn psde_equs[:G](psde.t₀, psde.q₀, psde.p₀, zero(psde.p₀))
    @test_nowarn psde_equs[:B](psde.t₀, psde.q₀, psde.p₀, zeros(eltype(psde.q₀), psde.d, psde.m))
    @test_nowarn psde_equs[:G](psde.t₀, psde.q₀, psde.p₀, zeros(eltype(psde.p₀), psde.d, psde.m))

    spsde_equs = get_function_tuple(spsde)
    @test_nowarn spsde_equs[:v](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.q₀))
    @test_nowarn spsde_equs[:f1](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:f2](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:B](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.q₀))
    @test_nowarn spsde_equs[:G1](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:G2](spsde.t₀, spsde.q₀, spsde.p₀, zero(spsde.p₀))
    @test_nowarn spsde_equs[:B](spsde.t₀, spsde.q₀, spsde.p₀, zeros(eltype(spsde.q₀), spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G1](spsde.t₀, spsde.q₀, spsde.p₀, zeros(eltype(spsde.p₀), spsde.d, spsde.m))
    @test_nowarn spsde_equs[:G2](spsde.t₀, spsde.q₀, spsde.p₀, zeros(eltype(spsde.p₀), spsde.d, spsde.m))

end
