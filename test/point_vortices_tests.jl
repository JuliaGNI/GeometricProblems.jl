
module PointVorticesTests

    using Test
    using GeometricIntegrators
    using GeometricIntegrators.Utils
    using GeometricProblems.PointVortices


    const Δt = 0.01
    const nt = 1000

    const ref = [0.18722529318641928, 0.38967432450068706, 0.38125332930294187, 0.4258020604293123]


    @testset "$(rpad("Point Vortices",80))" begin
        ode  = point_vortices_ode()
        iode = point_vortices_iode()
        idae = point_vortices_idae()

        int = Integrator(ode, getTableauGLRK(1), Δt)
        sol = integrate(int, nt)
        @test rel_err(sol.q, ref) < 3E-2

        int = IntegratorVPRKpSymmetric(iode, getTableauVPGLRK(1), Δt)
        sol = integrate(int, nt)
        @test rel_err(sol.q, ref) < 3E-2

        int = Integrator(idae, getTableauGLRKpSymmetric(1), Δt)
        sol = integrate(int, nt)
        @test rel_err(sol.q, ref) < 3E-2
    end

end


using .PointVorticesTests
