using Test
using GeometricIntegrators
using GeometricProblems.Pendulum
using GeometricSolutions


@testset "$(rpad("Pendulum",80))" begin

    @test_nowarn odeproblem()
    @test_nowarn podeproblem()
    @test_nowarn hodeproblem()
    @test_nowarn iodeproblem()
    @test_nowarn idaeproblem()

    @test_nowarn hodeensemble()


    ode  = odeproblem()
    pode = podeproblem()
    hode = hodeproblem()
    iode = iodeproblem()

    ode_sol  = integrate(ode,  Gauss(2))
    pode_sol = integrate(pode, Gauss(2))
    hode_sol = integrate(hode, Gauss(2))
    iode_sol = integrate(iode, MidpointProjection(VPRKGauss(2)))

    # ODE state is [θ, θ̇]; PODE has q=[θ], p=[θ̇]: compare final values
    @test ode_sol.q[end][1] ≈ pode_sol.q[end][1]
    @test ode_sol.q[end][2] ≈ pode_sol.p[end][1]

    # HODE shares the canonical vector fields with PODE: solutions must agree
    @test hode_sol.q[end][1] ≈ pode_sol.q[end][1]
    @test hode_sol.p[end][1] ≈ pode_sol.p[end][1]

    # IODE degenerate state is the same [θ, θ̇] as ODE
    @test relative_maximum_error(ode_sol.q, iode_sol.q) < 1E-4

    iode_sol2 = integrate(iode, SymmetricProjection(VPRKGauss(2)))
    @test relative_maximum_error(iode_sol.q, iode_sol2.q) < 1E-4

end


@testset "$(rpad("Pendulum HODE Ensemble",80))" begin

    # maximum relative energy error of a solution (the symplectic Gauss method
    # conserves the Hamiltonian up to a small bounded error)
    function max_energy_error(sol, params)
        h₀ = hamiltonian(sol[0].t, sol[0].q, sol[0].p, params)
        maximum(abs((hamiltonian(sol[n].t, sol[n].q, sol[n].p, params) - h₀) / h₀) for n in eachtimestep(sol))
    end

    q̄ = [1.0]
    p̄ = [0.5]

    # (i) different initial conditions, identical parameters
    ens = hodeensemble()
    @test length(ens) == 100
    @test all(p == default_parameters() for p in ens.parameters)            # parameters identical
    @test length(unique([(ic.q[1], ic.p[1]) for ic in ens.ics])) == 100     # initial conditions vary

    sols = integrate(ens, Gauss(2))
    @test length(sols) == length(ens)
    @test all(max_energy_error(sols[i], ens.parameters[i]) < 1E-4 for i in eachindex(ens.ics))


    # (ii) identical initial condition, different parameters (varying length l)
    lengths = [0.5, 1.0, 2.0, 4.0]
    params  = [(l = L, m = 1.0, g = 1.0) for L in lengths]

    ens = hodeensemble(q̄, p̄, params)
    @test length(ens) == length(lengths)
    @test all(ic.q == ens.ics[begin].q && ic.p == ens.ics[begin].p for ic in ens.ics)  # initial conditions identical
    @test [p.l for p in ens.parameters] == lengths                                     # lengths vary

    sols = integrate(ens, Gauss(2))
    @test length(sols) == length(ens)
    @test all(max_energy_error(sols[i], ens.parameters[i]) < 1E-4 for i in eachindex(ens.ics))
    # different lengths must produce different trajectories
    @test length(unique([sols[i].q[end][1] for i in eachindex(ens.ics)])) == length(lengths)


    # (iii) different initial conditions AND different parameters simultaneously
    #       a 2×2 grid of initial conditions paired with the four parameter sets
    ens = hodeensemble([0.5], [1.5], [-0.5], [0.5], [2], [2]; parameters=params)
    @test length(ens) == length(lengths)
    @test length(unique([(ic.q[1], ic.p[1]) for ic in ens.ics])) == length(lengths)    # initial conditions vary
    @test [p.l for p in ens.parameters] == lengths                                     # lengths vary

    sols = integrate(ens, Gauss(2))
    @test length(sols) == length(ens)
    @test all(max_energy_error(sols[i], ens.parameters[i]) < 1E-4 for i in eachindex(ens.ics))

end
