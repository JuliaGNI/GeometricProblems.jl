@doc raw"""
# Mathematical Pendulum

"""
module Pendulum

    using GeometricEquations

    # export odeproblem, podeproblem, iodeproblem, idaeproblem

    # function pendulum_ode_v(v, t, q, params)
    #     v[1] = q[2]
    #     v[2] = sin(q[1])
    #     nothing
    # end

    # function odeproblem(q₀=[acos(0.4), 0.0])
    #     ODE(pendulum_ode_v, q₀)
    # end


    # function pendulum_pode_v(v, t, q, p, params)
    #     v[1] = p[1]
    #     nothing
    # end

    # function pendulum_pode_f(f, t, q, p, params)
    #     f[1] = sin(q[1])
    #     nothing
    # end

    # function podeproblem(q₀=[acos(0.4)], p₀=[0.0])
    #     PODE(pendulum_pode_v, pendulum_pode_f, q₀, p₀)
    # end


    # function pendulum_iode_α(p, t, q, v, params)
    #     p[1] = q[2]
    #     p[2] = 0
    #     nothing
    # end

    # function pendulum_iode_f(f, t, q, v, params)
    #     f[1] = sin(q[1])
    #     f[2] = v[1] - q[2]
    #     nothing
    # end

    # function pendulum_iode_g(g, t, q, λ, params)
    #     g[1] = 0
    #     g[2] = λ[1]
    #     nothing
    # end

    # function pendulum_iode_v(v, t, q, p, params)
    #     v[1] = q[2]
    #     v[2] = sin(q[1])
    #     nothing
    # end

    # function iodeproblem(q₀=[acos(0.4), 0.0], p₀=[0.0, 0.0])
    #     IODE(pendulum_iode_α, pendulum_iode_f,
    #          pendulum_iode_g, q₀, p₀; v̄=pendulum_iode_v)
    # end


    # # function pendulum_pdae_v(t, q, p, v)
    # #     v[1] = 0
    # #     nothing
    # # end
    # #
    # # function pendulum_pdae_f(t, q, p, f)
    # #     f[1] = sin(q[1])
    # #     nothing
    # # end
    # #
    # # function pendulum_pdae_u(t, q, p, λ, u)
    # #     u[1] = λ[1]
    # #     nothing
    # # end
    # #
    # # function pendulum_pdae_g(t, q, p, λ, g)
    # #     g[1] = 0
    # #     nothing
    # # end
    # #
    # # # TODO
    # # function pendulum_pdae_ϕ(t, q, p, ϕ)
    # #     ϕ[1] = p[1] - q[2]
    # #     nothing
    # # end
    # #
    # # # TODO
    # # function pendulum_pdae(q₀=[acos(0.4)], p₀=[0.0], λ₀=[0.0, 0.0])
    # #     PDAE(pendulum_pdae_v, pendulum_pdae_f, pendulum_pdae_u, pendulum_pdae_g, pendulum_pdae_ϕ, q₀, p₀, λ₀)
    # # end


    # function pendulum_idae_u(u, t, q, p, λ, params)
    #     u[1] = λ[1]
    #     u[2] = λ[2]
    #     nothing
    # end

    # function pendulum_idae_g(g, t, q, p, λ, params)
    #     g[1] = 0
    #     g[2] = λ[1]
    #     nothing
    # end

    # function pendulum_idae_ϕ(ϕ, t, q, p, params)
    #     ϕ[1] = p[1] - q[2]
    #     ϕ[2] = p[2]
    #     nothing
    # end

    # function idaeproblem(q₀=[acos(0.4), 0.0], p₀=[0.0, 0.0], λ₀=[0.0, 0.0])
    #     IDAE(pendulum_iode_f, pendulum_iode_α,
    #          pendulum_idae_u, pendulum_idae_g,
    #          pendulum_idae_ϕ, q₀, p₀, λ₀;
    #          v̄=pendulum_iode_v)
    # end

    #maybe inaccording to the above
    using EulerLagrange
    using LinearAlgebra
    using Parameters


    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const default_parameters= (
        m = 1.0,
        ω = 0.5,
    )
    
    const p₀ = [0.5]
    const q₀ = [0.5]
    const x₀ = vcat(q₀, p₀)

    const timestep = 0.01
    const timespan = (0.0, 10.0)

    function hamiltonian(t,q,p,params)
        @unpack m,ω = params
        (p[1]^2)/(2*m) - ω^2 * cos(q[1])
    end

    function lagrangian(t,q,v,params)
        @unpack m,ω = params
        m*(v[1]^2)/2 + ω^2 * cos(q[1]) 
    end

    # function θ̇(t, q, p, params)
    #     @unpack m = params
    #     - U * sin(q[1])
    # end

    function θ̇(t, q, p, params) # \partial H / \partial p  = \dot q =: v
        @unpack m,ω = params
        p[1]/m
    end


    function θ̇(v, t, q, p, params)
        v[1] = θ̇(t, q, p, params)
        nothing
    end


    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(1)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, x, v = lagrangian_variables(1)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, timespan, timestep, q₀, p₀; v̄ = θ̇, parameters = parameters)
    end


end
