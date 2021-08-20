@doc raw"""
# Mathematical Pendulum

"""
module Pendulum

    using GeometricEquations

    export pendulum_ode, pendulum_pode, pendulum_iode, pendulum_idae

    function pendulum_ode_f(t, x, f)
        f[1] = x[2]
        f[2] = sin(x[1])
        nothing
    end

    function pendulum_ode(x₀=[acos(0.4), 0.0])
        ODE(pendulum_ode_f, x₀)
    end


    function pendulum_pode_v(t, q, p, v)
        v[1] = p[1]
        nothing
    end

    function pendulum_pode_f(t, q, p, f)
        f[1] = sin(q[1])
        nothing
    end

    function pendulum_pode(q₀=[acos(0.4)], p₀=[0.0])
        PODE(pendulum_pode_v, pendulum_pode_f, q₀, p₀)
    end


    function pendulum_iode_α(t, q, v, p)
        p[1] = q[2]
        p[2] = 0
        nothing
    end

    function pendulum_iode_f(t, q, v, f)
        f[1] = sin(q[1])
        f[2] = v[1] - q[2]
        nothing
    end

    function pendulum_iode_g(t, q, λ, g)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function pendulum_iode_v(t, q, p, v)
        v[1] = q[2]
        v[2] = sin(q[1])
        nothing
    end

    function pendulum_iode(q₀=[acos(0.4), 0.0], p₀=[0.0, 0.0])
        IODE(pendulum_iode_α, pendulum_iode_f,
             pendulum_iode_g, q₀, p₀; v̄=pendulum_iode_v)
    end


    # function pendulum_pdae_v(t, q, p, v)
    #     v[1] = 0
    #     nothing
    # end
    #
    # function pendulum_pdae_f(t, q, p, f)
    #     f[1] = sin(q[1])
    #     nothing
    # end
    #
    # function pendulum_pdae_u(t, q, p, λ, u)
    #     u[1] = λ[1]
    #     nothing
    # end
    #
    # function pendulum_pdae_g(t, q, p, λ, g)
    #     g[1] = 0
    #     nothing
    # end
    #
    # # TODO
    # function pendulum_pdae_ϕ(t, q, p, ϕ)
    #     ϕ[1] = p[1] - q[2]
    #     nothing
    # end
    #
    # # TODO
    # function pendulum_pdae(q₀=[acos(0.4)], p₀=[0.0], λ₀=[0.0, 0.0])
    #     PDAE(pendulum_pdae_v, pendulum_pdae_f, pendulum_pdae_u, pendulum_pdae_g, pendulum_pdae_ϕ, q₀, p₀, λ₀)
    # end


    function pendulum_idae_u(t, q, p, λ, u)
        u[1] = λ[1]
        u[2] = λ[2]
        nothing
    end

    function pendulum_idae_g(t, q, p, λ, g)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function pendulum_idae_ϕ(t, q, p, ϕ)
        ϕ[1] = p[1] - q[2]
        ϕ[2] = p[2]
        nothing
    end

    function pendulum_idae(q₀=[acos(0.4), 0.0], p₀=[0.0, 0.0], λ₀=[0.0, 0.0])
        IDAE(pendulum_iode_f, pendulum_iode_α,
             pendulum_idae_u, pendulum_idae_g,
             pendulum_idae_ϕ, q₀, p₀, λ₀;
             v̄=pendulum_iode_v)
    end

end
