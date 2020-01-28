module LotkaVolterra2dsingular

    using GeometricIntegrators.Equations
    using GeometricIntegrators.Solutions

    export ϑ₁, ϑ₂, hamiltonian
    export lotka_volterra_2d_ode, lotka_volterra_2d_iode, lotka_volterra_2d_idae
    export compute_energy_error, compute_momentum_error


    const A1=1.0
    const A2=3.0
    const B1=1.0
    const B2=5.0

    const X0=1.0
    const Y0=1.0


    function ϑ₁(t, q)
        q[2] + log(q[2]) / q[1]
    end

    function ϑ₂(t, q)
        zero(eltype(q))
    end


    function f₁(t, q, v)
        - v[1] * log(q[2]) / q[1]^2
    end

    function f₂(t, q, v)
        v[1] * (1 + 1 / (q[1] * q[2]))
    end


    const q₀=[X0, Y0]
    const p₀=[ϑ₁(0, q₀), ϑ₂(0, q₀)]


    function hamiltonian(t, q)
        A1*q[1] + A2*q[2] - B1*log(q[1]) - B2*log(q[2])
    end


    function lotka_volterra_2d_ode_v(t, q, v)
        v[1] = q[1] / (1 + q[1]*q[2]) * (A2 * q[2] - B2)
        v[2] = q[2] / (1 + q[1]*q[2]) * (B1 - A1 * q[1])
        nothing
    end

    function lotka_volterra_2d_ode(q₀=q₀)
        ODE(lotka_volterra_2d_ode_v, q₀; h=hamiltonian)
    end


    function lotka_volterra_2d_iode_ϑ(t, q, p)
        p[1] = ϑ₁(t,q)
        p[2] = ϑ₂(t,q)
        nothing
    end

    function lotka_volterra_2d_iode_ϑ(t, q, v, p)
        lotka_volterra_2d_iode_ϑ(t, q, p)
    end

    function lotka_volterra_2d_iode_f(t, q, v, f)
        f[1] = f₁(t,q,v) - A1 + B1 / q[1]
        f[2] = f₂(t,q,v) - A2 + B2 / q[2]
        nothing
    end

    function lotka_volterra_2d_iode_g(t, q, λ, g)
        g[1] = f₁(t,q,λ)
        g[2] = f₂(t,q,λ)
        nothing
    end

    function lotka_volterra_2d_iode_v(t, q, p, v)
        lotka_volterra_2d_ode_v(t, q, v)
    end

    function lotka_volterra_2d_iode(q₀=q₀, p₀=p₀)
        IODE(lotka_volterra_2d_iode_ϑ, lotka_volterra_2d_iode_f,
             lotka_volterra_2d_iode_g, q₀, p₀;
             h=hamiltonian, v=lotka_volterra_2d_iode_v)
    end

    function lotka_volterra_2d_idae_u(t, q, p, λ, u)
        u[1] = λ[1]
        u[2] = λ[2]
        nothing
    end

    function lotka_volterra_2d_idae_g(t, q, p, λ, g)
        g[1] = f₁(t,q,λ)
        g[2] = f₂(t,q,λ)
        nothing
    end

    function lotka_volterra_2d_idae_ϕ(t, q, p, ϕ)
        ϕ[1] = p[1] - ϑ₁(t,q)
        ϕ[2] = p[2] - ϑ₂(t,q)
        nothing
    end

    function lotka_volterra_2d_idae(q₀=q₀, p₀=p₀, λ₀=zero(q₀))
        IDAE(lotka_volterra_2d_iode_ϑ, lotka_volterra_2d_iode_f,
             lotka_volterra_2d_idae_u, lotka_volterra_2d_idae_g,
             lotka_volterra_2d_idae_ϕ, q₀, p₀, λ₀;
             v=lotka_volterra_2d_iode_v)
    end


    function compute_energy_error(t, q::DataSeries{T}) where {T}
        h = SDataSeries(T, q.nt)
        e = SDataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i])
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_momentum_error(t, q, p)
        p1_error = zeros(q.nt+1)
        p2_error = zeros(q.nt+1)

        for i in 1:(q.nt+1)
            p1_error[i] = p.d[1,i] - ϑ₁(t.t[i], q.d[:,i])
            p2_error[i] = p.d[2,i] - ϑ₂(t.t[i], q.d[:,i])
        end

        (p1_error, p2_error)
    end

end
