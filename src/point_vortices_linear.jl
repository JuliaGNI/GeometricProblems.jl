@doc raw"""
# Planar Point Vortices with linear one-form

"""
module PointVorticesLinear

    using GeometricEquations

    export odeproblem, iodeproblem, iodeproblem_dg, lodeproblem_formal_lagrangian,
           hamiltonian, angular_momentum, ϑ1, ϑ2, ϑ3, ϑ4,
           compute_energy, compute_energy_error, compute_angular_momentum_error,
           compute_momentum_error, compute_one_form

    const Δt = 0.01
    const nt = 1000
    const tspan = (0.0, Δt*nt)
    
    const γ₁ = 4.0
    const γ₂ = 2.0
    const d  = 1.0
    const q₀ = [γ₂*d/(γ₁+γ₂), 0.0, -γ₁*d/(γ₁+γ₂), 0.0]


    function hamiltonian(t,q)
        γ₁ * γ₂ * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (4π)
    end

    ϑ1(q) = - γ₁ * q[2] / 2
    ϑ2(q) = + γ₁ * q[1] / 2
    ϑ3(q) = - γ₂ * q[4] / 2
    ϑ4(q) = + γ₂ * q[3] / 2

    ϑ(q) = [ϑ1(q), ϑ2(q), ϑ3(q), ϑ4(q)]

    function angular_momentum(t,q)
        # γ₁ * (q[1]^2 + q[2]^2) * S(q[1],q[2]) +
        # γ₂ * (q[3]^2 + q[4]^2) * S(q[3],q[4])
        q[1] * ϑ2(q) - q[2] * ϑ1(q) +
        q[3] * ϑ4(q) - q[4] * ϑ3(q)
    end


    dϑ1d1(t, q) = zero(eltype(q))
    dϑ1d2(t, q) = - γ₁ / 2
    dϑ1d3(t, q) = zero(eltype(q))
    dϑ1d4(t, q) = zero(eltype(q))
    dϑ2d1(t, q) = + γ₁ / 2
    dϑ2d2(t, q) = zero(eltype(q))
    dϑ2d3(t, q) = zero(eltype(q))
    dϑ2d4(t, q) = zero(eltype(q))
    dϑ3d1(t, q) = zero(eltype(q))
    dϑ3d2(t, q) = zero(eltype(q))
    dϑ3d3(t, q) = zero(eltype(q))
    dϑ3d4(t, q) = - γ₂ / 2
    dϑ4d1(t, q) = zero(eltype(q))
    dϑ4d2(t, q) = zero(eltype(q))
    dϑ4d3(t, q) = + γ₂ / 2
    dϑ4d4(t, q) = zero(eltype(q))


    function ϑ(p, t, q)
        p[1] = ϑ1(q)
        p[2] = ϑ2(q)
        p[3] = ϑ3(q)
        p[4] = ϑ4(q)
    end

    function ω(Ω, t, q)
        Ω[1,1] = 0
        Ω[1,2] = dϑ1d2(t,q) - dϑ2d1(t,q)
        Ω[1,3] = dϑ1d3(t,q) - dϑ3d1(t,q)
        Ω[1,4] = dϑ1d4(t,q) - dϑ4d1(t,q)

        Ω[2,1] = dϑ2d1(t,q) - dϑ1d2(t,q)
        Ω[2,2] = 0
        Ω[2,3] = dϑ2d3(t,q) - dϑ3d2(t,q)
        Ω[2,4] = dϑ2d4(t,q) - dϑ4d2(t,q)

        Ω[3,1] = dϑ3d1(t,q) - dϑ1d3(t,q)
        Ω[3,2] = dϑ3d2(t,q) - dϑ2d3(t,q)
        Ω[3,3] = 0
        Ω[3,4] = dϑ3d4(t,q) - dϑ4d3(t,q)

        Ω[4,1] = dϑ4d1(t,q) - dϑ1d4(t,q)
        Ω[4,2] = dϑ4d2(t,q) - dϑ2d4(t,q)
        Ω[4,3] = dϑ4d3(t,q) - dϑ3d4(t,q)
        Ω[4,4] = 0

        nothing
    end


    f1(t, q, v) = + γ₁ * v[2] / 2
    f2(t, q, v) = - γ₁ * v[1] / 2
    f3(t, q, v) = + γ₂ * v[4] / 2
    f4(t, q, v) = - γ₂ * v[3] / 2


    dHd1(t, q) = + γ₁ * γ₂ * (q[1] - q[3]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π)
    dHd2(t, q) = + γ₁ * γ₂ * (q[2] - q[4]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π)
    dHd3(t, q) = - γ₁ * γ₂ * (q[1] - q[3]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π)
    dHd4(t, q) = - γ₁ * γ₂ * (q[2] - q[4]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π)

    function dH(dH, t, q, params)
        dH[1] = dHd1(t, q)
        dH[2] = dHd2(t, q)
        dH[3] = dHd3(t, q)
        dH[4] = dHd4(t, q)
        nothing
    end


    function point_vortices_v(v, t, q, params)
        v[1] = - dHd2(t,q) / γ₁
        v[2] = + dHd1(t,q) / γ₁
        v[3] = - dHd4(t,q) / γ₂
        v[4] = + dHd3(t,q) / γ₂
        nothing
    end

    function odeproblem(q₀=q₀; tspan = tspan, tstep = Δt)
        ODEProblem(point_vortices_v, tspan, tstep, q₀)
    end


    function point_vortices_ϑ(p, t, q, v, params)
        p[1] = ϑ1(q)
        p[2] = ϑ2(q)
        p[3] = ϑ3(q)
        p[4] = ϑ4(q)
        nothing
    end

    function point_vortices_f(f, t, q, v, params)
        f[1] = f1(t,q,v) - dHd1(t,q)
        f[2] = f2(t,q,v) - dHd2(t,q)
        f[3] = f3(t,q,v) - dHd3(t,q)
        f[4] = f4(t,q,v) - dHd4(t,q)
        nothing
    end

    function point_vortices_g(g, t, q, λ, params)
        g[1] = f1(t,q,λ)
        g[2] = f2(t,q,λ)
        g[3] = f3(t,q,λ)
        g[4] = f4(t,q,λ)
        nothing
    end

    point_vortices_g(g, t, q, p, λ, params) = point_vortices_g(g, t, q, λ, params)
    point_vortices_g(g, t, q, v, p, λ, params) = point_vortices_g(g, t, q, p, λ, params)

    function point_vortices_v(v, t, q, p, params)
        point_vortices_v(v, t, q, params)
    end

    function iodeproblem(q₀=q₀, p₀=ϑ(q₀); tspan = tspan, tstep = Δt)
        IODEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_g, tspan, tstep, q₀, p₀;
                    v̄=point_vortices_v)
    end

    function iodeproblem_dg(q₀=q₀; tspan = tspan, tstep = Δt)
        IODEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_g, tspan, tstep, q₀, q₀;
                    v=point_vortices_v)
    end

    function lodeproblem_formal_lagrangian(q₀=q₀, p₀=ϑ(q₀); tspan = tspan, tstep = Δt)
        LODEProblem(ϑ, point_vortices_f, point_vortices_g, tspan, tstep, q₀, p₀;
                    v̄=point_vortices_v, Ω=ω, ∇H=dH)
    end


    function compute_energy(t, q)
        h = zeros(q.nt+1)
        for i in 1:(q.nt+1)
            h[i] = hamiltonian(t.t[i], q.d[:,i])
        end
        return h
    end

    function compute_energy_error(t, q)
        h = zeros(q.nt+1)
        for i in 1:(q.nt+1)
            h[i] = hamiltonian(t.t[i], q.d[:,i])
        end
        h_error = (h .- h[1]) / h[1]
    end

    function compute_angular_momentum_error(t, q)
        P = zeros(q.nt+1)
        for i in 1:(q.nt+1)
            P[i] = angular_momentum(t.t[i], q.d[:,i])
        end
        P_error = (P .- P[1]) / P[1]
    end

    function compute_momentum_error(t, q, p)
        p1_error = zeros(q.nt+1)
        p2_error = zeros(q.nt+1)
        p3_error = zeros(q.nt+1)
        p4_error = zeros(q.nt+1)

        for i in 1:(q.nt+1)
            p1_error[i] = p.d[1,i] - ϑ1(q.d[:,i])
            p2_error[i] = p.d[2,i] - ϑ2(q.d[:,i])
            p3_error[i] = p.d[3,i] - ϑ3(q.d[:,i])
            p4_error[i] = p.d[4,i] - ϑ4(q.d[:,i])
        end

        (p1_error, p2_error, p3_error, p4_error)
    end

    function compute_one_form(t, q)
        p1 = zeros(q.nt+1)
        p2 = zeros(q.nt+1)
        p3 = zeros(q.nt+1)
        p4 = zeros(q.nt+1)

        for i in 1:(q.nt+1)
            p1[i] = ϑ1(q.d[:,i])
            p2[i] = ϑ2(q.d[:,i])
            p3[i] = ϑ3(q.d[:,i])
            p4[i] = ϑ4(q.d[:,i])
        end

        (p1, p2, p3, p4)
    end

end
