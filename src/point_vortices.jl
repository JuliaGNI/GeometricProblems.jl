@doc raw"""
# Planar Point Vortices

"""
module PointVortices

    using GeometricEquations
    using GeometricSolutions

    export odeproblem, iodeproblem, idaeproblem,
           iodeproblem_dg, lodeproblem_formal_lagrangian,
           hamiltonian, angular_momentum, ϑ1, ϑ2, ϑ3, ϑ4
    
    export compute_energy_error, compute_angular_momentum_error

    const Δt = 0.01
    const nt = 1000
    const timespan = (0.0, Δt*nt)
    
    const reference_solution = [0.18722529318641928, 0.38967432450068706, 0.38125332930294187, 0.4258020604293123]

    const γ₁ = +0.5
    const γ₂ = +0.5
    const X0 = +0.5
    const Y0 = +0.1
    const X1 = +0.5
    const Y1 = -0.1


    S(x::T, y::T) where {T} = 1 + x^2 + y^2
    S2(x::T, y::T) where {T} = 1 + 2x^2 + 2y^2
    
    dSdx(x::T, y::T) where {T} = 2x
    dSdy(x::T, y::T) where {T} = 2y


    function hamiltonian(t, q, params)
        γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π)
    end

    ϑ1(q) = - γ₁ * q[2] * S(q[1], q[2]) / 2
    ϑ2(q) = + γ₁ * q[1] * S(q[1], q[2]) / 2
    ϑ3(q) = - γ₂ * q[4] * S(q[3], q[4]) / 2
    ϑ4(q) = + γ₂ * q[3] * S(q[3], q[4]) / 2

    ϑ(q) = [ϑ1(q), ϑ2(q), ϑ3(q), ϑ4(q)]

    function angular_momentum(t, q, params)
        q[1] * ϑ2(q) - q[2] * ϑ1(q) +
        q[3] * ϑ4(q) - q[4] * ϑ3(q)
    end


    const q₀ = [X0, Y0, X1, Y1]
    const p₀ = ϑ(q₀)


    dϑ1d1(t, q) = - γ₁ * q[2] * dSdx(q[1], q[2]) / 2
    dϑ1d2(t, q) = - γ₁ * q[2] * dSdy(q[1], q[2]) / 2 - γ₁ * S(q[1], q[2]) / 2
    dϑ1d3(t, q) = zero(eltype(q))
    dϑ1d4(t, q) = zero(eltype(q))
    dϑ2d1(t, q) = + γ₁ * q[1] * dSdx(q[1], q[2]) / 2 + γ₁ * S(q[1], q[2]) / 2
    dϑ2d2(t, q) = + γ₁ * q[1] * dSdy(q[1], q[2]) / 2
    dϑ2d3(t, q) = zero(eltype(q))
    dϑ2d4(t, q) = zero(eltype(q))
    dϑ3d1(t, q) = zero(eltype(q))
    dϑ3d2(t, q) = zero(eltype(q))
    dϑ3d3(t, q) = - γ₂ * q[4] * dSdx(q[3], q[4]) / 2
    dϑ3d4(t, q) = - γ₂ * q[4] * dSdy(q[3], q[4]) / 2 - γ₂ * S(q[3], q[4]) / 2
    dϑ4d1(t, q) = zero(eltype(q))
    dϑ4d2(t, q) = zero(eltype(q))
    dϑ4d3(t, q) = + γ₂ * q[3] * dSdx(q[3], q[4]) / 2 + γ₂ * S(q[3], q[4]) / 2
    dϑ4d4(t, q) = + γ₂ * q[3] * dSdy(q[3], q[4]) / 2



    function ϑ(p, t, q, params)
        p[1] = ϑ1(q)
        p[2] = ϑ2(q)
        p[3] = ϑ3(q)
        p[4] = ϑ4(q)
    end

    function ω(Ω, t, q, params)
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


    function f1(t, q, v)
        γ₁ * ( dSdx(q[1],q[2]) * (q[1] * v[2] - q[2] * v[1]) + v[2] * S(q[1], q[2]) ) / 2
    end

    function f2(t, q, v)
        γ₁ * ( dSdy(q[1],q[2]) * (q[1] * v[2] - q[2] * v[1]) - v[1] * S(q[1], q[2]) ) / 2
    end

    function f3(t, q, v)
        γ₂ * ( dSdx(q[3],q[4]) * (q[3] * v[4] - q[4] * v[3]) + v[4] * S(q[3], q[4]) ) / 2
    end

    function f4(t, q, v)
        γ₂ * ( dSdy(q[3],q[4]) * (q[3] * v[4] - q[4] * v[3]) - v[3] * S(q[3], q[4]) ) / 2
    end


    function dHd1(t, q)
        + γ₁ * γ₂ * dSdx(q[1],q[2]) * S(q[3],q[4]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) +
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[1] - q[3]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end

    function dHd2(t, q)
        + γ₁ * γ₂ * dSdy(q[1],q[2]) * S(q[3],q[4]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) +
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[2] - q[4]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end

    function dHd3(t, q)
        + γ₁ * γ₂ * dSdx(q[3],q[4]) * S(q[1],q[2]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) -
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[1] - q[3]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end

    function dHd4(t, q)
        + γ₁ * γ₂ * dSdy(q[3],q[4]) * S(q[1],q[2]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π) -
          γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * (q[2] - q[4]) / ( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / π
    end

    function dH(dH, t, q, params)
        dH[1] = dHd1(t, q)
        dH[2] = dHd2(t, q)
        dH[3] = dHd3(t, q)
        dH[4] = dHd4(t, q)
        nothing
    end


    function point_vortices_v(v, t, q, params)
        denominator1 = 1 / (γ₁ * S2(q[1], q[2]))
        denominator2 = 1 / (γ₂ * S2(q[3], q[4]))

        v[1] = - dHd2(t,q) * denominator1
        v[2] = + dHd1(t,q) * denominator1
        v[3] = - dHd4(t,q) * denominator2
        v[4] = + dHd3(t,q) * denominator2

        nothing
    end

    function point_vortices_v(v, t, q, p, params)
        point_vortices_v(v, t, q, params)
    end


    function point_vortices_ϑ(p, t, q, params)
        p[1] = ϑ1(q)
        p[2] = ϑ2(q)
        p[3] = ϑ3(q)
        p[4] = ϑ4(q)
        nothing
    end

    function point_vortices_ϑ(p, t, q, v, params)
        point_vortices_ϑ(p, t, q, params)
    end

    function point_vortices_f(f, t, q, v, params)
        f[1] = f1(t,q,v) - dHd1(t,q)
        f[2] = f2(t,q,v) - dHd2(t,q)
        f[3] = f3(t,q,v) - dHd3(t,q)
        f[4] = f4(t,q,v) - dHd4(t,q)
        nothing
    end

    point_vortices_f(f, t, q, v, p, params) = point_vortices_f(f, t, q, v, params)

    function point_vortices_g(g, t, q, λ, params)
        g[1] = f1(t,q,λ)
        g[2] = f2(t,q,λ)
        g[3] = f3(t,q,λ)
        g[4] = f4(t,q,λ)
        nothing
    end

    point_vortices_g(g, t, q, p, λ, params) = point_vortices_g(g, t, q, λ, params)
    point_vortices_g(g, t, q, v, p, λ, params) = point_vortices_g(g, t, q, p, λ, params)

    function point_vortices_u(u, t, q, v, params)
        u[1] = v[1]
        u[2] = v[2]
        u[3] = v[3]
        u[4] = v[4]
        nothing
    end

    point_vortices_u(u, t, q, p, λ, params) = point_vortices_u(u, t, q, λ, params)
    point_vortices_u(u, t, q, v, p, λ, params) = point_vortices_u(u, t, q, p, λ, params)

    function point_vortices_ϕ(ϕ, t, q, p, params)
        ϕ[1] = p[1] - ϑ1(q)
        ϕ[2] = p[2] - ϑ2(q)
        ϕ[3] = p[3] - ϑ3(q)
        ϕ[4] = p[4] - ϑ4(q)
        nothing
    end

    point_vortices_ϕ(ϕ, t, q, v, p, params) = point_vortices_ϕ(ϕ, t, q, p, params)


    function odeproblem(q₀=q₀; timespan = timespan, timestep = Δt)
        ODEProblem(point_vortices_v, timespan, timestep, q₀)
    end

    function iodeproblem(q₀=q₀, p₀=ϑ(q₀); timespan = timespan, timestep = Δt)
        IODEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_g, timespan, timestep, q₀, p₀;
                    v̄=point_vortices_v)
    end

    function idaeproblem(q₀=q₀, p₀=ϑ(q₀), λ₀=zero(q₀); timespan = timespan, timestep = Δt)
        IDAEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_u, point_vortices_g,
                    point_vortices_ϕ, timespan, timestep, q₀, p₀, λ₀;
                    v̄=point_vortices_v)
    end

    function idoeproblem_dg(q₀=q₀; timespan = timespan, timestep = Δt)
        IODEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_g, timespan, timestep, q₀, q₀;
                    v=point_vortices_v)
    end

    function lodeproblem_formal_lagrangian(q₀=q₀, p₀=ϑ(q₀); timespan = timespan, timestep = Δt)
        LODEProblem(ϑ, point_vortices_f, point_vortices_g, timespan, timestep, q₀, p₀;
                    v̄=point_vortices_v, Ω=ω, ∇H=dH)
    end


    function compute_energy_error(t, q::DataSeries{T}, params) where {T}
        h = DataSeries(T, q.nt)
        e = DataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i], params)
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_angular_momentum_error(t, q::DataSeries{T}, params) where {T}
        m = DataSeries(T, q.nt)
        e = DataSeries(T, q.nt)

        for i in axes(q,2)
            m[i] = angular_momentum(t[i], q[:,i], params)
            e[i] = (m[i] - m[0]) / m[0]
        end

        (m, e)
    end

end
