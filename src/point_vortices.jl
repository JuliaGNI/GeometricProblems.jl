@doc raw"""
# Planar Point Vortices

"""
module PointVortices

    using GeometricEquations
    using GeometricSolutions

    export point_vortices_ode, point_vortices_iode, point_vortices_idae,
           point_vortices_dg, point_vortices_formal_lagrangian,
           hamiltonian, angular_momentum, ϑ1, ϑ2, ϑ3, ϑ4
    
    export compute_energy_error, compute_angular_momentum_error

    const Δt = 0.01
    const nt = 1000
    const tspan = (0.0, Δt*nt)
    
    const reference_solution = [0.18722529318641928, 0.38967432450068706, 0.38125332930294187, 0.4258020604293123]

    const γ₁ = +0.5
    const γ₂ = +0.5
    const X0 = +0.5
    const Y0 = +0.1
    const X1 = +0.5
    const Y1 = -0.1


    function S(x::T, y::T) where {T}
        1 + x^2 + y^2
    end

    function S2(x::T, y::T) where {T}
        1 + 2x^2 + 2y^2
    end

    function dSdx(x::T, y::T) where {T}
        2x
    end

    function dSdy(x::T, y::T) where {T}
        2y
    end


    function hamiltonian(t, q, params)
        γ₁ * γ₂ * S(q[1],q[2]) * S(q[3],q[4]) * log( (q[1] - q[3])^2 + (q[2] - q[4])^2 ) / (2π)
    end

    function angular_momentum(t, q, params)
        q[1] * ϑ2(t, q) - q[2] * ϑ1(t, q) +
        q[3] * ϑ4(t, q) - q[4] * ϑ3(t, q)
    end


    function ϑ1(t, q)
        - γ₁ * q[2] * S(q[1], q[2]) / 2
    end

    function ϑ2(t, q)
        + γ₁ * q[1] * S(q[1], q[2]) / 2
    end

    function ϑ3(t, q)
        - γ₂ * q[4] * S(q[3], q[4]) / 2
    end

    function ϑ4(t, q)
        + γ₂ * q[3] * S(q[3], q[4]) / 2
    end


    const q₀ = [X0, Y0, X1, Y1]
    const p₀ = [ϑ1(0, q₀), ϑ2(0, q₀), ϑ3(0, q₀), ϑ4(0, q₀)]


    function dϑ1d1(t, q)
        - γ₁ * q[2] * dSdx(q[1], q[2]) / 2
    end

    function dϑ1d2(t, q)
        - γ₁ * q[2] * dSdy(q[1], q[2]) / 2 - γ₁ * S(q[1], q[2]) / 2
    end

    function dϑ1d3(t, q)
        zero(eltype(q))
    end

    function dϑ1d4(t, q)
        zero(eltype(q))
    end

    function dϑ2d1(t, q)
        + γ₁ * q[1] * dSdx(q[1], q[2]) / 2 + γ₁ * S(q[1], q[2]) / 2
    end

    function dϑ2d2(t, q)
        + γ₁ * q[1] * dSdy(q[1], q[2]) / 2
    end

    function dϑ2d3(t, q)
        zero(eltype(q))
    end

    function dϑ2d4(t, q)
        zero(eltype(q))
    end

    function dϑ3d1(t, q)
        zero(eltype(q))
    end

    function dϑ3d2(t, q)
        zero(eltype(q))
    end

    function dϑ3d3(t, q)
        - γ₂ * q[4] * dSdx(q[3], q[4]) / 2
    end

    function dϑ3d4(t, q)
        - γ₂ * q[4] * dSdy(q[3], q[4]) / 2 - γ₂ * S(q[3], q[4]) / 2
    end

    function dϑ4d1(t, q)
        zero(eltype(q))
    end

    function dϑ4d2(t, q)
        zero(eltype(q))
    end

    function dϑ4d3(t, q)
        + γ₂ * q[3] * dSdx(q[3], q[4]) / 2 + γ₂ * S(q[3], q[4]) / 2
    end

    function dϑ4d4(t, q)
        + γ₂ * q[3] * dSdy(q[3], q[4]) / 2
    end


    function ϑ(t, q, p, params)
        p[1] = ϑ1(t, q)
        p[2] = ϑ2(t, q)
        p[3] = ϑ3(t, q)
        p[4] = ϑ4(t, q)
    end

    function ω(t, q, Ω, params)
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

    function dH(t, q, dH, params)
        dH[1] = dHd1(t, q)
        dH[2] = dHd2(t, q)
        dH[3] = dHd3(t, q)
        dH[4] = dHd4(t, q)
        nothing
    end


    function point_vortices_p₀(q₀, t₀, params)
        p₀ = zero(q₀)
        tq = zeros(eltype(q₀), size(q₀,1))
        tp = zeros(eltype(p₀), size(p₀,1))

        if ndims(q₀) == 1
            ϑ(t₀, q₀, p₀, params)
        else
            for i in 1:size(q₀,2)
                simd_copy_xy_first!(tq, q₀, i)
                ϑ(t₀, tq, tp, params)
                simd_copy_yx_first!(tp, p₀, i)
            end
        end
        p₀
    end


    function point_vortices_v(t, q, v, params)
        denominator1 = 1 / (γ₁ * S2(q[1], q[2]))
        denominator2 = 1 / (γ₂ * S2(q[3], q[4]))

        v[1] = - dHd2(t,q) * denominator1
        v[2] = + dHd1(t,q) * denominator1
        v[3] = - dHd4(t,q) * denominator2
        v[4] = + dHd3(t,q) * denominator2

        nothing
    end

    function point_vortices_v(t, q, p, v, params)
        point_vortices_v(t, q, v, params)
    end


    function point_vortices_ϑ(t, q, p, params)
        p[1] = ϑ1(t,q)
        p[2] = ϑ2(t,q)
        p[3] = ϑ3(t,q)
        p[4] = ϑ4(t,q)
        nothing
    end

    function point_vortices_ϑ(t, q, v, p, params)
        point_vortices_ϑ(t, q, p, params)
    end

    function point_vortices_f(t, q, v, f, params)
        f[1] = f1(t,q,v) - dHd1(t,q)
        f[2] = f2(t,q,v) - dHd2(t,q)
        f[3] = f3(t,q,v) - dHd3(t,q)
        f[4] = f4(t,q,v) - dHd4(t,q)
        nothing
    end

    function point_vortices_f(t, q, v, p, f, params)
        point_vortices_f(t, q, v, f, params)
    end

    function point_vortices_g(t, q, λ, g, params)
        g[1] = f1(t,q,λ)
        g[2] = f2(t,q,λ)
        g[3] = f3(t,q,λ)
        g[4] = f4(t,q,λ)
        nothing
    end

    function point_vortices_g(t, q, p, λ, g, params)
        point_vortices_g(t, q, λ, g, params)
    end

    function point_vortices_u(t, q, v, u, params)
        u[1] = v[1]
        u[2] = v[2]
        u[3] = v[3]
        u[4] = v[4]
        nothing
    end

    function point_vortices_u(t, q, p, v, u, params)
        point_vortices_u(t, q, v, u, params)
    end

    function point_vortices_ϕ(t, q, p, ϕ, params)
        ϕ[1] = p[1] - ϑ1(t,q)
        ϕ[2] = p[2] - ϑ2(t,q)
        ϕ[3] = p[3] - ϑ3(t,q)
        ϕ[4] = p[4] - ϑ4(t,q)
        nothing
    end


    function point_vortices_ode(q₀=q₀; tspan = tspan, tstep = Δt)
        ODEProblem(point_vortices_v, tspan, tstep, q₀)
    end

    function point_vortices_iode(q₀=q₀, p₀=point_vortices_p₀(q₀, 0.0, NullParameters()); tspan = tspan, tstep = Δt)
        IODEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_g, tspan, tstep, q₀, p₀;
                    v̄=point_vortices_v)
    end

    function point_vortices_idae(q₀=q₀, p₀=point_vortices_p₀(q₀, 0.0, NullParameters()), λ₀=zero(q₀); tspan = tspan, tstep = Δt)
        IDAEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_u, point_vortices_g,
                    point_vortices_ϕ, tspan, tstep, q₀, p₀, λ₀;
                    v̄=point_vortices_v)
    end

    function point_vortices_dg(q₀=q₀; tspan = tspan, tstep = Δt)
        IODEProblem(point_vortices_ϑ, point_vortices_f,
                    point_vortices_g, tspan, tstep, q₀, q₀;
                    v=point_vortices_v)
    end

    function point_vortices_formal_lagrangian(q₀=q₀; tspan = tspan, tstep = Δt)
        p₀ = point_vortices_p₀(q₀, 0.0, NullParameters())
        LODEProblem(ϑ, point_vortices_f, point_vortices_g, tspan, tstep, q₀, p₀;
                    v̄=point_vortices_v, Ω=ω, ∇H=dH)
    end


    function compute_energy_error(t, q::AbstractDataSeries{T}, params) where {T}
        h = DataSeries(T, q.nt)
        e = DataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i], params)
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_angular_momentum_error(t, q::AbstractDataSeries{T}, params) where {T}
        m = DataSeries(T, q.nt)
        e = DataSeries(T, q.nt)

        for i in axes(q,2)
            m[i] = angular_momentum(t[i], q[:,i], params)
            e[i] = (m[i] - m[0]) / m[0]
        end

        (m, e)
    end


end
