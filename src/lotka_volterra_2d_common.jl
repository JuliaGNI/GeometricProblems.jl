
using Parameters

import NaNMath: log

export ϑ, ω, hamiltonian


function ϑ(Θ::AbstractVector, t, q::AbstractVector)
    Θ[1] = ϑ₁(t, q)
    Θ[2] = ϑ₂(t, q)
    nothing
end

function ϑ(t, q::AbstractVector)
    [ϑ₁(t, q), ϑ₂(t, q)]
end

function ϑ(t, q::AbstractVector, k::Int)
    if k == 1
        ϑ₁(t, q)
    elseif k == 2
        ϑ₂(t, q)
    else
        throw(BoundsError(ϑ, k))
    end
end

ϑ(t, q::AbstractVector, params::NamedTuple, k::Int) = ϑ(t, q, k)

function lotka_volterra_2d_pᵢ(qᵢ, tᵢ=0)
    pᵢ = zero(qᵢ)

    if ndims(qᵢ) == 1
        ϑ(pᵢ, tᵢ, qᵢ)
    else
        for i in axes(qᵢ, 2)
            ϑ((@view pᵢ[:, i]), tᵢ, (@view qᵢ[:, i]))
        end
    end
    pᵢ
end

function ω(Ω, t, q)
    Ω[1, 1] = 0
    Ω[1, 2] = dϑ₁dx₂(t, q) - dϑ₂dx₁(t, q)

    Ω[2, 1] = dϑ₂dx₁(t, q) - dϑ₁dx₂(t, q)
    Ω[2, 2] = 0

    nothing
end


function hamiltonian(t, q, params)
    @unpack a₁, a₂, b₁, b₂ = params
    a₁ * q[1] + a₂ * q[2] + b₁ * log(q[1]) + b₂ * log(q[2])
end

hamiltonian_iode(t, q, params) = hamiltonian(t, q, params) # This is a workaround. It should be removed asap.
hamiltonian_iode(t, q, v, params) = hamiltonian(t, q, params)

hamiltonian_pode(t, q, p, params) = hamiltonian(t, q, params)

function lagrangian(t, q, v, params)
    ϑ₁(t, q) * v[1] + ϑ₂(t, q) * v[2] - hamiltonian(t, q, params)
end


function dHd₁(t, q, params)
    @unpack a₁, b₁ = params
    a₁ + b₁ / q[1]
end

function dHd₂(t, q, params)
    @unpack a₂, b₂ = params
    a₂ + b₂ / q[2]
end


function v₁(t, q, params)
    @unpack a₁, a₂, b₁, b₂ = params
    +q[1] * (a₂ * q[2] + b₂)
end

function v₂(t, q, params)
    @unpack a₁, a₂, b₁, b₂ = params
    -q[2] * (a₁ * q[1] + b₁)
end


f₁(t, q, v) = dϑ₁dx₁(t, q) * v[1] + dϑ₂dx₁(t, q) * v[2]
f₂(t, q, v) = dϑ₁dx₂(t, q) * v[1] + dϑ₂dx₂(t, q) * v[2]

g₁(t, q, v) = dϑ₁dx₁(t, q) * v[1] + dϑ₁dx₂(t, q) * v[2]
g₂(t, q, v) = dϑ₂dx₁(t, q) * v[1] + dϑ₂dx₂(t, q) * v[2]


lotka_volterra_2d_ϑ(Θ, t, q, params) = ϑ(Θ, t, q)
lotka_volterra_2d_ϑ(Θ, t, q, v, params) = ϑ(Θ, t, q)
lotka_volterra_2d_ω(Ω, t, q, params) = ω(Ω, t, q)
lotka_volterra_2d_ω(Ω, t, q, v, params) = ω(Ω, t, q)


function lotka_volterra_2d_dH(dH, t, q, params)
    dH[1] = dHd₁(t, q, params)
    dH[2] = dHd₂(t, q, params)
    nothing
end

function lotka_volterra_2d_v(v, t, q, params)
    v[1] = v₁(t, q, params)
    v[2] = v₂(t, q, params)
    nothing
end

function lotka_volterra_2d_v(v, t, q, p, params)
    lotka_volterra_2d_v(v, t, q, params)
end

function lotka_volterra_2d_v_ham(v, t, q, p, params)
    v .= 0
    nothing
end

function lotka_volterra_2d_v_dae(v, t, q, params)
    v[1] = v[3]
    v[2] = v[4]
    v[3] = 0
    v[4] = 0
    nothing
end

function lotka_volterra_2d_f(f::AbstractVector, t, q::AbstractVector, v::AbstractVector, params)
    f[1] = f₁(t, q, v) - dHd₁(t, q, params)
    f[2] = f₂(t, q, v) - dHd₂(t, q, params)
    nothing
end

function lotka_volterra_2d_f_ham(f::AbstractVector, t, q::AbstractVector, params)
    f[1] = -dHd₁(t, q, params)
    f[2] = -dHd₂(t, q, params)
    nothing
end

function lotka_volterra_2d_f_ham(f::AbstractVector, t, q::AbstractVector, v::AbstractVector, params)
    lotka_volterra_2d_f_ham(f, t, q, params)
end

function lotka_volterra_2d_g(g::AbstractVector, t, q::AbstractVector, v::AbstractVector, params)
    g[1] = f₁(t, q, v)
    g[2] = f₂(t, q, v)
    nothing
end

lotka_volterra_2d_g(g, t, q, p, λ, params) = lotka_volterra_2d_g(g, t, q, λ, params)
lotka_volterra_2d_g(g, t, q, v, p, λ, params) = lotka_volterra_2d_g(g, t, q, p, λ, params)

function lotka_volterra_2d_ḡ(g::AbstractVector, t, q::AbstractVector, v::AbstractVector, params)
    g[1] = g₁(t, q, v)
    g[2] = g₂(t, q, v)
    nothing
end

lotka_volterra_2d_ḡ(g, t, q, p, λ, params) = lotka_volterra_2d_ḡ(g, t, q, λ, params)
lotka_volterra_2d_ḡ(g, t, q, v, p, λ, params) = lotka_volterra_2d_ḡ(g, t, q, p, λ, params)

function lotka_volterra_2d_u_dae(u, t, q, λ, params)
    u[1] = 0
    u[2] = 0
    u[3] = λ[1]
    u[4] = λ[2]
    nothing
end

function lotka_volterra_2d_u(u, t, q, λ, params)
    u .= λ
    nothing
end

lotka_volterra_2d_u(u, t, q, p, λ, params) = lotka_volterra_2d_u(u, t, q, λ, params)
lotka_volterra_2d_u(u, t, q, v, p, λ, params) = lotka_volterra_2d_u(u, t, q, p, λ, params)

function lotka_volterra_2d_ū(u, t, q, λ, params)
    u .= λ
    nothing
end

lotka_volterra_2d_ū(u, t, q, p, λ, params) = lotka_volterra_2d_ū(u, t, q, λ, params)
lotka_volterra_2d_ū(u, t, q, v, p, λ, params) = lotka_volterra_2d_ū(u, t, q, p, λ, params)

function lotka_volterra_2d_ϕ_dae(ϕ, t, q, params)
    ϕ[1] = q[3] - v₁(t, q, params)
    ϕ[2] = q[4] - v₂(t, q, params)
    nothing
end

function lotka_volterra_2d_ϕ(ϕ, t, q, p, params)
    ϕ[1] = p[1] - ϑ₁(t, q)
    ϕ[2] = p[2] - ϑ₂(t, q)
    nothing
end

lotka_volterra_2d_ϕ(ϕ, t, q, v, p, params) = lotka_volterra_2d_ϕ(ϕ, t, q, p, params)

function lotka_volterra_2d_ψ(ψ, t, q, p, q̇, ṗ, params)
    ψ[1] = ṗ[1] - g₁(t, q, q̇)
    ψ[2] = ṗ[2] - g₂(t, q, q̇)
    nothing
end

lotka_volterra_2d_ψ(ψ, t, q, v, p, q̇, ṗ, params) = lotka_volterra_2d_ψ(ψ, t, q, p, q̇, ṗ, params)

function lotka_volterra_2d_ψ_lode(ψ, t, q, v, p, q̇, ṗ, params)
    ψ[1] = f₁(t, q, v) - g₁(t, q, v) - dHd₁(t, q, params)
    ψ[2] = f₂(t, q, v) - g₂(t, q, v) - dHd₂(t, q, params)
    nothing
end
