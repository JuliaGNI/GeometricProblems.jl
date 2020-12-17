
using Parameters

export ϑ, ω, hamiltonian


function ϑ(t::Number, q::AbstractVector, Θ::AbstractVector)
    Θ[1] = ϑ₁(t,q)
    Θ[2] = ϑ₂(t,q)
    nothing
end

function ϑ(t::Number, q::AbstractVector)
    [ϑ₁(t,q), ϑ₂(t,q)]
end

function ϑ(t::Number, q::AbstractVector, k::Int)
    if k == 1
        ϑ₁(t, q)
    elseif k == 2
        ϑ₂(t, q)
    else
        nothing
    end
end
    
ϑ(t::Number, q::AbstractVector, params::NamedTuple, k::Int) = ϑ(t,q,k)

function lotka_volterra_2d_pᵢ(qᵢ, tᵢ=0)
    pᵢ = zero(qᵢ)

    if ndims(qᵢ) == 1
        ϑ(tᵢ, qᵢ, pᵢ)
    else
        for i in 1:size(qᵢ,2)
            ϑ(tᵢ, (@view qᵢ[:,i]), (@view pᵢ[:,i]))
        end
    end
    pᵢ
end

function ω(t, q, Ω)
    Ω[1,1] = 0
    Ω[1,2] = dϑ₁dx₂(t,q) - dϑ₂dx₁(t,q)

    Ω[2,1] = dϑ₂dx₁(t,q) - dϑ₁dx₂(t,q)
    Ω[2,2] = 0

    nothing
end


function hamiltonian(t, q, params)
    @unpack a₁, a₂, b₁, b₂ = params
    a₁*q[1] + a₂*q[2] + b₁*log(q[1]) + b₂*log(q[2])
end

hamiltonian(t, q, p, params) = hamiltonian(t, q, params)

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
    + q[1] * (a₂*q[2] + b₂)
end

function v₂(t, q, params)
    @unpack a₁, a₂, b₁, b₂ = params
    - q[2] * (a₁*q[1] + b₁)
end


f₁(t, q, v) = dϑ₁dx₁(t,q) * v[1] + dϑ₂dx₁(t,q) * v[2]
f₂(t, q, v) = dϑ₁dx₂(t,q) * v[1] + dϑ₂dx₂(t,q) * v[2]

g₁(t, q, v) = dϑ₁dx₁(t,q) * v[1] + dϑ₁dx₂(t,q) * v[2]
g₂(t, q, v) = dϑ₂dx₁(t,q) * v[1] + dϑ₂dx₂(t,q) * v[2]


lotka_volterra_2d_ϑ(t, q, Θ, params) = ϑ(t, q, Θ)
lotka_volterra_2d_ϑ(t, q, v, Θ, params) = ϑ(t, q, Θ)
lotka_volterra_2d_ω(t, q, Ω, params) = ω(t, q, Ω)


function lotka_volterra_2d_dH(t, q, dH, params)
    dH[1] = dHd₁(t, q, params)
    dH[2] = dHd₂(t, q, params)
    nothing
end

function lotka_volterra_2d_v(t, q, v, params)
    v[1] = v₁(t, q, params)
    v[2] = v₂(t, q, params)
    nothing
end

function lotka_volterra_2d_v(t, q, p, v, params)
    lotka_volterra_2d_v(t, q, v, params)
end

function lotka_volterra_2d_v_ham(t, q, p, v, params)
    v .= 0
    nothing
end

function lotka_volterra_2d_v_dae(t, q, v, params)
    v[1] = v[3]
    v[2] = v[4]
    v[3] = 0
    v[4] = 0
    nothing
end

function lotka_volterra_2d_f(t::Real, q::Vector, v::Vector, f::Vector, params)
    f[1] = f₁(t,q,v) - dHd₁(t, q, params)
    f[2] = f₂(t,q,v) - dHd₂(t, q, params)
    nothing
end

function lotka_volterra_2d_f_ham(t::Real, q::Vector, f::Vector, params)
    f[1] = - dHd₁(t, q, params)
    f[2] = - dHd₂(t, q, params)
    nothing
end

function lotka_volterra_2d_f_ham(t::Real, q::Vector, v::Vector, f::Vector, params)
    lotka_volterra_2d_f_ham(t, q, f, params)
end

function lotka_volterra_2d_g(t::Real, q::Vector, v::Vector, g::Vector, params)
    g[1] = f₁(t,q,v)
    g[2] = f₂(t,q,v)
    nothing
end

function lotka_volterra_2d_g(t::Real, q::Vector, p::Vector, v::Vector, g::Vector, params)
    lotka_volterra_2d_g(t, q, v, g, params)
end

function lotka_volterra_2d_g̅(t::Real, q::Vector, v::Vector, g::Vector, params)
    g[1] = g₁(t,q,v)
    g[2] = g₂(t,q,v)
    nothing
end

function lotka_volterra_2d_g̅(t::Real, q::Vector, p::Vector, v::Vector, g::Vector, params)
    lotka_volterra_2d_g̅(t, q, v, g, params)
end

function lotka_volterra_2d_u_dae(t, q, λ, u, params)
    u[1] = 0
    u[2] = 0
    u[3] = λ[1]
    u[4] = λ[2]
    nothing
end

function lotka_volterra_2d_u(t, q, v, u, params)
    u .= v
    nothing
end

function lotka_volterra_2d_u(t, q, p, v, u, params)
    lotka_volterra_2d_u(t, q, v, u, params)
end

function lotka_volterra_2d_u̅(t, q, v, u, params)
    u[1] = v[1]
    u[2] = v[2]
    nothing
end

function lotka_volterra_2d_u̅(t, q, p, v, u, params)
    lotka_volterra_2d_u̅(t, q, v, u, params)
end

function lotka_volterra_2d_ϕ_dae(t, q, ϕ, params)
    ϕ[1] = q[3] - v₁(t,q,params)
    ϕ[2] = q[4] - v₂(t,q,params)
    nothing
end

function lotka_volterra_2d_ϕ(t, q, p, ϕ, params)
    ϕ[1] = p[1] - ϑ₁(t,q)
    ϕ[2] = p[2] - ϑ₂(t,q)
    nothing
end

function lotka_volterra_2d_ψ(t, q, p, v, f, ψ, params)
    ψ[1] = f[1] - g₁(t,q,v)
    ψ[2] = f[2] - g₂(t,q,v)
    nothing
end
