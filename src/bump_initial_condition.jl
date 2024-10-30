"""
Third-degree spline that is used as a basis to construct the initial conditions. 
"""
function h(x::T) where T <: Real
    if 0 ≤ x ≤ 1
        1 - 3 * x ^ 2 / 2 + 3 * x ^ 3 / 4 +2
    elseif 1 < x ≤ 2
        (2 - x) ^ 3 / 4 +2
    else
        zero(T) +2 
    end 
end

function ∂h(x::T) where T <: Real
    if 0 ≤ x ≤ 1
        -3x + 9 * x ^ 2 / 4
    elseif 1 < x ≤ 2
        - 3 * (2 - x) ^ 2 / 4
    else
        zero(T)
    end
end

function s(ξ::Real, μ::Real)
    20μ * abs(ξ + μ / 2)
end

function ∂s(ξ::T, μ::T) where T <: Real
    ξ + μ / 2 ≥ 0 ? 20μ : -20μ
end

function s(ξ::AbstractVector{T}, μ::T) where T <: Real
    s_closure(ξ_scalar) = s(ξ_scalar, μ)
    s_closure.(ξ)
end

function ∂s(ξ::AbstractVector{T}, μ::T) where T <: Real
    ∂s_closure(ξ_scalar) = s(ξ_scalar, μ)
    ∂s_closure.(ξ)
end

u₀(ξ::Real, μ::Real) = h(s(ξ, μ))

function u₀(ξ::AbstractVector{T}, μ::T) where T <: Real
    h.(s(ξ, μ))
end

function ∂u₀(ξ::T, μ::T) where T <: Real
    ∂h(s(ξ, μ)) * ∂s(ξ, μ)
end

function ∂u₀(ξ::AbstractVector{T}, μ::T) where T <: Real 
    ∂u₀_closure(ξ_scalar) = ∂u₀(ξ_scalar, μ)
    ∂u₀_closure.(ξ)
end

function compute_domain(N::Integer, T=Float64)
    Δx = 1. / (N - 1)
    T(-0.5) : Δx : T(0.5)
end

function compute_p₀(ξ::T, μ::T) where T <: Real
    - μ * ∂u₀(ξ, μ)
end

function compute_p₀(Ω::AbstractVector{T}, μ::T) where T 
    p₀_closure(ξ::T) = compute_p₀(ξ, μ)
    p₀_closure.(Ω)
end

@doc raw"""
Produces initial conditions for the bump function. Here the ``p``-part is initialized with zeros. 
"""
function compute_initial_condition(μ::T, N::Integer) where T 
    Ω = compute_domain(N, T)
    (q =u₀(Ω, μ), p = zero(Ω))
end 

@doc raw"""
Produces initial condition for the bump function. Here the ``p``-part is initialized as ``-\mu\partial_\xi u_0(\xi, \mu)``. 
"""
function compute_initial_condition2(μ::T, N::Integer) where T 
    Ω = compute_domain(N, T)
    (q =u₀(Ω, μ), p = compute_p₀(Ω, μ))
end 

function compute_initial_condition3(μ::T, N::Integer) where T 
    Ω = compute_domain(N, T)
    (q =cos.(2pi*Ω), p = 2pi * sin.(2pi*Ω))
end 