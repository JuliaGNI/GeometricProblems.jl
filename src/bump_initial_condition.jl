using ForwardDiff

"""
Third-degree spline that is used as a basis to construct the initial conditions. 
"""
function h(x::T) where T 
    if 0 ≤ x ≤ 1
        1 - 3 * x ^ 2 / 2 + 3 * x ^ 3 / 4 
    elseif T(1) < x ≤ 2
        (2 - x) ^ 3 / 4
    else
        T(0)
    end 
end

function s(ξ::Union{T, ForwardDiff.Dual}, μ::T) where T <: Real
    20 * μ * abs(ξ + μ / 2)
end

function s(ξ::AbstractVector{T}, μ::T) where T 
    s_closure(ξ_scal) = s(ξ_scal, μ)
    s_closure.(ξ)
end

function u₀(ξ::AbstractVector{T}, μ::T) where T 
    h.(s(ξ, μ))
end

u₀(ξ::Union{T, ForwardDiff.Dual}, μ::T) where T <: Real = h(s(ξ, μ))

function get_domain(N::Integer, T=Float64)
    Δx = 1. / (N - 1)
    T(-0.5) : Δx : T(0.5)
end

function get_p₀(ξ::T, μ::T) where T <: Real
    - μ * ForwardDiff.derivative(ξ -> u₀(ξ, μ), ξ)
end

function get_p₀(Ω::AbstractVector{T}, μ::T) where T 
    p₀_closure(ξ::T) = get_p₀(ξ, μ)
    p₀_closure.(Ω)
end

@doc raw"""
Produces initial conditions for the bump function. Here the ``p``-part is initialized with zeros. 
"""
function get_initial_condition(μ::T, N::Integer) where T 
    Ω = get_domain(N, T)
    (q =u₀(Ω, μ), p = zero(Ω))
end 

@doc raw"""
Produces initial condition for the bump function. Here the ``p``-part is initialized as ``-\mu\partial_\xi u_0(\xi, \mu)``. 
"""
function get_initial_condition2(μ::T, N::Integer) where T 
    Ω = get_domain(N, T)
    (q =u₀(Ω, μ), p = get_p₀(Ω, μ))
end 