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

function s(ξ::T, μ::T) where T 
    28 * μ * abs(ξ + μ / 2)
end

function s(ξ::AbstractVector{T}, μ::T) where T 
    s_closure(ξ_scal) = s(ξ_scal, μ)
    s_closure.(ξ)
end

function u₀(ξ::AbstractVector{T}, μ::T) where T 
    h.(s(ξ, μ))
end

u₀(ξ, μ) = h(s(ξ, μ))

function get_initial_condition(μ::T, Δx::T) where T 
    Ω = T(-.5):Δx:T(0.5)
    (q=u₀(Ω,μ), p=zero(Ω))
end 

"""
If you call this function with a Float and an Integer, then the integer will be interpreted as the number of nodes, i.e. degrees of freedom.
"""
function get_initial_condition(μ::T, Ñ::Integer) where T 
    get_initial_condition(μ, T( 1 / (Ñ - 1)))
end