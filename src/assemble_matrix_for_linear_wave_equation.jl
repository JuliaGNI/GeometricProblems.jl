using OffsetArrays

function assemble_matrix(μ::T, Δx::T, Ñ::Integer) where T 
    K = zeros(T, Ñ + 2, Ñ + 2)
    K = OffsetArray(K, OffsetArrays.Origin(0, 0))
    fac = μ ^ 2 / (2 * Δx)
    for i in 1:Ñ
        K[i,i] += fac 
        K[i,i+1] -= fac 
        K[i,i-1] -= fac 
        K[i-1,i-1] += fac/T(2)
        K[i+1,i+1] += fac/T(2)
    end

    K 
end

function assemble_matrix(μ::T, Ñ::Integer) where T
    assemble_matrix(μ, one(T) / (Ñ - 1), Ñ)
end