module Diagnostics

using GeometricSolutions
using GeometricSolutions: compute_invariant, compute_invariant_error, compute_error_drift

export compute_one_form, compute_momentum_error
export compute_invariant, compute_invariant_error, compute_error_drift


"""
Compute the one-form (symplectic potential) for the solution of a Lagrangian system.

Arguments: `(t::TimeSeries, q::DataSeries, one_form::Base.Callable)`

The `one_form` function needs to take three arguments `(t,q,k)` where `k` is the index
of the one-form component.

Returns a DataSeries similar to `q` holding the time series of the one-form.
"""
function compute_one_form(t::TimeSeries, q::DataSeries, one_form::Base.Callable)
    ϑ = similar(q)
    try
        for i in axes(p, 2)
            for k in axes(p, 1)
                ϑ[k, i] = one_form(t[i], q[:, i], k)
            end
        end
    catch ex
        if isa(ex, DomainError)
            @warn("DOMAIN ERROR: One-form diagnostics crashed.")
        else
            throw(ex)
        end
    end
    return ϑ
end

"""
Computes the difference of the momentum and the one-form of an implicit ODE or DAE system.

Arguments: `(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, one_form::Function)`

The `one_form` function needs to take three arguments `(t,q,k)` where `k` is the index
of the one-form component.

Returns a DataSeries similar to `p` holding the time series of the difference between the momentum and the one-form.
"""

function compute_momentum_error(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, one_form::Function) where {T}
    e = similar(p)

    for n in axes(p, 1)
        for k in eachindex(p[n])
            e[n][k] = p[n][k] - one_form(t[n], q[n], k)
        end
    end

    return e
end

"""
Computes the difference of the momentum and the one-form of an implicit ODE or DAE system.

Arguments: `(p::DataSeries{DT}, ϑ::DataSeries{DT})`

Returns a DataSeries similar to `p` holding the time series of the difference between `p` and `ϑ`.
"""
function compute_momentum_error(p::DataSeries{DT}, ϑ::DataSeries{DT}) where {DT}
    @assert axes(p) == axes(ϑ)

    e = similar(p)
    parent(e) .= parent(p) .- parent(ϑ)

    return e
end

end
