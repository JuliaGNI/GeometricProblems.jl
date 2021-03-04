module Diagnostics

    using GeometricIntegrators.Common
    using GeometricIntegrators.Solutions

    export compute_one_form, compute_invariant, compute_invariant_error, compute_momentum_error, compute_error_drift

    """
    Takes a 1d DataSeries holding an invariant and computes the relative error `(inv(t)-inv(0))/inv(0)`.

    Returns a 1d DataSeries similar to the argument holding the time series of the relativ errors.
    """
    function compute_relative_error(invds::DataSeries{T,1}) where {T}
        errds = similar(invds)
        for i in eachindex(errds, invds)
            errds[i] = (invds[i] - invds[0]) / invds[0]
        end
        return errds
    end

    """
    Compute the one-form (symplectic potential) for the solution of a Lagrangian system.

    Arguments: `(t::TimeSeries, q::DataSeries, one_form::Function)`
    
    The `one_form` function needs to take three arguments `(t,q,k)` where `k` is the index
    of the one-form component.

    Returns a DataSeries similar to `q` holding the time series of the one-form.
    """
    function compute_one_form(t::TimeSeries, q::DataSeries, one_form::Function)
        ϑ = similar(q)
        try
            for i in axes(p,2)
                for k in axes(p,1)
                    ϑ[k,i] = one_form(t[i], q[:,i], k)
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
    Compute an invariant for the solution of an ODE or DAE system.

    Arguments: `(t::TimeSeries, q::DataSeries{T}, invariant::Function)`

    The `invariant` functions needs to take two arguments `(t,q)` and return the 
    corresponding value of the invariant.

    Returns a 1d DataSeries holding the time series of the invariant.
    """
    function compute_invariant(t::TimeSeries, q::DataSeries{<:AbstractArray{T}}, invariant::Function) where {T}
        invds = SDataSeries(T, ntime(q))
        try
            for i in eachindex(invds)
                invds[i] = invariant(t[i], q[i])
            end
        catch ex
            if isa(ex, DomainError)
                @warn("DOMAIN ERROR: Invariant diagnostics crashed.")
            else
                throw(ex)
            end
        end
        return invds
    end

    """
    Compute an invariant for the solution of a partitioned ODE or DAE system.

    Arguments: `(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, invariant::Function)`

    The `invariant` functions needs to take three arguments `(t,q,p)` and return the 
    corresponding value of the invariant.

    Returns a 1d DataSeries holding the time series of the invariant.
    """
    function compute_invariant(t::TimeSeries, q::DataSeries{<:AbstractArray{T}}, p::DataSeries{<:AbstractArray{T}}, invariant::Function) where {T}
        invds = SDataSeries(T, ntime(q))
        try
            for i in eachindex(invds)
                invds[i] = invariant(t[i], q[i], p[i])
            end
        catch ex
            if isa(ex, DomainError)
                @warn("DOMAIN ERROR: Invariant diagnostics crashed.")
            else
                throw(ex)
            end
        end
        return invds
    end

    """
    Compute the relative error of an invariant for the solution of an ODE or DAE system.

    Arguments: `(t::TimeSeries, q::DataSeries{T}, invariant::Function)`

    The `invariant` functions needs to take two arguments `(t,q)` and return the 
    corresponding value of the invariant.

    Returns a tuple of two 1d DataSeries holding the time series of the invariant and the relativ error, respectively.
    """
    function compute_invariant_error(t::TimeSeries, q::DataSeries{T}, invariant::Function) where {T}
        invds = compute_invariant(t, q, invariant)
        errds = compute_relative_error(invds)
        (invds, errds)
    end

    """
    Compute the relative error of an invariant for the solution of a partitioned ODE or DAE system.

    Arguments: `(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, invariant::Function)`

    The `invariant` functions needs to take three arguments `(t,q,p)` and return the 
    corresponding value of the invariant.

    Returns a tuple of two 1d DataSeries holding the time series of the invariant and the relativ error, respectively.
    """
    function compute_invariant_error(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, invariant::Function) where {T}
        invds = compute_invariant(t, q, p, invariant)
        errds = compute_relative_error(invds)
        (invds, errds)
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

        for i in axes(p,3)
            for n in axes(p,2)
                for k in axes(p,1)
                    e[k,n,i] = p[k,n,i] - one_form(t[n], q[:,n,i], k)
                end
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

    """
    Computes the drift in an invariant error.

    Arguments: `(t::TimeSeries, invariant_error::DataSeries{T,1}, interval_length=100)`

    The time series of the solution is split into intervals of `interval_length` time steps.
    In each interval, the maximum of the absolute value of the invariant error is determined.
    Returns a tuple of a TimeSeries that holds the centers of all intervals and a 1d DataSeries
    that holds the maxima.

    This is useful to detect drifts in invariants that are not preserved exactly but whose error
    is oscillating such as the energy error of Hamiltonian systems with symplectic integrators.
    """
    function compute_error_drift(t::TimeSeries, invariant_error::DataSeries{T,1}, interval_length=100) where {T}
        @assert ntime(t) == ntime(invariant_error)
        @assert mod(t.n, interval_length) == 0

        nint   = div(ntime(invariant_error), interval_length)
        Tdrift = TimeSeries(nint, timestep(t)*interval_length)
        Idrift = SDataSeries(T, nint)

        Tdrift[0] = t[0]

        for i in 1:nint
            i1 = interval_length*(i-1)+1
            i2 = interval_length*i
            Idrift[i] = maximum(abs.(invariant_error[i1:i2]))
            Tdrift[i] = div(t[i1] + t[i2], 2)
        end

        (Tdrift, Idrift)
    end

end
