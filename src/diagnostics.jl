module Diagnostics

    using GeometricIntegrators.Solutions

    export compute_one_form, compute_invariant, compute_energy_error, compute_momentum_error, compute_error_drift


    subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))


    function compute_one_form(t::TimeSeries, q::DataSeries, one_form::Function=ϑ)
        p = similar(q)

        for i in axes(p,2)
            for k in axes(p,1)
                p[k,i] = one_form(t[i], q[:,i], k)
            end
        end

        return p
    end

    function compute_invariant(t::TimeSeries, q::DataSeries{T}, invariant::Function) where {T}
        invds = SDataSeries(T, 2, q.nt)
        for i in axes(q,2)
            invds[1,i] = invariant(t[i], q[:,i])
            invds[2,i] = (invds[1,i] - invds[1,0]) / invds[1,0]
        end
        return invds
    end

    function compute_energy_error(t::TimeSeries, q::DataSeries{T}, energy::Function) where {T}
        h = SDataSeries(T, q.nt)
        e = SDataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = energy(t[i], q[:,i])
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_momentum_error(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, one_form::Function) where {T}
        e = SDataSeries(T, p.nd, p.nt)

        for i in axes(p,2)
            for k in axes(p,1)
                e[k,i] = p[k,i] - one_form(t[i], q[:,i], k)
            end
        end

        return e
    end

    function compute_momentum_error(p::DataSeries{DT}, ϑ::DataSeries{DT}) where {DT}
        @assert p.nd == ϑ.nd
        @assert p.nt == ϑ.nt
        @assert p.ni == ϑ.ni

        err = SDataSeries(DT, p.nd, p.nt, p.ni)

        for k in 1:p.ni
            for j in 0:p.nt
                for i in 1:p.nd
                    err[i,j,k] = p[i,j,k] - ϑ[i,j,k]
                end
            end
        end

        return err
    end

    function compute_error_drift(t::TimeSeries, e::DataSeries{T,1}, lint=100) where {T}
        @assert mod(t.n, lint) == 0

        nint   = div(e.nt, lint)
        Tdrift = TimeSeries(nint, t.Δt*lint)
        Edrift = SDataSeries(T, nint)

        compute_timeseries!(Tdrift, t[0])

        for i in 1:nint
            i1 = lint*(i-1)+1
            i2 = lint*i
            Edrift[i] = maximum(e[i1:i2])
        end

        (Tdrift, Edrift)
    end

end
