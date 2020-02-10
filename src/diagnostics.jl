module Diagnostics

    using RecipesBase
    using GeometricIntegrators.Solutions

    export compute_energy_error, compute_momentum_error, compute_error_drift


    function compute_one_form(t::TimeSeries, q::DataSeries, one_form::Function=ϑ)
        p = similar(q)

        for i in axes(p,2)
            for k in axes(p,1)
                p[k,i] = one_form(t[i], q[:,i], k)
            end
        end

        return e
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

    function compute_error_drift(e::DataSeries{T,1}, lint=100) where {T}
        nint = div(e.nt, lint)

        # d = SDataSeries(T, nint)
        #
        # for i in 1:nint
        #     i1 = lint*(i-1)+1
        #     i2 = lint*i
        #     d[i] = maximum(e[i1:i2])
        # end
        #
        # d .-= d[0]

        drift = zeros(T, nint)

        for i in 1:nint
            i1 = lint*(i-1)+1
            i2 = lint*i
            drift[i] = maximum(e[i1:i2])
        end

        drift .- drift[1]
    end


    @userplot PlotEnergyError
    @recipe function f(p::PlotEnergyError)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Energy error plot should be given a solution. Got: $(typeof(p.args))")
        end

        sol = p.args[1]

        legend := :none
        size := (800,400)

        H, ΔH = compute_energy_error(sol.t, sol.q)

        @series begin
            xlabel := "t"
            ylabel := "[H(t) - H(0)] / H(0)"
            sol.t, ΔH
        end
    end


end
