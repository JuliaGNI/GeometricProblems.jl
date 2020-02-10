module Diagnostics

    using RecipesBase
    using GeometricIntegrators.Solutions

    export compute_energy_error, compute_momentum_error


    function compute_energy_error(t::TimeSeries, q::DataSeries{T}, energy::Function=hamiltonian) where {T}
        h = SDataSeries(T, q.nt)
        e = SDataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i])
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_momentum_error(t::TimeSeries, q::DataSeries{T}, p::DataSeries{T}, one_form::Function=ϑ) where {T}
        e = SDataSeries(T, p.nd, p.nt)

        for i in axes(p,2)
            for k in axes(p,1)
                e[k,i] = p[k,i] - ϑ(t[i], q[:,i], k)
            end
        end

        return e
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
