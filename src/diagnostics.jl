module Diagnostics

    using Plots
    using Plots.PlotMeasures
    # using Plots: convert_sci_unicode
    using Showoff
    using Showoff: format_fixed_scientific
    using RecipesBase

    using GeometricIntegrators.Solutions

    export compute_energy_error, compute_momentum_error, compute_error_drift


    subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))


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


    function energy_ticks(ΔH)
        ticks = optimize_ticks(extrema(ΔH)...; k_min=4)[1]
        # labels = Showoff.showoff(ticks, :scientific)
        labels = format_fixed_scientific.(ticks, 2, false)

        # if backend_name() == :gr
        #     # for GR
        #     yticks = (ticks, convert_sci_unicode.(labels))
        # else
        #     # for most backends
        #     yticks = (ticks, labels)
        # end
        # yticks

        (ticks, labels)
    end


    @userplot PlotEnergyError
    @recipe function f(p::PlotEnergyError)
        if length(p.args) == 1 && typeof(p.args[1]) <: Solution
            sol = p.args[1]
            t = sol.t
            H, ΔH = compute_energy_error(sol.t, sol.q)
        elseif length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
            t  = p.args[1]
            ΔH = p.args[2]
            @assert length(t) == length(ΔH)
        else
            error("Energy error plot should be given a solution or a timeseries and a data series. Got: $(typeof(p.args))")
        end

        legend := :none
        size := (800,400)

        @series begin
            xlabel := "t"
            ylabel := "[H(t) - H(0)] / H(0)"
            yticks := energy_ticks(ΔH)
            right_margin := 10mm
            xlims  := (t[0], Inf)
            yformatter := :scientific
            # widen := false
            t, ΔH
        end
    end


    @userplot PlotEnergyDrift
    @recipe function f(p::PlotEnergyDrift)
        if length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
            t = p.args[1]
            d = p.args[2]
            @assert length(t) == length(d)
        else
            error("Energy drift plot should be given a timeseries and a data series. Got: $(typeof(p.args))")
        end

        legend := :none
        size := (800,400)

        @series begin
            xlabel := "t"
            ylabel := "ΔH Drift"
            # yticks := energy_ticks(ΔH)
            right_margin := 10mm
            xlims  := (t[0], Inf)
            yformatter := :scientific
            # widen := false
            t, d
        end
    end


    @userplot PlotMomentumError
    @recipe function f(p::PlotMomentumError)
        if length(p.args) == 1 && typeof(p.args[1]) <: Solution
            sol = p.args[1]
            t = sol.t
            Δp = compute_momentum_error(sol.t, sol.q, sol.p)
        elseif length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
            t  = p.args[1]
            Δp = p.args[2]
            @assert t.n == Δp.nt
        else
            error("Momentum error plots should be given a solution or a timeseries and a data series. Got: $(typeof(p.args))")
        end

        legend := :none
        size := (800,200*Δp.nd)

        # traces
        layout := (Δp.nd,1)

        for i in 1:Δp.nd
            @series begin
                subplot := i
                if i == Δp.nd
                    xlabel := "t"
                end
                ylabel := "p" * subscript(i) * "(t) - ϑ" * subscript(i) * "(t)"
                xlims  := (t[0], Inf)
                right_margin := 10mm
                t, Δp[i,:]
            end
        end
    end


    @userplot PlotLagrangeMultiplier
    @recipe function f(p::PlotLagrangeMultiplier)
        if length(p.args) == 1 && typeof(p.args[1]) <: Solution
            sol = p.args[1]
            t = sol.t
            λ = sol.λ
        elseif length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
            t = p.args[1]
            λ = p.args[2]
            @assert t.n == λ.nt
        else
            error("Lagrange multiplier plots should be given a solution or a timeseries and a data series. Got: $(typeof(p.args))")
        end

        legend := :none
        size := (800,200*λ.nd)

        # traces
        layout := (λ.nd,1)

        for i in 1:λ.nd
            @series begin
                subplot := i
                if i == λ.nd
                    xlabel := "t"
                end
                ylabel := "λ" * subscript(i) * "(t)"
                xlims  := (t[0], Inf)
                right_margin := 10mm
                t, λ[i,:]
            end
        end
    end


end
