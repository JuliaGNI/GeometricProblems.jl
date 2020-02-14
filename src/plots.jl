module Plots

    using Plots
    using Plots.PlotMeasures
    using LaTeXStrings
    using RecipesBase

    using GeometricIntegrators.Solutions: DataSeries, TimeSeries, Solution
    using ..Diagnostics


    @userplot PlotEnergyError
    @recipe function f(p::PlotEnergyError; energy=nothing)
        if length(p.args) == 1 && typeof(p.args[1]) <: Solution
            sol = p.args[1]
            t = sol.t
            H, ΔH = compute_energy_error(sol.t, sol.q, energy)
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
            xlabel := L"t"
            ylabel := L"[H(t) - H(0)] / H(0)"
            xlims  := (t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
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
            xlabel := L"t"
            ylabel := L"\Delta H"
            xlims  := (t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
            t[1:end], d[1:end]
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
                    xlabel := L"t"
                else
                    xaxis := false
                end
                ylabel := latexstring("p_$i (t) - \\vartheta_$i (t)")
                xlims  := (t[0], Inf)
                yformatter := :scientific
                guidefont := font(18)
                tickfont := font(12)
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
                    xlabel := L"t"
                else
                    xaxis := false
                end
                ylabel := latexstring("\\lambda_$i (t)")
                xlims  := (t[0], Inf)
                guidefont := font(18)
                tickfont := font(12)
                right_margin := 10mm
                t, λ[i,:]
            end
        end
    end


end
