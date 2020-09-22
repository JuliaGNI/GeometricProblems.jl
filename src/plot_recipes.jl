module PlotRecipes

    using Plots
    using Plots.PlotMeasures
    using RecipesBase
    using LaTeXStrings

    using GeometricIntegrators.Solutions
    using ..Diagnostics


    subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))


    @userplot PlotEnergyError
    @recipe function f(p::PlotEnergyError; energy=nothing, nplot=1, latex=true)
        if length(p.args) == 1 && typeof(p.args[1]) <: Solution
            sol = p.args[1]
            t = sol.t
            if typeof(sol) <: Union{SolutionPODE, SolutionPDAE}
                H, ΔH = compute_invariant_error(sol.t, sol.q, sol.p, energy)
            else
                H, ΔH = compute_invariant_error(sol.t, sol.q, energy)
            end
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
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
            t[0:nplot:end], ΔH[0:nplot:end]
        end
    end


    @userplot PlotEnergyDrift
    @recipe function f(p::PlotEnergyDrift; latex=true)
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
            seriestype := :scatter
            if latex
                xguide := L"t"
                yguide := L"\Delta H"
            else
                xguide := "t"
                yguide := "ΔH"
            end
            xlims  := (t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
            t[1:end], d[1:end]
        end
    end


    @userplot PlotMomentumError
    @recipe function f(p::PlotMomentumError; nplot=1, k=0, latex=true)
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

        ntrace = (k == 0 ?    Δp.nd  :    1 )
        trange = (k == 0 ? (1:Δp.nd) : (k:k))

        size   := (800,200*ntrace)
        legend := :none
        layout := (ntrace,1)

        for i in trange
            @series begin
                if k == 0
                    subplot := i
                end
                if i == Δp.nd || k ≠ 0
                    if latex
                        xguide := L"t"
                    else
                        xguide := "t"
                    end
                else
                    xaxis := false
                end
                if latex
                    yguide := latexstring("p_$i (t) - \\vartheta_$i (t)")
                else
                    yguide := "p" * subscript(i) * "(t) - ϑ" * subscript(i) * "(t)"
                end
                xlims  := (t[0], Inf)
                yformatter := :scientific
                guidefont := font(18)
                tickfont := font(12)
                right_margin := 24mm
                right_margin := 12mm
                t[0:nplot:end], Δp[i,0:nplot:end]
            end
        end
    end


    @userplot PlotLagrangeMultiplier
    @recipe function f(p::PlotLagrangeMultiplier; k=0, nplot=1, latex=true)
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

        ntrace = (k == 0 ?    λ.nd  :    1 )
        trange = (k == 0 ? (1:λ.nd) : (k:k))

        size   := (800,200*ntrace)
        legend := :none
        layout := (ntrace,1)

        right_margin  := 10mm

        if ntrace == 1 && backend() == Plots.GRBackend()
            left_margin := 10mm
        end

        guidefont := font(12)
        tickfont  := font(10)

        for i in trange
            @series begin
                if k == 0
                    subplot := i
                end
                if i == λ.nd || k ≠ 0
                    if latex
                        xguide := L"t"
                    else
                        xguide := "t"
                    end
                    bottom_margin := 10mm
                else
                    xaxis := false
                end
                # if i < λ.nd && k == 0
                #     if backend() == Plots.PGFPlotsXBackend()
                #         bottom_margin := -6mm
                #     elseif backend() == Plots.GRBackend()
                #         bottom_margin := -3mm
                #     end
                # end
                # if i > 1 && k == 0
                #     if backend() == Plots.PGFPlotsXBackend()
                #         top_margin := -6mm
                #     elseif backend() == Plots.GRBackend()
                #         top_margin := -3mm
                #     end
                # end
                if latex
                    yguide := latexstring("\\lambda_$i (t)")
                else
                    yguide := "λ" * subscript(i) * "(t)"
                end
                xlims  := (t[0], Inf)
                t[0:nplot:end], λ[i,0:nplot:end]
            end
        end
    end


end
