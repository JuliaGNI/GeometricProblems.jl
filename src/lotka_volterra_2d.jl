module LotkaVolterra2d

    using Plots
    using Plots.PlotMeasures
    using RecipesBase
    using LaTeXStrings
    using Reexport

    using GeometricIntegrators.Solutions

    @reexport using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem

    import ..Diagnostics: compute_invariant_error, compute_momentum_error
    import GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: hamiltonian, ϑ, ϑ₁, ϑ₂, ω

    export hamiltonian, ϑ, ϑ₁, ϑ₂, ω
    export compute_energy_error, compute_momentum_error
    export f


    compute_energy_error(t,q) = compute_invariant_error(t,q,hamiltonian)
    compute_momentum_error(t,q,p) = compute_momentum_error(t,q,p,ϑ)


    function ϑ(t::Number, q::AbstractVector, k::Int)
        if k == 1
            ϑ₁(t, q)
        elseif k == 2
            ϑ₂(t, q)
        else
            nothing
        end
    end


    @userplot PlotLotkaVolterra2d
    @recipe function f(p::PlotLotkaVolterra2d; nplot=1, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Lotka-Volterra plots should be given a solution. Got: $(typeof(p.args))")
        end
        sol = p.args[1]

        H, ΔH = compute_energy_error(sol.t, sol.q);

        size   := (1000,400)
        layout := @layout [solPlot{0.3w,1.0h} EPlot]
        legend := :none

        guidefont := font(18)
        tickfont  := font(12)

        @series begin
            subplot := 1

            if sol.nt ≤ 200
                markersize := 5
            else
                markersize  := 1
                markercolor := 1
                linecolor   := 1
                markerstrokewidth := 1
                markerstrokecolor := 1
            end

            # seriestype := :scatter
            if latex
                xlabel := L"x_1"
                ylabel := L"x_2"
            else
                xlabel := "x₁"
                ylabel := "x₂"
            end
            xlims  := xlims
            ylims  := ylims
            aspectratio := 1
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end

        @series begin
            subplot := 2
            if latex
                xlabel := L"t"
                ylabel := L"[H(t) - H(0)] / H(0)"
            else
                xlabel := "t"
                ylabel := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            right_margin := 10mm
            sol.t[0:nplot:end], ΔH[0:nplot:end]
        end
    end


    @userplot PlotLotkaVolterra2dSolution
    @recipe function f(p::PlotLotkaVolterra2dSolution; nplot=1, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Lotka-Volterra plots should be given a solution. Got: $(typeof(p.args))")
        end
        sol = p.args[1]

        if sol.nt ≤ 200
            markersize := 5
        else
            markersize  := 1
            markercolor := 1
            linecolor   := 1
            markerstrokewidth := 1
            markerstrokecolor := 1
        end

        legend := :none
        size := (400,400)

        # solution
        @series begin
            # seriestype := :scatter
            if latex
                xlabel := L"x_1"
                ylabel := L"x_2"
            else
                xlabel := "x₁"
                ylabel := "x₂"
            end
            xlims  := xlims
            ylims  := ylims
            aspectratio := 1
            guidefont := font(18)
            tickfont := font(12)
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end
    end


    @userplot PlotLotkaVolterra2dTraces
    @recipe function f(p::PlotLotkaVolterra2dTraces; latex=true)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Lotka-Volterra plots should be given a solution. Got: $(typeof(p.args))")
        end
        sol = p.args[1]

        H, ΔH = compute_energy_error(sol.t, sol.q);

        size   := (800,600)
        legend := :none
        guidefont := font(18)
        tickfont  := font(12)

        # traces
        layout := @layout [x₁Plot
                           x₂Plot
                           EPlot]

        if latex
            ylabels = (L"x_1", L"x_2")
        else
            ylabels = ("x₁", "x₂")
        end

        for i in 1:2
            @series begin
                subplot := i
                ylabel := ylabels[i]
                xlims  := (sol.t[0], Inf)
                xaxis := false
                right_margin := 10mm
                sol.t, sol.q[i,:]
            end
        end

        @series begin
            subplot := 3
            if latex
                xlabel := L"t"
                ylabel := L"[H(t) - H(0)] / H(0)"
            else
                xlabel := "t"
                ylabel := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
            sol.t, ΔH
        end
    end

end
