module LotkaVolterra2dPlots

    using ..Diagnostics
    using GeometricIntegrators.Equations
    using GeometricIntegrators.Solutions

    using LaTeXStrings
    using Plots
    using Plots.PlotMeasures
    # using RecipesBase


    """
    Plots the solution of a 2D Lotka-Volterra model together with the energy error.

    Arguments:
    * `sol <: Solution`
    * `equ <: Equation`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `xlims=:auto`: xlims for solution plot
    * `ylims=:auto`: ylims for solution plot
    * `latex=true`: use LaTeX guides
    """
    @userplot Plot_Lotka_Volterra_2d
    @recipe function f(p::Plot_Lotka_Volterra_2d; nplot=1, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: Equation)
            error("Lotka-Volterra plots should be given two arguments: a solution and an equation. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        equ = p.args[2]
        params = equ.parameters

        H, ΔH = compute_invariant_error(sol.t, sol.q, (t,q) -> equ.h(t,q,params))

        size   := (1000,400)
        layout := @layout [solPlot{0.3w,1.0h} EPlot]
        legend := :none

        right_margin := 10mm
        bottom_margin := 10mm

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
                xguide := L"x_1"
                yguide := L"x_2"
            else
                xguide := "x₁"
                yguide := "x₂"
            end
            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end

        @series begin
            subplot := 2
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            ylims  := :auto
            yformatter := :scientific
            right_margin := 10mm
            sol.t[0:nplot:end], ΔH[0:nplot:end]
        end
    end


    """
    Plots the solution of a 2D Lotka-Volterra model.

    Arguments:
    * `sol <: Solution`
    * `equ <: Equation`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `xlims=:auto`: xlims for solution plot
    * `ylims=:auto`: ylims for solution plot
    * `latex=true`: use LaTeX guides
    """
    @userplot Plot_Lotka_Volterra_2d_Solution
    @recipe function f(p::Plot_Lotka_Volterra_2d_Solution; nplot=1, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: Equation)
            error("Lotka-Volterra plots should be given two arguments: a solution and an equation. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        equ = p.args[2]
        params = equ.parameters

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
                xguide := L"x_1"
                yguide := L"x_2"
            else
                xguide := "x₁"
                yguide := "x₂"
            end
            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            guidefont := font(18)
            tickfont := font(12)
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end
    end


    """
    Plots time traces of the solution of a 2D Lotka-Volterra model and its energy error.

    Arguments:
    * `sol <: Solution`
    * `equ <: Equation`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `latex=true`: use LaTeX guides
    """
    @userplot Plot_Lotka_Volterra_2d_Traces
    @recipe function f(p::Plot_Lotka_Volterra_2d_Traces; nplot=1, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: Equation)
            error("Lotka-Volterra plots should be given two arguments: a solution and an equation. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        equ = p.args[2]
        params = equ.parameters

        H, ΔH = compute_invariant_error(sol.t, sol.q, (t,q) -> equ.h(t,q,params))

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
                yguide := ylabels[i]
                xlims  := (sol.t[0], Inf)
                xaxis := false
                right_margin := 10mm
                sol.t[0:nplot:end], sol.q[i,0:nplot:end]
            end
        end

        @series begin
            subplot := 3
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
            sol.t[0:nplot:end], ΔH[0:nplot:end]
        end
    end

end
