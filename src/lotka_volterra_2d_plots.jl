module LotkaVolterra2dPlots

    using GeometricIntegrators.Common
    using GeometricIntegrators.Equations
    using GeometricIntegrators.Solutions
    using GeometricProblems.Diagnostics

    using LaTeXStrings
    using Measures: mm
    using RecipesBase


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
    @recipe function f(p::Plot_Lotka_Volterra_2d; nplot=1, nt=:auto, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: Equation)
            error("Lotka-Volterra plot should be given two arguments: a solution and an equation. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        equ = p.args[2]
        params = equ.parameters

        if nt == :auto
            nt = ntime(sol)
        end

        if nt > ntime(sol)
            nt = ntime(sol)
        end

        H, ΔH = compute_invariant_error(sol.t, sol.q, (t,q) -> equ.invariants[:h](t,q,params))

        size   := (1000,400)
        layout := (1,2)
        legend := :none

        right_margin := 10mm
        bottom_margin := 10mm

        guidefontsize := 18
        tickfontsize  := 12

        @series begin
            subplot := 1
            left_margin := 10mm

            if ntime(sol) ≤ 200
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
            sol.q[1,0:nplot:nt], sol.q[2,0:nplot:nt]
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
            sol.t[0:nplot:nt], ΔH[0:nplot:nt]
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
    @recipe function f(p::Plot_Lotka_Volterra_2d_Solution; nplot=1, nt=:auto, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: Equation)
            error("Lotka-Volterra solution plot should be given two arguments: a solution and an equation. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        equ = p.args[2]
        params = equ.parameters

        if nt == :auto
            nt = ntime(sol)
        end

        if nt > ntime(sol)
            nt = ntime(sol)
        end

        if ntime(sol) ≤ 200
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
            guidefontsize := 18
            tickfontsize  := 12
            sol.q[1,0:nplot:nt], sol.q[2,0:nplot:nt]
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
    @recipe function f(p::Plot_Lotka_Volterra_2d_Traces; nplot=1, nt=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: Equation)
            error("Lotka-Volterra traces plot should be given two arguments: a solution and an equation. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        equ = p.args[2]
        params = equ.parameters

        if nt == :auto
            nt = ntime(sol)
        end

        if nt > ntime(sol)
            nt = ntime(sol)
        end

        H, ΔH = compute_invariant_error(sol.t, sol.q, (t,q) -> equ.invariants[:h](t,q,params))

        size   := (800,600)
        legend := :none
        guidefontsize := 18
        tickfontsize  := 12

        # traces
        
        layout := (3,1)

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
                sol.t[0:nplot:nt], sol.q[i,0:nplot:nt]
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
            guidefontsize := 18
            tickfontsize  := 12
            right_margin := 10mm
            sol.t[0:nplot:nt], ΔH[0:nplot:nt]
        end
    end

end
