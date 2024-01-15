module MasslessChargedParticlePlots

    using GeometricEquations
    using GeometricSolutions
    using GeometricProblems.Diagnostics

    using LaTeXStrings
    using Measures: mm
    using RecipesBase

    """
    Plots the solution of a massless charged particle together with the energy error.

    Arguments:
    * `sol <: GeometricSolution`
    * `equ <: GeometricProblem`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `xlims=:auto`: xlims for solution plot
    * `ylims=:auto`: ylims for solution plot
    """
    @userplot Plot_Massless_Charged_Particle
    @recipe function f(p::Plot_Massless_Charged_Particle; nplot=1, nt=:auto, xlims=:auto, ylims=:auto, elims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: GeometricSolution) || !(typeof(p.args[2]) <: GeometricProblem)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
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

        H, ΔH = compute_invariant_error(sol.t, sol.q, (t,q) -> invariants(equ)[:h](t,q,params))

        size   := (800,300)
        layout := (1,2)
        legend := :none

        right_margin := 10mm
        bottom_margin := 10mm

        guidefontsize := 14
        tickfontsize  := 10

        @series begin
            subplot := 1

            if ntime(sol) ≤ 200
                markersize := 5
            else
                markersize  := .5
                markercolor := 1
                linecolor   := 1
                markerstrokewidth := .5
                markerstrokecolor := 1
            end
            if nplot > 1
                seriestype := :scatter
            end

            if latex
                xguide := L"x_1"
                yguide := L"x_2"
            else
            # if backend() == Plots.PlotlyJSBackend()
                xguide := "x₁"
                yguide := "x₂"
            end
            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            sol.q[0:nplot:end, 1], sol.q[0:nplot:end, 2]
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
            ylims  := elims
            yformatter := :scientific
            sol.t[0:nplot:nt], ΔH[0:nplot:nt]
        end
    end


    """
    Plots the solution of a massless charged particle.

    Arguments:
    * `sol <: GeometricSolution`
    * `equ <: GeometricProblem`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `xlims=:auto`: xlims for solution plot
    * `ylims=:auto`: ylims for solution plot
    """
    @userplot Plot_Massless_Charged_Particle_Solution
    @recipe function f(p::Plot_Massless_Charged_Particle_Solution; nplot=1, nt=:auto, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: GeometricSolution) || !(typeof(p.args[2]) <: GeometricProblem)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
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
            markersize  := .5
            markercolor := 1
            linecolor   := 1
            markerstrokewidth := .5
            markerstrokecolor := 1
        end
        if nplot > 1
            seriestype := :scatter
        end

        legend := :none
        size := (400,400)

        # if backend() == Plots.PlotlyJSBackend() || backend() == Plots.GRBackend()
            guidefontsize := 18
            tickfontsize  := 12
        # else
        #     guidefontsize := 14
        #     tickfontsize  := 10
        # end

        # solution
        @series begin
            if latex
                xguide := L"x_1"
                yguide := L"x_2"
            else
            # if backend() == Plots.PlotlyJSBackend()
                xguide := "x₁"
                yguide := "x₂"
            end

            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            sol.q[0:nplot:nt, 1], sol.q[0:nplot:nt, 2]
        end
    end


    """
    Plots time traces of the energy error of a massless charged particle.

    Arguments:
    * `sol <: GeometricSolution`
    * `equ <: GeometricProblem`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    """
    @userplot Plot_Massless_Charged_Particle_Energy_Error
    @recipe function f(p::Plot_Massless_Charged_Particle_Energy_Error; nplot=1, nt=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: GeometricSolution) || !(typeof(p.args[2]) <: GeometricProblem)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
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

        H, ΔH = compute_invariant_error(sol.t, sol.q, (t,q) -> invariants(equ)[:h](t,q,params))

        size   := (800,200)
        legend := :none
        right_margin := 10mm

        # if backend() == Plots.GRBackend()
            left_margin := 10mm
            bottom_margin := 10mm
        # end

        guidefontsize := 12
        tickfontsize  := 10

        @series begin
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            # if backend() == Plots.PlotlyJSBackend()
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end

            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            sol.t[0:nplot:nt], ΔH[0:nplot:nt]
        end
    end


    """
    Plots time traces of the momentum error of a massless charged particle.

    Arguments:
    * `sol <: GeometricSolution`
    * `equ <: GeometricProblem`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `k=0`: index of momentum component (0 plots all components)
    """
    @userplot Plot_Massless_Charged_Particle_Momentum_Error
    @recipe function f(p::Plot_Massless_Charged_Particle_Momentum_Error; k=0, nplot=1, nt=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: GeometricSolution) || !(typeof(p.args[2]) <: GeometricProblem)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
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

        Δp = compute_momentum_error(sol.t, sol.q, sol.p, params)

        ntrace = (k == 0 ?    Δp.nd  :    1 )
        trange = (k == 0 ? (1:Δp.nd) : (k:k))

        size   := (800,200*ntrace)
        legend := :none
        layout := (ntrace,1)

        right_margin := 10mm

        if ntrace == 1 #&& backend() == Plots.GRBackend()
            left_margin := 10mm
        end

        guidefontsize := 12
        tickfontsize  := 10

        for i in trange
            @series begin
                if k == 0
                    subplot := i
                end
                if i == Δp.nd || k ≠ 0
                    if latex
                        xguide := L"t"
                    # if backend() == Plots.PlotlyJSBackend()
                    else
                        xguide := "t"
                    end
                else
                    xaxis := false
                    bottom_margin := 10mm
                end
                # if i < Δp.nd && k == 0
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
                    yguide := latexstring("p_$i (t) - \\vartheta_$i (t)")
                # if backend() == Plots.PlotlyJSBackend()
                else
                    yguide := "p" * subscript(i) * "(t) - ϑ" * subscript(i) * "(t)"
                end
                xlims  := (sol.t[0], Inf)
                yformatter := :scientific
                sol.t[0:nplot:nt], Δp[0:nplot:nt, i]
            end
        end
    end


    """
    Plots time traces of the solution of a massless charged particle trajectory and its energy error.

    Arguments:
    * `sol <: GeometricSolution`
    * `equ <: GeometricProblem`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    """
    @userplot Plot_Massless_Charged_Particle_Traces
    @recipe function f(p::Plot_Massless_Charged_Particle_Traces; nplot=1, nt=:auto, xlims=:auto, ylims=:auto, elims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: GeometricSolution) || !(typeof(p.args[2]) <: GeometricProblem)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
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

        H, ΔH = compute_invariant_error(sol.t, sol.q, (t,q) -> invariants(equ)[:h](t,q,params))

        size   := (800,600)
        legend := :none
        right_margin := 10mm

        # traces
        layout := (3,1)

        # if backend() == Plots.PlotlyJSBackend() || backend() == Plots.GRBackend()
            guidefontsize := 14
            tickfontsize  := 10
        # else
        #     guidefontsize := 12
        #     tickfontsize  := 10
        # end

        if latex
            ylabels = (L"x_1", L"x_2")
        # if backend() == Plots.PlotlyJSBackend()
        else
            ylabels = ("x₁", "x₂")
        end

        lims = (xlims, ylims, elims)

        for i in 1:2
            @series begin
                subplot := i
                yguide := ylabels[i]
                xlims  := (sol.t[0], Inf)
                ylims  := lims[i]
                xaxis  := false

                # if i == 1 && backend() == Plots.PGFPlotsXBackend()
                #     bottom_margin := -4mm
                # elseif i == 1 && backend() == Plots.GRBackend()
                #     bottom_margin := -2mm
                # end
                # if i == 2 && backend() == Plots.PGFPlotsXBackend()
                #     top_margin := -2mm
                #     bottom_margin := -2mm
                # elseif i == 2 && backend() == Plots.GRBackend()
                #     top_margin := -1mm
                #     bottom_margin := -1mm
                # end

                sol.t[0:nplot:nt], sol.q[0:nplot:nt, i]
            end
        end

        @series begin
            subplot := 3
            # if backend() == Plots.PGFPlotsXBackend()
            #     top_margin := -4mm
            # elseif backend() == Plots.GRBackend()
            #     top_margin := -2mm
            # end
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            # if backend() == Plots.PlotlyJSBackend()
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end
            xlims := (sol.t[0], Inf)
            ylims := elims
            yformatter := :scientific
            sol.t[0:nplot:nt], ΔH[0:nplot:nt]
        end
    end

end
