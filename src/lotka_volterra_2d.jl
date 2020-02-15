module LotkaVolterra2d

    using Plots
    using Plots.PlotMeasures
    using Plots: @layout
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
    @recipe function f(p::PlotLotkaVolterra2d)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Lotka-Volterra plots should be given a solution. Got: $(typeof(p.args))")
        end
        sol = p.args[1]

        legend := :none
        size := (800,400)

        H, ΔH = compute_energy_error(sol.t, sol.q);

        # solution
        layout := @layout [solPlot{0.4w,1.0h} EPlot]

        @series begin
            subplot := 1
            # seriestype := :scatter
            xlabel := L"x_1"
            ylabel := L"x_2"
            aspectratio := 1
            guidefont := font(18)
            tickfont := font(12)
            sol.q[1,:], sol.q[2,:]
        end

        @series begin
            subplot := 2
            xlabel := L"t"
            ylabel := L"[H(t) - H(0)] / H(0)"
            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
            sol.t, ΔH
        end
    end


    @userplot PlotLotkaVolterra2dSolution
    @recipe function f(p::PlotLotkaVolterra2dSolution)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Lotka-Volterra plots should be given a solution. Got: $(typeof(p.args))")
        end
        sol = p.args[1]

        legend := :none
        size := (400,400)

        # solution
        @series begin
            # seriestype := :scatter
            xlabel := L"x_1"
            ylabel := L"x_2"
            aspectratio := 1
            guidefont := font(18)
            tickfont := font(12)
            sol.q[1,:], sol.q[2,:]
        end
    end


    @userplot PlotLotkaVolterra2dTraces
    @recipe function f(p::PlotLotkaVolterra2dTraces)
        if length(p.args) != 1 || !(typeof(p.args[1]) <: Solution)
            error("Lotka-Volterra plots should be given a solution. Got: $(typeof(p.args))")
        end
        sol = p.args[1]

        legend := :none
        size := (800,600)

        H, ΔH = compute_energy_error(sol.t, sol.q);

        # traces
        layout := @layout [x₁Plot
                           x₂Plot
                           EPlot]

        ylabels = (L"x_1", L"x_2")

        for i in 1:2
            @series begin
                subplot := i
                ylabel := ylabels[i]
                xlims  := (sol.t[0], Inf)
                xaxis := false
                guidefont := font(18)
                tickfont := font(12)
                right_margin := 10mm
                sol.t, sol.q[i,:]
            end
        end

        @series begin
            subplot := 3
            xlabel := L"t"
            ylabel := L"[H(t) - H(0)] / H(0)"
            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            guidefont := font(18)
            tickfont := font(12)
            right_margin := 10mm
            sol.t, ΔH
        end
    end

end
