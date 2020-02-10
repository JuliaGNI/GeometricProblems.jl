module LotkaVolterra2d

    using GeometricIntegrators.Solutions
    using Plots: @layout
    using RecipesBase
    using Reexport

    @reexport using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem

    import ..Diagnostics: compute_energy_error, compute_momentum_error
    import GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: hamiltonian, ϑ, ϑ₁, ϑ₂, ω

    export hamiltonian, ϑ, ϑ₁, ϑ₂, ω
    export compute_energy_error, compute_momentum_error
    export f


    compute_energy_error(t,q) = compute_energy_error(t,q,hamiltonian)
    compute_momentum_error(t,q) = compute_momentum_error(t,q,ϑ)


    function ϑ(t::Number, q::Vector, k::Int)
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
            xlabel := "x₁"
            ylabel := "x₂"
            aspectratio := 1
            subplot := 1
            sol.q[1,:], sol.q[2,:]
        end

        @series begin
            xlabel := "t"
            ylabel := "[H(t) - H(0)] / H(0)"
            subplot := 2
            sol.t, ΔH
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

        ylabels = ("x₁", "x₂")

        for i in 1:2
            @series begin
                # xlabel := "t"
                ylabel := ylabels[i]
                subplot := i
                sol.t, sol.q[i,:]
            end
        end

        @series begin
            xlabel := "t"
            ylabel := "[H(t) - H(0)] / H(0)"
            subplot := 3
            sol.t, ΔH
        end
    end

end
