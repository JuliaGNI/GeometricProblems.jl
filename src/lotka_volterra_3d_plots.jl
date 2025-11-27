module LotkaVolterra3dPlots

    using GeometricSolutions
    using LaTeXStrings
    using Measures: mm
    using RecipesBase

    # export compute_energy_error, compute_momentum_error


    compute_energy_error(t,q,params) = compute_invariant_error(t,q, (t,q) -> hamiltonian(t,q,params))
    # compute_momentum_error(t,q,p) = compute_momentum_error(t,q,p,ϑ)


    @userplot PlotLotkaVolterra3d
    @recipe function f(p::PlotLotkaVolterra3d; latex=true)
        # if length(p.args) != 2 || !(typeof(p.args[1][1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
        #     error("Lotka-Volterra plots should be given a solution. Got: $(typeof(p.args))")
        # end
        sol = p.args[1]
        params = p.args[2]

        H, ΔH = compute_energy_error(sol.t, sol.q, params);

        size   := (800,1200)
        legend := :none
        guidefontsize := 12
        tickfontsize  := 8
        left_margin := 12mm
        right_margin := 8mm

        # traces
        layout := (4,1)

        if latex
            ylabels = (L"x_1", L"x_2", L"x_3")
        else
            ylabels = ("x₁", "x₂", "x₃")
        end

        for i in 1:3
            @series begin
                subplot := i
                yguide := ylabels[i]
                xlims  := (sol.t[0], Inf)
                xaxis := false
                sol.t, sol.q[i,:]
            end
        end

        @series begin
            subplot := 4
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            sol.t, ΔH
        end
    end

end
