module PlotRecipes

using GeometricSolutions
using LaTeXStrings
using Measures: mm
using RecipesBase


subscript(i::Integer) = i < 0 ? error("$i is negative") : join('₀' + d for d in reverse(digits(i)))


@userplot PlotEnergyError
@recipe function f(p::PlotEnergyError; energy=nothing, nplot=1, nt=:auto, latex=true)
    if length(p.args) == 1 && typeof(p.args[1]) <: GeometricSolution
        sol = p.args[1]
        t = sol.t
        if typeof(sol) <: Union{SolutionPODE,SolutionPDAE}
            H, ΔH = compute_invariant_error(sol.t, sol.q, sol.p, energy)
        else
            H, ΔH = compute_invariant_error(sol.t, sol.q, energy)
        end
    elseif length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
        t = p.args[1]
        ΔH = p.args[2]
        @assert length(t) == length(ΔH)
    else
        error("Energy error plot should be given a solution or a timeseries and a data series. Got: $(typeof(p.args))")
    end

    if nt == :auto
        nt = ntime(t)
    end

    if nt > ntime(t)
        nt = ntime(t)
    end

    legend := :none
    size := (800, 400)

    @series begin
        if latex
            xguide := L"t"
            yguide := L"[H(t) - H(0)] / H(0)"
        else
            xguide := "t"
            yguide := "[H(t) - H(0)] / H(0)"
        end
        xlims := (t[0], Inf)
        yformatter := :scientific
        guidefontsize := 18
        tickfontsize := 12
        right_margin := 12mm
        t[0:nplot:nt], ΔH[0:nplot:nt]
    end
end


@userplot PlotEnergyDrift
@recipe function f(p::PlotEnergyDrift; nt=:auto, latex=true)
    if length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
        t = p.args[1]
        d = p.args[2]
        @assert length(t) == length(d)
    else
        error("Energy drift plot should be given a timeseries and a dataseries. Got: $(typeof(p.args))")
    end

    if nt == :auto
        nt = ntime(t)
    end

    if nt > ntime(t)
        nt = ntime(t)
    end

    legend := :none
    size := (800, 400)

    @series begin
        seriestype := :scatter
        if latex
            xguide := L"t"
            yguide := L"\Delta H"
        else
            xguide := "t"
            yguide := "ΔH"
        end
        xlims := (t[0], Inf)
        yformatter := :scientific
        guidefontsize := 18
        tickfontsize := 12
        right_margin := 12mm
        t[1:nt], d[1:nt]
    end
end


@userplot PlotConstraintError
@recipe function f(p::PlotConstraintError; nplot=1, nt=:auto, k=0, latex=true, plot_title=nothing)
    if length(p.args) == 1 && typeof(p.args[1]) <: GeometricSolution
        sol = p.args[1]
        t = sol.t
        Δp = compute_momentum_error(sol)
        # TODO: fix
    elseif length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
        t = p.args[1]
        Δp = p.args[2]
        @assert ntime(t) == ntime(Δp)
    else
        error("Constraint error plots should be given a solution or a timeseries and a dataseries. Got: $(typeof(p.args))")
    end

    if nt == :auto
        nt = ntime(t)
    end

    if nt > ntime(t)
        nt = ntime(t)
    end

    nd = length(Δp[begin])

    ntrace = (k == 0 ? nd : 1)
    trange = (k == 0 ? (1:nd) : (k:k))

    size := (800, 200 * ntrace)
    legend := :none
    layout := (ntrace, 1)

    for i in trange
        @series begin
            if k == 0
                subplot := i
            end
            if plot_title !== nothing
                if i == trange[begin]
                    title := plot_title
                else
                    title := "   "
                end
            end
            if i == nd || k ≠ 0
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
            xlims := (t[0], Inf)
            yformatter := :scientific
            guidefontsize := 18
            tickfontsize := 12
            right_margin := 24mm
            right_margin := 12mm
            t[0:nplot:nt], Δp[i, 0:nplot:nt]
        end
    end
end


@userplot PlotLagrangeMultiplier
@recipe function f(p::PlotLagrangeMultiplier; nplot=1, nt=:auto, k=0, latex=true, plot_title=nothing)
    if length(p.args) == 1 && typeof(p.args[1]) <: GeometricSolution
        sol = p.args[1]
        t = sol.t
        λ = sol.λ
    elseif length(p.args) == 2 && typeof(p.args[1]) <: TimeSeries && typeof(p.args[2]) <: DataSeries
        t = p.args[1]
        λ = p.args[2]
        @assert ntime(t) == ntime(λ)
    else
        error("Lagrange multiplier plots should be given a solution or a timeseries and a data series. Got: $(typeof(p.args))")
    end

    if nt == :auto
        nt = ntime(t)
    end

    if nt > ntime(t)
        nt = ntime(t)
    end

    nd = length(λ[begin])

    ntrace = (k == 0 ? nd : 1)
    trange = (k == 0 ? (1:nd) : (k:k))

    size := (800, 200 * ntrace)
    legend := :none
    layout := (ntrace, 1)

    right_margin := 10mm

    if ntrace == 1 && backend() == Plots.GRBackend()
        left_margin := 10mm
    end

    guidefontsize := 12
    tickfontsize := 10

    for i in trange
        @series begin
            if k == 0
                subplot := i
            end
            if plot_title !== nothing
                if i == trange[begin]
                    title := plot_title
                else
                    title := "   "
                end
            end
            if i == nd || k ≠ 0
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
            xlims := (t[0], Inf)
            t[0:nplot:nt], λ[i, 0:nplot:nt]
        end
    end
end

end
