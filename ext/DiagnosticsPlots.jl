module DiagnosticsPlots

using Makie
using GeometricEquations: invariants, parameters
# Import types/functions explicitly: `using Makie` also exports a `TimeSeries`
# (its recipe), which would otherwise clash with `GeometricSolutions.TimeSeries`.
using GeometricSolutions: TimeSeries, ScalarDataSeries, DataSeries, GeometricSolution,
                          SolutionPODE, SolutionPDAE, ntime,
                          compute_invariant_error, compute_momentum_error,
                          compute_error_drift

import GeometricProblems.Diagnostics
import GeometricProblems.Diagnostics: subscript

# Downsampled index range 0:nplot:nt (nt = :auto → all stored steps).
function _steprange(t, nplot, nt)
    n = nt === :auto ? ntime(t) : min(nt, ntime(t))
    return 0:nplot:n
end

# Compute the relative energy error of a solution from its `:h` invariant.
function _energy_error(sol; energy = nothing)
    h = energy === nothing ? invariants(sol.problem)[:h] : energy
    params = parameters(sol.problem)
    if sol isa Union{SolutionPODE, SolutionPDAE}
        _, ΔH = compute_invariant_error(sol.t, sol.q, sol.p, params, h)
    else
        _, ΔH = compute_invariant_error(sol.t, sol.q, params, h)
    end
    return ΔH
end

# Stack one axis per component of a per-DOF time series (constraint error, λ, …).
function _plot_components(t, d, label; nplot = 1, nt = :auto, k = 0, latex = true, plot_title = nothing)
    r  = _steprange(t, nplot, nt)
    ts = [t[j] for j in r]
    nd = length(d[begin])
    trange = k == 0 ? (1:nd) : (k:k)

    fig = Figure(size = (800, 200 * length(trange)))
    for (row, i) in enumerate(trange)
        islast = i == last(trange)
        ax = Axis(fig[row, 1];
            xlabel = islast ? (latex ? L"t" : "t") : "",
            ylabel = label(i, latex),
            xticklabelsvisible = islast,
        )
        lines!(ax, ts, [d[j][i] for j in r])
        if plot_title !== nothing && row == 1
            ax.title = plot_title
        end
        xlims!(ax, t[begin], t[end])
    end
    return fig
end

"""
    plot_energy_error(t, ΔH; nplot, nt, latex)
    plot_energy_error(sol::GeometricSolution; energy, nplot, nt, latex)

Plot the relative energy error `[H(t) - H(0)] / H(0)` as a function of time.

Either pass a precomputed error series (`t::TimeSeries`, `ΔH::DataSeries`), or a
`GeometricSolution`, in which case the error is computed from the solution's `:h`
invariant (override with the `energy` keyword). Returns a Makie `Figure`.
"""
function Diagnostics.plot_energy_error(t::Union{TimeSeries, ScalarDataSeries}, ΔH::DataSeries;
        nplot = 1, nt = :auto, latex = true)
    r = _steprange(t, nplot, nt)
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"[H(t) - H(0)] / H(0)" : "[H(t) - H(0)] / H(0)",
    )
    lines!(ax, [t[j] for j in r], [ΔH[j] for j in r])
    xlims!(ax, t[begin], t[end])
    return fig
end

function Diagnostics.plot_energy_error(sol::GeometricSolution; energy = nothing, kwargs...)
    Diagnostics.plot_energy_error(sol.t, _energy_error(sol; energy); kwargs...)
end

"""
    plot_energy_drift(t, d; nt, latex)

Scatter plot of the energy drift `ΔH` (maximum absolute energy error per interval,
see `GeometricSolutions.compute_drift`) as a function of time. Returns a `Figure`.
"""
function Diagnostics.plot_energy_drift(t::Union{TimeSeries, ScalarDataSeries}, d::DataSeries;
        nt = :auto, latex = true)
    # drift data is interval-based; the first entry (index 0) is not part of it.
    r = 1:(nt === :auto ? ntime(t) : min(nt, ntime(t)))
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"\Delta H" : "ΔH",
    )
    scatter!(ax, [t[j] for j in r], [d[j] for j in r])
    xlims!(ax, t[begin], t[end])
    return fig
end
function Diagnostics.plot_energy_drift(sol::GeometricSolution; energy = nothing, kwargs...)
    interval = div(ntime(sol), 10)
    ΔH = _energy_error(sol; energy)
    Diagnostics.plot_energy_drift(compute_error_drift(sol.t, ΔH, interval)...; kwargs...)
end

"""
    plot_constraint_error(sol::GeometricSolution; nplot, nt, k, latex, plot_title)
    plot_constraint_error(t, Δp; nplot, nt, k, latex, plot_title)

Plot the constraint (momentum) error `pᵢ(t) - ϑᵢ(t)` for an implicit ODE/DAE
solution, one stacked axis per degree of freedom (`k = 0`), or a single component
(`k = i`). Returns a Makie `Figure`.
"""
function Diagnostics.plot_constraint_error(t::Union{TimeSeries, ScalarDataSeries}, Δp::DataSeries; kwargs...)
    _plot_components(t, Δp,
        (i, latex) -> latex ? L"p_%$i(t) - \vartheta_%$i(t)" :
                              "p" * subscript(i) * "(t) - ϑ" * subscript(i) * "(t)";
        kwargs...)
end

function Diagnostics.plot_constraint_error(sol::GeometricSolution; kwargs...)
    Diagnostics.plot_constraint_error(sol.t, compute_momentum_error(sol); kwargs...)
end

"""
    plot_lagrange_multiplier(sol::GeometricSolution; nplot, nt, k, latex, plot_title)
    plot_lagrange_multiplier(t, λ; nplot, nt, k, latex, plot_title)

Plot the Lagrange multipliers `λᵢ(t)` of a DAE solution, one stacked axis per
component (`k = 0`) or a single component (`k = i`). Returns a Makie `Figure`.
"""
function Diagnostics.plot_lagrange_multiplier(t::Union{TimeSeries, ScalarDataSeries}, λ::DataSeries; kwargs...)
    _plot_components(t, λ,
        (i, latex) -> latex ? L"\lambda_%$i(t)" : "λ" * subscript(i) * "(t)";
        kwargs...)
end

function Diagnostics.plot_lagrange_multiplier(sol::GeometricSolution; kwargs...)
    Diagnostics.plot_lagrange_multiplier(sol.t, sol.λ; kwargs...)
end

end
