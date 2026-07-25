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

# Compute the relative error of one of a solution's invariants. `invariant` is either a key
# into `invariants(sol.problem)` or the invariant function itself.
function _invariant_error(sol, invariant)
    f = invariant isa Symbol ? invariants(sol.problem)[invariant] : invariant
    params = parameters(sol.problem)
    if sol isa Union{SolutionPODE, SolutionPDAE}
        _, Δ = compute_invariant_error(sol.t, sol.q, sol.p, params, f)
    else
        _, Δ = compute_invariant_error(sol.t, sol.q, params, f)
    end
    return Δ
end

# Compute the relative energy error of a solution from its `:h` invariant.
_energy_error(sol; energy = nothing) = _invariant_error(sol, energy === nothing ? :h : energy)

# Default axis labels for the error and the drift of an invariant named `s`, e.g. `:h` → "H".
# A non-`Symbol` invariant (a bare function) has no name to show, so it falls back to "I".
_invariant_name(s::Symbol) = uppercase(string(s))
_invariant_name(_) = "I"

function _error_label(s, latex)
    n = _invariant_name(s)
    latex ? L"[%$n(t) - %$n(0)] / %$n(0)" : "[$n(t) - $n(0)] / $n(0)"
end

function _drift_label(s, latex)
    n = _invariant_name(s)
    latex ? L"\Delta %$n" : "Δ$n"
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
        xlims!(ax, ts[begin], ts[end])
    end
    return fig
end

"""
    plot_invariant_error(t, Δ; ylabel, nplot, nt, latex)
    plot_invariant_error(sol::GeometricSolution; invariant, ylabel, nplot, nt, latex)

Plot the relative error `[I(t) - I(0)] / I(0)` of an invariant as a function of time.

Either pass a precomputed error series (`t::TimeSeries`, `Δ::DataSeries`), or a
`GeometricSolution`, in which case the error is computed from the invariant named by the
`invariant` keyword (a key into `invariants(sol.problem)`, `:h` by default; a callable is
also accepted). The axis label is derived from that name and can be overridden with
`ylabel`. Returns a Makie `Figure`.
"""
function Diagnostics.plot_invariant_error(t::Union{TimeSeries, ScalarDataSeries}, Δ::DataSeries;
        ylabel = nothing, invariant = :h, nplot = 1, nt = :auto, latex = true)
    r  = _steprange(t, nplot, nt)
    ts = [t[j] for j in r]
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = ylabel === nothing ? _error_label(invariant, latex) : ylabel,
    )
    lines!(ax, ts, [Δ[j] for j in r])
    xlims!(ax, ts[begin], ts[end])
    return fig
end

function Diagnostics.plot_invariant_error(sol::GeometricSolution; invariant = :h, kwargs...)
    Diagnostics.plot_invariant_error(sol.t, _invariant_error(sol, invariant); invariant, kwargs...)
end

"""
    plot_invariant_drift(t, d; ylabel, nt, latex)
    plot_invariant_drift(sol::GeometricSolution; invariant, ylabel, nt, latex)

Scatter plot of the drift of an invariant (the maximum absolute error per interval, see
`GeometricSolutions.compute_error_drift`) as a function of time. Keywords are those of
[`plot_invariant_error`](@ref). Returns a `Figure`.
"""
function Diagnostics.plot_invariant_drift(t::Union{TimeSeries, ScalarDataSeries}, d::DataSeries;
        ylabel = nothing, invariant = :h, nt = :auto, latex = true)
    # drift data is interval-based; the first entry (index 0) is not part of it.
    r  = 1:(nt === :auto ? ntime(t) : min(nt, ntime(t)))
    ts = [t[j] for j in r]
    fig = Figure(size = (800, 400))
    ax = Axis(fig[1, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = ylabel === nothing ? _drift_label(invariant, latex) : ylabel,
    )
    scatter!(ax, ts, [d[j] for j in r])
    xlims!(ax, ts[begin], ts[end])
    return fig
end

function Diagnostics.plot_invariant_drift(sol::GeometricSolution; invariant = :h, kwargs...)
    interval = div(ntime(sol), 10)
    Δ = _invariant_error(sol, invariant)
    Diagnostics.plot_invariant_drift(compute_error_drift(sol.t, Δ, interval)...; invariant, kwargs...)
end

"""
    plot_energy_error(t, ΔH; nplot, nt, latex)
    plot_energy_error(sol::GeometricSolution; energy, nplot, nt, latex)

Plot the relative energy error `[H(t) - H(0)] / H(0)` as a function of time.

Either pass a precomputed error series (`t::TimeSeries`, `ΔH::DataSeries`), or a
`GeometricSolution`, in which case the error is computed from the solution's `:h`
invariant (override with the `energy` keyword). Returns a Makie `Figure`.

This is [`plot_invariant_error`](@ref) for the `:h` invariant.
"""
function Diagnostics.plot_energy_error(t::Union{TimeSeries, ScalarDataSeries}, ΔH::DataSeries; kwargs...)
    Diagnostics.plot_invariant_error(t, ΔH; invariant = :h, kwargs...)
end

function Diagnostics.plot_energy_error(sol::GeometricSolution; energy = nothing, kwargs...)
    Diagnostics.plot_energy_error(sol.t, _energy_error(sol; energy); kwargs...)
end

"""
    plot_energy_drift(t, d; nt, latex)

Scatter plot of the energy drift `ΔH` (maximum absolute energy error per interval,
see `GeometricSolutions.compute_error_drift`) as a function of time. Returns a `Figure`.

This is [`plot_invariant_drift`](@ref) for the `:h` invariant.
"""
function Diagnostics.plot_energy_drift(t::Union{TimeSeries, ScalarDataSeries}, d::DataSeries; kwargs...)
    Diagnostics.plot_invariant_drift(t, d; invariant = :h, kwargs...)
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

"""
    plot_convergence(h, ε; order, ylabel, latex)

Log-log plot of an error `ε` against the time step `h`, for a convergence study over a
sequence of time steps. Passing `order = p` adds a reference line of slope `p`, anchored
one factor of two above the first data point, so that the observed slope can be compared
against the expected order of the method. Returns a Makie `Figure`.

Data points with a vanishing error are dropped: they carry no information on a logarithmic
axis and would make the axis limits degenerate. Methods that are exact for the diagnostic
in question can therefore produce an empty plot.
"""
function Diagnostics.plot_convergence(h::AbstractVector, ε::AbstractVector;
        order = nothing, ylabel = nothing, latex = true)
    length(h) == length(ε) || throw(ArgumentError("h and ε must have the same length"))

    keep = findall(>(0), abs.(ε))

    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1];
        xscale = log10,
        yscale = log10,
        xlabel = latex ? L"h" : "h",
        ylabel = ylabel === nothing ? (latex ? L"\varepsilon" : "ε") : ylabel,
    )

    # The reference slope is anchored at twice the first *kept* error, matching the offset
    # the study in GeometricExamples has always used, and spans the full range of h.
    if order !== nothing && !isempty(keep)
        offset = 2 * abs(ε[first(keep)])
        h₀ = h[first(keep)]
        lines!(ax, h, [offset * (hᵢ / h₀)^order for hᵢ in h];
            linestyle = :dash, label = latex ? L"O(h^{%$order})" : "O(h^$order)")
    end

    if !isempty(keep)
        scatterlines!(ax, h[keep], abs.(ε[keep]); label = latex ? L"\varepsilon" : "ε")
        order === nothing || axislegend(ax; position = :rb)
    end

    return fig
end

"""
    plot_order(h, p; latex)

Plot the observed order of convergence `p` (typically `log2(ε(h) / ε(h/2))`) against the
time step `h` on a logarithmic `h` axis. Returns a Makie `Figure`.
"""
function Diagnostics.plot_order(h::AbstractVector, p::AbstractVector; latex = true)
    length(h) == length(p) || throw(ArgumentError("h and p must have the same length"))

    fig = Figure(size = (600, 600))
    ax = Axis(fig[1, 1];
        xscale = log10,
        xlabel = latex ? L"h" : "h",
        ylabel = latex ? L"p" : "p",
    )
    scatterlines!(ax, h, p)
    return fig
end

end
