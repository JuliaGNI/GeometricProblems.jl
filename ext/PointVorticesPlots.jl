module PointVorticesPlots

using Makie
using GeometricEquations: invariants, parameters
# Import GeometricSolutions symbols explicitly: `using Makie` also exports names
# such as `TimeSeries`, which would otherwise clash with GeometricSolutions.
using GeometricSolutions: ntime, compute_invariant_error

import GeometricProblems.PointVortices
import GeometricProblems.PointVorticesLinear

# Downsampled integer time indices 0:nplot:nt (nt = :auto → all stored steps).
_indices(sol, nplot, nt) = 0:nplot:(nt === :auto ? ntime(sol) : min(nt, ntime(sol)))

_energy_error(sol, equ) = compute_invariant_error(sol.t, sol.q, parameters(equ), invariants(equ)[:h])[2]

# Both vortices live in the same plane, so the state (q₁, q₂, q₃, q₄) is read as the two
# positions (q₁, q₂) and (q₃, q₄). `_vortex(sol, idx, k)` returns the x/y coordinate pair of
# vortex `k`.
function _vortex(sol, idx, k)
    i = 2k - 1
    ([sol.q[j][i] for j in idx], [sol.q[j][i+1] for j in idx])
end

# The plot functions depend only on the solution `sol`, the equation `equ` (for the `:h`
# invariant) and the keyword arguments, so the same implementations serve the deformed
# `PointVortices` model and the linearised `PointVorticesLinear` one.

function _plot_phase_portrait(sol; nplot = 1, nt = :auto, latex = true)
    idx = _indices(sol, nplot, nt)
    fig = Figure(size = (400, 400))
    ax = Axis(fig[1, 1];
        aspect = 1,
        xlabel = latex ? L"x" : "x",
        ylabel = latex ? L"y" : "y",
    )
    for k in 1:2
        x, y = _vortex(sol, idx, k)
        lines!(ax, x, y; label = latex ? L"\gamma_%$k" : "γ" * string(k))
    end
    axislegend(ax)
    return fig
end

function _plot_solution(sol, equ; nplot = 1, nt = :auto, latex = true)
    idx = _indices(sol, nplot, nt)
    ts  = [sol.t[j] for j in idx]
    ΔH  = _energy_error(sol, equ)

    fig = Figure(size = (800, 300))

    ax_phase = Axis(fig[1, 1];
        aspect = 1,
        xlabel = latex ? L"x" : "x",
        ylabel = latex ? L"y" : "y",
    )
    for k in 1:2
        x, y = _vortex(sol, idx, k)
        lines!(ax_phase, x, y)
    end

    ax_energy = Axis(fig[1, 2];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"[H(t) - H(0)] / H(0)" : "[H(t) - H(0)] / H(0)",
    )
    lines!(ax_energy, ts, [ΔH[j] for j in idx])
    xlims!(ax_energy, ts[begin], ts[end])

    return fig
end

function _plot_traces(sol, equ; nplot = 1, nt = :auto, latex = true)
    idx = _indices(sol, nplot, nt)
    ts  = [sol.t[j] for j in idx]
    ΔH  = _energy_error(sol, equ)

    # (x, y) of the first vortex, then of the second one.
    ylabels = latex ? (L"x_1", L"y_1", L"x_2", L"y_2") : ("x₁", "y₁", "x₂", "y₂")

    fig = Figure(size = (800, 1000))
    for i in 1:4
        ax = Axis(fig[i, 1]; ylabel = ylabels[i], xticklabelsvisible = false)
        lines!(ax, ts, [sol.q[j][i] for j in idx])
        xlims!(ax, ts[begin], ts[end])
    end
    ax_energy = Axis(fig[5, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"[H(t) - H(0)] / H(0)" : "[H(t) - H(0)] / H(0)",
    )
    lines!(ax_energy, ts, [ΔH[j] for j in idx])
    xlims!(ax_energy, ts[begin], ts[end])
    return fig
end

# Attach the generic implementations to both problems' plot-function stubs.
@doc """
    plot_phase_portrait(sol; nplot, nt, latex)

Plot the trajectories of both point vortices in the `x`–`y` plane. Returns a Makie
`Figure`.

# Arguments
- `sol <: GeometricSolution`

# Keyword arguments
- `nplot = 1`: plot every `nplot`-th time step
- `nt = :auto`: last time step to plot
- `latex = true`: use LaTeX axis labels
"""
PointVortices.plot_phase_portrait(sol; kwargs...) = _plot_phase_portrait(sol; kwargs...)
@doc """
    plot_solution(sol, equ; nplot, nt, latex)

Plot the trajectories of both point vortices next to the relative energy error of the
solution. Returns a Makie `Figure`.
"""
PointVortices.plot_solution(sol, equ; kwargs...) = _plot_solution(sol, equ; kwargs...)
@doc """
    plot_traces(sol, equ; nplot, nt, latex)

Plot the time traces of the four state components together with the relative energy error,
stacked vertically. Returns a Makie `Figure`.
"""
PointVortices.plot_traces(sol, equ; kwargs...) = _plot_traces(sol, equ; kwargs...)

# The same implementations, and hence the same documentation as above.
PointVorticesLinear.plot_phase_portrait(sol; kwargs...) = _plot_phase_portrait(sol; kwargs...)
PointVorticesLinear.plot_solution(sol, equ; kwargs...) = _plot_solution(sol, equ; kwargs...)
PointVorticesLinear.plot_traces(sol, equ; kwargs...) = _plot_traces(sol, equ; kwargs...)

end
