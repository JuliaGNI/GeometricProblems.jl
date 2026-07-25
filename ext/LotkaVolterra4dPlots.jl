module LotkaVolterra4dPlots

using Makie
using GeometricEquations: invariants, parameters
# Import GeometricSolutions symbols explicitly: `using Makie` also exports names
# such as `TimeSeries`, which would otherwise clash with GeometricSolutions.
using GeometricSolutions: ntime, compute_invariant_error

import GeometricProblems.LotkaVolterra4d

_energy_error(sol, equ) = compute_invariant_error(sol.t, sol.q, parameters(equ), invariants(equ)[:h])[2]

"""
    plot_phase_portrait(sol; nplot, nt, latex)

Plot the projection of a 4D Lotka-Volterra trajectory onto the `(x₁, x₂)` plane.
Returns a Makie `Figure`.

# Arguments
- `sol <: GeometricSolution`

# Keyword arguments
- `nplot = 1`: plot every `nplot`-th time step
- `nt = :auto`: last time step to plot
- `latex = true`: use LaTeX axis labels
"""
function LotkaVolterra4d.plot_phase_portrait(sol; nplot = 1, nt = :auto, latex = true)
    idx = 0:nplot:(nt === :auto ? ntime(sol) : min(nt, ntime(sol)))
    fig = Figure(size = (400, 400))
    ax = Axis(fig[1, 1];
        aspect = 1,
        xlabel = latex ? L"x_1" : "x₁",
        ylabel = latex ? L"x_2" : "x₂",
    )
    lines!(ax, [sol.q[k][1] for k in idx], [sol.q[k][2] for k in idx])
    return fig
end

"""
    plot_traces(sol, equ; nplot, nt, latex)

Plot the time traces `x₁(t)`, …, `x₄(t)` of a 4D Lotka-Volterra solution together
with its relative energy error, stacked vertically. Returns a Makie `Figure`.

# Arguments
- `sol <: GeometricSolution`
- `equ`: the problem `sol` was obtained from, which supplies the `:h` invariant and its parameters

# Keyword arguments
- `nplot = 1`: plot every `nplot`-th time step
- `nt = :auto`: last time step to plot
- `latex = true`: use LaTeX axis labels
"""
function LotkaVolterra4d.plot_traces(sol, equ; nplot = 1, nt = :auto, latex = true)
    idx = 0:nplot:(nt === :auto ? ntime(sol) : min(nt, ntime(sol)))
    ts  = [sol.t[k] for k in idx]
    ΔH  = _energy_error(sol, equ)

    ylabels = latex ? (L"x_1", L"x_2", L"x_3", L"x_4") : ("x₁", "x₂", "x₃", "x₄")

    fig = Figure(size = (800, 1200))
    for i in 1:4
        ax = Axis(fig[i, 1]; ylabel = ylabels[i], xticklabelsvisible = false)
        lines!(ax, ts, [sol.q[k][i] for k in idx])
        xlims!(ax, ts[begin], ts[end])
    end
    ax_energy = Axis(fig[5, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"[H(t) - H(0)] / H(0)" : "[H(t) - H(0)] / H(0)",
    )
    lines!(ax_energy, ts, [ΔH[k] for k in idx])
    xlims!(ax_energy, ts[begin], ts[end])
    return fig
end

end
