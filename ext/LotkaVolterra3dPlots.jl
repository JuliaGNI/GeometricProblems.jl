module LotkaVolterra3dPlots

using Makie
# Import GeometricSolutions symbols explicitly: `using Makie` also exports names
# such as `TimeSeries`, which would otherwise clash with GeometricSolutions.
using GeometricSolutions: ntime, compute_invariant_error

import GeometricProblems.LotkaVolterra3d

"""
    plot_phase_portrait(sol; nplot, nt, latex)

Plot the 3D trajectory `(x₁, x₂, x₃)` of a 3D Lotka-Volterra solution on the
Casimir surface. Returns a Makie `Figure`.

# Arguments
- `sol <: GeometricSolution`

# Keyword arguments
- `nplot = 1`: plot every `nplot`-th time step
- `nt = :auto`: last time step to plot
- `latex = true`: use LaTeX axis labels
"""
function LotkaVolterra3d.plot_phase_portrait(sol; nplot = 1, nt = :auto, latex = true)
    idx = 0:nplot:(nt === :auto ? ntime(sol) : min(nt, ntime(sol)))
    fig = Figure(size = (600, 600))
    ax = Axis3(fig[1, 1];
        xlabel = latex ? L"x_1" : "x₁",
        ylabel = latex ? L"x_2" : "x₂",
        zlabel = latex ? L"x_3" : "x₃",
    )
    lines!(ax, [sol.q[k][1] for k in idx], [sol.q[k][2] for k in idx], [sol.q[k][3] for k in idx])
    return fig
end

"""
    plot_traces(sol, params; nplot, nt, latex)

Plot the time traces `x₁(t)`, `x₂(t)`, `x₃(t)` of a 3D Lotka-Volterra solution
together with its relative energy error, stacked vertically. Returns a Makie
`Figure`.

# Arguments
- `sol <: GeometricSolution`
- `params`: the Hamiltonian parameters (named tuple)

# Keyword arguments
- `nplot = 1`: plot every `nplot`-th time step
- `nt = :auto`: last time step to plot
- `latex = true`: use LaTeX axis labels
"""
function LotkaVolterra3d.plot_traces(sol, params; nplot = 1, nt = :auto, latex = true)
    idx = 0:nplot:(nt === :auto ? ntime(sol) : min(nt, ntime(sol)))
    ts  = [sol.t[k] for k in idx]
    _, ΔH = compute_invariant_error(sol.t, sol.q, params, LotkaVolterra3d.hamiltonian)

    ylabels = latex ? (L"x_1", L"x_2", L"x_3") : ("x₁", "x₂", "x₃")

    fig = Figure(size = (800, 1000))
    for i in 1:3
        ax = Axis(fig[i, 1]; ylabel = ylabels[i], xticklabelsvisible = false)
        lines!(ax, ts, [sol.q[k][i] for k in idx])
        xlims!(ax, ts[begin], ts[end])
    end
    ax_energy = Axis(fig[4, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"[H(t) - H(0)] / H(0)" : "[H(t) - H(0)] / H(0)",
    )
    lines!(ax_energy, ts, [ΔH[k] for k in idx])
    xlims!(ax_energy, ts[begin], ts[end])
    return fig
end

end
