module MasslessChargedParticlePlots

using Makie
using GeometricEquations: invariants, parameters
# Import GeometricSolutions symbols explicitly: `using Makie` also exports names
# such as `TimeSeries`, which would otherwise clash with GeometricSolutions.
using GeometricSolutions: ntime, compute_invariant_error

import GeometricProblems.MasslessChargedParticle

# Downsampled integer time indices 0:nplot:nt (nt = :auto → all stored steps).
_indices(sol, nplot, nt) = 0:nplot:(nt === :auto ? ntime(sol) : min(nt, ntime(sol)))

_energy_error(sol, equ) = compute_invariant_error(sol.t, sol.q, parameters(equ), invariants(equ)[:h])[2]

"""
    plot_phase_portrait(sol; nplot, nt, latex)

Plot the `x₁`–`x₂` trajectory of a massless charged particle. Returns a Makie
`Figure`.

# Arguments
- `sol <: GeometricSolution`

# Keyword arguments
- `nplot = 1`: plot every `nplot`-th time step
- `nt = :auto`: last time step to plot
- `latex = true`: use LaTeX axis labels
"""
function MasslessChargedParticle.plot_phase_portrait(sol; nplot = 1, nt = :auto, latex = true)
    idx = _indices(sol, nplot, nt)
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
    plot_solution(sol, equ; nplot, nt, latex)

Plot the `x₁`–`x₂` trajectory of a massless charged particle next to its relative
energy error. Returns a Makie `Figure`.
"""
function MasslessChargedParticle.plot_solution(sol, equ; nplot = 1, nt = :auto, latex = true)
    idx = _indices(sol, nplot, nt)
    ΔH  = _energy_error(sol, equ)

    fig = Figure(size = (800, 300))

    ax_phase = Axis(fig[1, 1];
        aspect = 1,
        xlabel = latex ? L"x_1" : "x₁",
        ylabel = latex ? L"x_2" : "x₂",
    )
    lines!(ax_phase, [sol.q[k][1] for k in idx], [sol.q[k][2] for k in idx])

    ax_energy = Axis(fig[1, 2];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"[H(t) - H(0)] / H(0)" : "[H(t) - H(0)] / H(0)",
    )
    lines!(ax_energy, [sol.t[k] for k in idx], [ΔH[k] for k in idx])

    return fig
end

"""
    plot_traces(sol, equ; nplot, nt, latex)

Plot the time traces `x₁(t)`, `x₂(t)` of a massless charged particle trajectory
together with its relative energy error, stacked vertically. Returns a Makie
`Figure`.
"""
function MasslessChargedParticle.plot_traces(sol, equ; nplot = 1, nt = :auto, latex = true)
    idx = _indices(sol, nplot, nt)
    ts  = [sol.t[k] for k in idx]
    ΔH  = _energy_error(sol, equ)

    ylabels = latex ? (L"x_1", L"x_2") : ("x₁", "x₂")

    fig = Figure(size = (800, 600))
    for i in 1:2
        ax = Axis(fig[i, 1]; ylabel = ylabels[i], xticklabelsvisible = false)
        lines!(ax, ts, [sol.q[k][i] for k in idx])
    end
    ax_energy = Axis(fig[3, 1];
        xlabel = latex ? L"t" : "t",
        ylabel = latex ? L"[H(t) - H(0)] / H(0)" : "[H(t) - H(0)] / H(0)",
    )
    lines!(ax_energy, ts, [ΔH[k] for k in idx])
    return fig
end

end
