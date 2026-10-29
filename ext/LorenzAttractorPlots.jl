module LorenzAttractorPlots

using Makie
# Import GeometricSolutions symbols explicitly: `using Makie` also exports names
# such as `TimeSeries`, which would otherwise clash with GeometricSolutions.
using GeometricSolutions: ntime

import GeometricProblems.LorenzAttractor

"""
    plot_phase_portrait(sol; nplot, nt, latex)

Plot the 3D trajectory `(x, y, z)` of a Lorenz attractor solution. Returns a Makie
`Figure`.

# Arguments
- `sol <: GeometricSolution`

# Keyword arguments
- `nplot = 1`: plot every `nplot`-th time step
- `nt = :auto`: last time step to plot
- `latex = true`: use LaTeX axis labels
"""
function LorenzAttractor.plot_phase_portrait(sol; nplot = 1, nt = :auto, latex = true)
    idx = 0:nplot:(nt === :auto ? ntime(sol) : min(nt, ntime(sol)))
    fig = Figure(size = (600, 600))
    ax = Axis3(fig[1, 1];
        xlabel = latex ? L"x" : "x",
        ylabel = latex ? L"y" : "y",
        zlabel = latex ? L"z" : "z"
    )
    lines!(ax, [sol.q[k][1] for k in idx], [sol.q[k][2] for k in idx],
        [sol.q[k][3] for k in idx]; linewidth = 0.6)
    return fig
end

end
