module HarmonicOscillatorPlots

using Makie
import GeometricProblems.HarmonicOscillator

function _spring(x; rod = 0.15, top = 2.0, n = 100, w = 10, r = 0.1)
    distance = top - x - 2 * rod
    diststep = distance / n / w
    spring_x = [0.0, 0.0, [r * cos(2π * i / n) for i in 1:(n * w)]..., 0.0, 0.0]
    spring_y = [0.0, 0.0, [r * sin(2π * i / n) for i in 1:(n * w)]..., 0.0, 0.0]
    spring_z = [x, x + rod, [x + rod + diststep * i for i in 1:(n * w)]..., top - rod, top]
    return (spring_x, spring_y, spring_z)
end

"""
    plot_spring(t, x[, fig, axfig]; rod, top, n, w, r)

Plot a spring-mass visualization at position `x`. When called with only `t` and `x`,
creates a new standalone figure. When `fig` and `axfig` are provided the spring is
embedded into an existing figure at position `axfig`.

# Arguments
- `t`: current time (unused, kept for interface compatibility)
- `x`: current position of the mass
- `fig`: existing `Figure` to embed into (optional)
- `axfig`: grid position within `fig` for the 3D axis (optional)

# Keyword arguments
- `rod = 0.15`: length of rigid rod sections at top and bottom
- `top = 2.0`: fixed z-coordinate of the ceiling attachment point
- `n = 100`: number of coils in the spring
- `w = 10`: number of turns per coil
- `r = 0.1`: spring coil radius
"""
function HarmonicOscillator.plot_spring(
        t, x,
        fig = Figure(size = (100, 300), backgroundcolor = :transparent),
        axfig = fig[1, 1];
        rod = 0.15, top = 2.0, n = 100, w = 10, r = 0.1)
    ax3 = Axis3(
        axfig;
        aspect = (0.5, 0.5, 3),
        viewmode = :fitzoom,
        protrusions = 0,
        xspinesvisible = false,
        yspinesvisible = false,
        zspinesvisible = false,
    )
    xlims!(ax3, -0.25, +0.25)
    ylims!(ax3, -0.25, +0.25)
    zlims!(ax3, -1, +2)
    hidedecorations!(ax3)
    lines!(ax3, _spring(x; rod, top, n, w, r)..., color = :black)
    scatter!(ax3, Observable(Point3f[(0.0, 0.0, x)]), markersize = 25)
    return fig
end

"""
    plot_solution(i, t, q, p, h, labels[, hmod]; params, trange, rod, top)

Plot a snapshot of the harmonic oscillator solution at time index `i`. The figure
contains a time trace of the position, a time trace of the energy, a phase-space
portrait with Hamiltonian level sets, and a spring-mass cartoon.

# Arguments
- `i`: time index of the current frame (used for slicing `t[begin:i]`, etc.)
- `t`, `q`, `p`, `h`: time, position, momentum, and energy arrays (0- or 1-indexed)
- `labels`: named tuple `(t, q, p, h)` of axis label strings
- `hmod`: optional modified-Hamiltonian function `(q, p) -> value`

# Keyword arguments
- `params`: Hamiltonian parameters (default: `HarmonicOscillator.default_parameters()`)
- `trange = 100.0`: width of the sliding time window
- `rod = 0.15`: rigid rod length for the spring cartoon
- `top = 2.0`: ceiling z-coordinate for the spring cartoon
"""
function HarmonicOscillator.plot_solution(
        i, t, q, p, h, labels, hmod = nothing;
        params = HarmonicOscillator.default_parameters(),
        trange = 100.0,
        rod = 0.15, top = 2.0)
    nq = 600
    np = 600
    qs = range(-2.5, +2.5, nq)
    ps = range(-2.0, +2.0, np)
    hs = [HarmonicOscillator.hamiltonian(0, _q, _p, params) for _q in qs, _p in ps]

    if hmod !== nothing
        hmods = [hmod(_q, _p) for _q in qs, _p in ps]
        hmodt = hmod.(q, p)
    end

    levels = [0.0, 0.05, 0.25, 0.65, 1.2, collect(2:7)...]

    tlims = [
        max(0.0, t[i] - trange),
        t[end] < trange ? t[end] : max(t[i], trange),
    ]

    xlim = min(5, maximum(abs.(q)))
    hlim = ceil(max(1, min(5, maximum(h))))

    if xlim < 2
        xticks = [-1, 0, +1]
    elseif xlim < 4
        xticks = [-2, 0, +2]
    else
        xticks = [-4, -2, 0, +2, +4]
    end

    hticks = 0:Int(floor(max(1, min(5, maximum(h)))))

    fig = Figure(size = (1000, 400), figure_padding = 5, fontsize = 24)

    ax_solution = Axis(
        fig[1, 1:3];
        aspect = 3,
        ylabel = "$(labels.q)($(labels.t))",
        yticks = xticks,
        xticklabelsvisible = false,
    )
    xlims!(ax_solution, tlims...)
    ylims!(ax_solution, min(-1, -xlim), max(+1, +xlim))
    lines!(ax_solution, t[begin:i], q[begin:i]; linewidth = 2)

    ax_hamiltonian = Axis(
        fig[2, 1:3];
        aspect = 3,
        xlabel = "$(labels.t)",
        ylabel = "$(labels.h)($(labels.t))",
        yticks = hticks,
    )
    xlims!(ax_hamiltonian, tlims...)
    ylims!(ax_hamiltonian, 0, hlim)
    lines!(ax_hamiltonian, t[begin:i], h[begin:i]; linewidth = 2, label = "$(labels.h)($(labels.t))")

    if hmod !== nothing
        lines!(ax_hamiltonian, t[begin:i], hmodt[begin:i]; linewidth = 2, label = "H̄(t)")
        axislegend(ax_hamiltonian; position = :lt, orientation = :horizontal)
    end

    ax_phasespace = Axis(
        fig[1:2, 4:5];
        aspect = 1,
        xlabel = "$(labels.q)",
        ylabel = "$(labels.p)",
        xticks = [-2, -1, 0, +1, +2],
    )
    contourf!(ax_phasespace, qs, ps, hs; levels = levels)
    if hmod !== nothing
        contour!(ax_phasespace, qs, ps, hmods; levels = levels, color = :black, linewidth = 2)
    end
    xlims!(ax_phasespace, -2.5, +2.5)
    ylims!(ax_phasespace, -2.0, +2.0)
    lines!(ax_phasespace, q[begin:i], p[begin:i]; color = :red, linewidth = 2)

    HarmonicOscillator.plot_spring(t[i], q[i], fig, fig[1:2, 6]; rod, top)

    return fig
end

end
