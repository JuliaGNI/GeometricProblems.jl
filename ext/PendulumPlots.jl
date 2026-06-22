module PendulumPlots

using Makie
import GeometricProblems.Pendulum

"""
    plot_pendulum(t, q[, fig, axfig])

Plot a pendulum visualization at angle `q`. When called with only `t` and `q`,
creates a new standalone figure. When `fig` and `axfig` are provided the pendulum
is embedded into an existing figure at position `axfig`.

# Arguments
- `t`: current time (unused, kept for interface compatibility)
- `q`: current angle of the pendulum (rad; q=0 is inverted, q=π is hanging)
- `fig`: existing `Figure` to embed into (optional)
- `axfig`: grid position within `fig` for the 2D axis (optional)
"""
function Pendulum.plot_pendulum(
        t, q,
        fig = Figure(size = (200, 200), backgroundcolor = :transparent),
        axfig = fig[1, 1])
    ax = Axis(
        axfig;
        aspect = 1,
        limits = (-1.5, 1.5, -1.5, 1.5),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    pivot = Point2f(0, 0)
    bob = Point2f(sin(q), cos(q))

    # Horizontal support bar
    lines!(ax, [-0.3, 0.3], [0.0, 0.0]; color = :darkgray, linewidth = 4)
    # Rod
    lines!(ax, [pivot, bob]; color = :black, linewidth = 2)
    # Pivot
    scatter!(ax, [pivot]; markersize = 10, color = :darkgray)
    # Bob
    scatter!(ax, [bob]; markersize = 25, color = :steelblue)

    return fig
end

"""
    plot_solution(i, t, q, p, h, labels[, hmod]; trange)

Plot a snapshot of the pendulum solution at time index `i`. The figure contains
a time trace of the angle, a time trace of the energy, a phase-space portrait
with Hamiltonian level sets, and a pendulum cartoon.

# Arguments
- `i`: time index of the current frame (used for slicing `t[begin:i]`, etc.)
- `t`, `q`, `p`, `h`: time, angle, momentum, and energy arrays
- `labels`: named tuple `(t, q, p, h)` of axis label strings
- `hmod`: optional modified-Hamiltonian function `(q, p) -> value`

# Keyword arguments
- `trange = 100.0`: width of the sliding time window
"""
function Pendulum.plot_solution(
        i, t, q, p, h, labels, hmod = nothing;
        trange = 100.0)
    nq = 600
    np = 600
    qs = range(0.0, 2π, nq)
    ps = range(-3.0, 3.0, np)
    hs = [Pendulum.hamiltonian(0, _q, _p) for _q in qs, _p in ps]

    if hmod !== nothing
        hmods = [hmod(_q, _p) for _q in qs, _p in ps]
        hmodt = hmod.(q, p)
    end

    # Levels covering the full range; -1.1 is below the global minimum (-1) so the
    # entire domain is filled. The separatrix H = 1 separates libration from rotation.
    levels = vcat([-1.1], collect(-0.9:0.3:0.9), [1.0, 2.0, 6.0])

    tlims = [
        max(0.0, t[i] - trange),
        t[end] < trange ? t[end] : max(t[i], trange),
    ]

    hmin = floor(min(0.0, minimum(h)))
    hlim = ceil(max(1, min(5, maximum(h))))

    fig = Figure(size = (1000, 400), figure_padding = 5, fontsize = 24)

    ax_solution = Axis(
        fig[1, 1:3];
        aspect = 3,
        ylabel = "$(labels.q)($(labels.t))",
        xticklabelsvisible = false,
    )
    xlims!(ax_solution, tlims...)
    lines!(ax_solution, t[begin:i], q[begin:i]; linewidth = 2)

    ax_hamiltonian = Axis(
        fig[2, 1:3];
        aspect = 3,
        xlabel = "$(labels.t)",
        ylabel = "$(labels.h)($(labels.t))",
    )
    xlims!(ax_hamiltonian, tlims...)
    ylims!(ax_hamiltonian, hmin, hlim)
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
        xticks = ([0, π, 2π], ["0", "π", "2π"]),
    )
    contourf!(ax_phasespace, qs, ps, hs; levels)
    if hmod !== nothing
        contour!(ax_phasespace, qs, ps, hmods; levels, color = :black, linewidth = 2)
    end
    # Highlight the separatrix
    contour!(ax_phasespace, qs, ps, hs; levels = [1.0], color = :white, linewidth = 2)
    xlims!(ax_phasespace, 0.0, 2π)
    ylims!(ax_phasespace, -3.0, 3.0)

    q_plot = mod.(q[begin:i], 2π)
    lines!(ax_phasespace, q_plot, p[begin:i]; color = :red, linewidth = 2)

    Pendulum.plot_pendulum(t[i], q[i], fig, fig[1:2, 6])

    return fig
end

end
