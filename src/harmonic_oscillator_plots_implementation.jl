# Physical constants (redundant with main plot module but needed for standalone)
const k = 0.5        # Spring constant
const m = 1.0        # Mass
const ω = sqrt(k/m)  # Angular frequency

# Visualization parameters
const n = 100        # Number of coils in the spring
const w = 10         # Winding density
const r = 0.1        # Radius of spring
const h = 2π / ω / n  # Time step for visualization
const A = 1.0        # Amplitude

# Spring geometry
function xpos(i)
    sin(r*h*i)
end

function ypos(i)
    sin(r*h*i)
end

function zpos(i)
    A*cos(ω*h*i)
end

function ϑ(i)
    -ω*A*sin(ω*h*i)
end

# Spring dimensions
const rod = 0.15        # Length of rigid rods at spring ends
const top = zpos(0) + 1  # Top position

"""
Generate spring geometry for a given frame.

Arguments:
- `i`: Frame index

Returns:
- Tuple of (x, y, z) coordinates for the spring
"""
function spring(i)
    bottom = zpos(i)
    distance = top - bottom - 2 * rod
    diststep = distance / n / w
    spring_x = [0.0, 0.0, [r * cos(2π*i/n) for i in 1:n*w]..., 0.0, 0.0]
    spring_y = [0.0, 0.0, [r * sin(2π*i/n) for i in 1:n*w]..., 0.0, 0.0]
    spring_z = [bottom, bottom + rod, [bottom + rod + diststep * i for i in 1:n*w]..., top - rod, top]
    (spring_x, spring_y, spring_z)
end

"""
Plot a single frame of the spring-mass animation.

Arguments:
- `i`: Frame index
- `q`: Position of the mass
- `fig`: Optional figure (default: new figure)
- `axfig`: Optional axis position in the figure

Returns:
- The figure with the plot
"""
function plot_spring(i, q; fig = nothing, axfig = [1, 1])
    if isnothing(fig)
        fig = Figure(size = (100, 300), backgroundcolor = :transparent)
    end
    ax = fig[axfig...] = Axis3(fig,
        aspect = (0.5, 0.5, 3),
        viewmode = :fitzoom,
        protrusions = 0,
        xspinesvisible = false,
        yspinesvisible = false,
        zspinesvisible = false,
    )

    xlims!(ax, -0.25, +0.25)
    ylims!(ax, -0.25, +0.25)
    zlims!(ax, -1, +2)

    hidedecorations!(ax)

    lines!(ax, spring(i)..., color = :black)
    scatter!(ax, Point3f[(0., 0., q)], markersize = 25)

    return fig
end

"""
Plot the harmonic oscillator solution in phase space.

Arguments:
- `sol`: GeometricSolution from integrating the harmonic oscillator
- `ntime`: Number of time steps to plot (default: all)

Returns:
- The figure with the plot
"""
function plot_harmonic_oscillator(sol::GeometricSolution; ntime = length(sol.t))
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1,1], xlabel = "Position q", ylabel = "Velocity v", title = "Harmonic Oscillator Phase Space")

    lines!(ax, sol.q[1:ntime], sol.v[1:ntime], color = :blue, linewidth = 2)
    scatter!(ax, sol.q[1], sol.v[1], color = :red, markersize = 15, marker = :circle)
    scatter!(ax, sol.q[ntime], sol.v[ntime], color = :green, markersize = 15, marker = :circle)

    return fig
end

"""
Create an animation of the spring-mass system.

Arguments:
- `frames`: Range of frame indices to animate
- `output_dir`: Directory to save frames (default: "harmonic-oscillator")
- `create_dir`: Whether to create output directory if it doesn't exist

Returns:
- Nothing, but saves PNG files for each frame
"""
function plot_spring_animation(frames = 0:100, output_dir = "harmonic-oscillator", create_dir = true)
    # Create output directory if it doesn't exist
    if create_dir && !isdir(output_dir)
        mkdir(output_dir)
    end

    # Generate and save each frame
    for i in frames
        fig = plot_spring(i, zpos(i))
        save(joinpath(output_dir, "harmonic-oscillator-$(string(i, pad = 3)).png"), fig)
        destroy!(fig)  # Free memory
    end
end

"""
Create an animation from a harmonic oscillator solution.

Arguments:
- `sol`: GeometricSolution from harmonic oscillator
- `frames_to_plot`: Which time steps to animate (default: every 10th)
- `output_dir`: Directory to save frames (default: "harmonic-oscillator-sol")

Returns:
- Nothing, but saves PNG files for each frame
"""
function plot_harmonic_oscillator_animation(sol::GeometricSolution; frames_to_plot = 1:10:length(sol.t), output_dir = "harmonic-oscillator-sol", create_dir = true)
    # Create output directory if it doesn't exist
    if create_dir && !isdir(output_dir)
        mkdir(output_dir)
    end

    # Get scale factor for visualization
    max_amplitude = maximum(abs, sol.q)
    scale_factor = 1.0  # Could be parameterized

    # Generate and save each frame
    for (frame_idx, time_idx) in enumerate(frames_to_plot)
        q_raw = sol.q[time_idx]
        # Scale the position appropriately for visualization
        q_scaled = q_raw * scale_factor

        fig = plot_spring(frame_idx, q_scaled)
        save(joinpath(output_dir, "harmonic-oscillator-sol-$(string(frame_idx, pad = 3)).png"), fig)
        destroy!(fig)  # Free memory
    end
end
