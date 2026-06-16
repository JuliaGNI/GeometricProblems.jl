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
function plot_harmonic_oscillator(sol; ntime = hasproperty(sol, :t) ? length(sol.t) : length(sol.q))
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1,1], xlabel = "Position q", ylabel = "Velocity v", title = "Harmonic Oscillator Phase Space")

    q_data = hasproperty(sol, :q) ? sol.q : sol.q_values
    v_data = hasproperty(sol, :v) ? sol.v : sol.v_values

    lines!(ax, q_data[1:ntime], v_data[1:ntime], color = :blue, linewidth = 2)
    scatter!(ax, q_data[1], v_data[1], color = :red, markersize = 15, marker = :circle)
    scatter!(ax, q_data[ntime], v_data[ntime], color = :green, markersize = 15, marker = :circle)

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
function plot_harmonic_oscillator_animation(sol; frames_to_plot = nothing, output_dir = "harmonic-oscillator-sol", create_dir = true)
    # Create output directory if it doesn't exist
    if create_dir && !isdir(output_dir)
        mkdir(output_dir)
    end

    # Handle different solution types
    if hasproperty(sol, :t)
        t_data = sol.t
    else
        t_data = getproperty(sol, :time_values, 1:length(getproperty(sol, :q_values, [1.0])))
    end

    if hasproperty(sol, :q)
        q_data = sol.q
    else
        q_data = getproperty(sol, :q_values, zeros(length(t_data)))
    end

    # Default frames: every 10th point if not specified
    if isnothing(frames_to_plot)
        frames_to_plot = 1:10:length(t_data)
    elseif frames_to_plot isa UnitRange{Int}
        # If it's a range of indices, use as-is
    else
        # Assume it's a range of times, find nearest indices
        frames_to_plot = map(frames_to_plot) do time_val
            argmin(abs.(t_data .- time_val))
        end
    end

    # Get scale factor for visualization
    max_amplitude = maximum(abs, q_data)
    scale_factor = 1.0  # Could be parameterized

    # Generate and save each frame
    for (frame_idx, time_idx) in enumerate(frames_to_plot)
        if time_idx <= length(q_data) && time_idx >= 1
            q_raw = q_data[time_idx]
            # Scale the position appropriately for visualization
            q_scaled = q_raw * scale_factor

            fig = plot_spring(frame_idx, q_scaled)
            save(joinpath(output_dir, "harmonic-oscillator-sol-$(string(frame_idx, pad = 3)).png"), fig)
            destroy!(fig)  # Free memory
        end
    end
end
