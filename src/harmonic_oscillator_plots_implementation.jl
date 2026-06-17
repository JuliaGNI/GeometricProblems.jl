# Implementation file for harmonic oscillator plotting
# This file is only included when CairoMakie is available via Requires.jl

using CairoMakie

# Spring geometry generation function
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
    plot_spring(i, q; fig=nothing, axfig=[1,1])

Plot a single frame of the harmonic oscillator spring-mass system.
"""
function plot_spring(i::Integer, q::Real; fig::Union{Figure,Nothing}=nothing, axfig::Vector=[1, 1])
    if fig === nothing
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
    plot_spring_animation(frames=0:100, output_dir="harmonic-oscillator", create_dir=true)

Create an animation of the spring-mass system.
"""
function plot_spring_animation(frames::AbstractRange=0:100,
                                output_dir::String="harmonic-oscillator",
                                create_dir::Bool=true)
    if create_dir && !isdir(output_dir)
        mkdir(output_dir)
    end

    for i in frames
        fig = plot_spring(i, zpos(i))
        save(joinpath(output_dir, "harmonic-oscillator-$(string(i, pad = 3)).png"), fig)
    end

    return nothing
end

"""
    plot_harmonic_oscillator(sol; ntime=nothing)

Plot the harmonic oscillator solution in phase space.
"""
function plot_harmonic_oscillator(sol; ntime::Union{Nothing,Integer}=nothing)
    # Handle different solution types
    if ntime === nothing
        ntime = hasproperty(sol, :t) ? length(sol.t) : length(sol.q)
    end

    # Extract position and velocity data
    q_data = hasproperty(sol, :q) ? sol.q : getproperty(sol, :q_values, [1.0])
    v_data = hasproperty(sol, :v) ? sol.v : getproperty(sol, :v_values, [0.0])

    # Create figure
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1,1], xlabel = "Position q", ylabel = "Velocity v",
             title = "Harmonic Oscillator Phase Space")

    # Plot trajectory
    lines!(ax, q_data[1:ntime], v_data[1:ntime], color = :blue, linewidth = 2)

    # Mark start and end points
    scatter!(ax, q_data[1], v_data[1], color = :red,
            markersize = 15, marker = :circle, label = "Start")
    scatter!(ax, q_data[ntime], v_data[ntime], color = :green,
            markersize = 15, marker = :circle, label = "End")

    # Add legend
    axislegend(ax, position = :lt)

    return fig
end
