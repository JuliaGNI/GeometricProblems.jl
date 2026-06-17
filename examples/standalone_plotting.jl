#!/usr/bin/env julia
#
# Standalone Harmonic Oscillator Plotting Script
# This works independently of the main GeometricProblems.jl package structure
# To avoid the complex Julia module integration issues I encountered.

println("🎨 Harmonic Oscillator Standalone Plotting Demo")
println("=" ^ 50)

# Ensure CairoMakie is available
try
    using CairoMakie
    println("✅ CairoMakie loaded")
catch
    println("📦 Installing CairoMakie...")
    using Pkg
    Pkg.add("CairoMakie")
    using CairoMakie
end

# Load harmonic oscillator module (if available)
try
    using GeometricProblems.HarmonicOscillator: odeproblem
    println("✅ GeometricProblems.HarmonicOscillator loaded")
catch
    println("⚠️  GeometricProblems.HarmonicOscillator not available")
    println("   Will work with standalone data")
end

# Define the plotting functions directly (no module integration issues)
#
# Physical constants (matching the original script)
const k = 0.5        # Spring constant
const m = 1.0        # Mass
const ω = sqrt(k/m)  # Angular frequency

# Visualization parameters
const n = 100        # Number of coils in the spring
const w = 10         # Winding density
const r = 0.1        # Radius of spring
const h = 2π / ω / n  # Time step for visualization
const A = 1.0        # Amplitude

# Spring geometry functions
xpos(i) = sin(r*h*i)
ypos(i) = sin(r*h*i)
zpos(i) = A*cos(ω*h*i)
ϑ(i) = -ω*A*sin(ω*h*i)

# Spring dimensions
const rod = 0.15        # Length of rigid rods at spring ends
const top = zpos(0) + 1  # Top position

# Spring geometry generation
function spring(i)
    bottom = zpos(i)
    distance = top - bottom - 2 * rod
    diststep = distance / n / w
    spring_x = [0.0, 0.0, [r * cos(2π*i/n) for i in 1:n*w]..., 0.0, 0.0]
    spring_y = [0.0, 0.0, [r * sin(2π*i/n) for i in 1:n*w]..., 0.0, 0.0]
    spring_z = [bottom, bottom + rod, [bottom + rod + diststep * i for i in 1:n*w]..., top - rod, top]
    (spring_x, spring_y, spring_z)
end

# Single frame plotting
function plot_spring(i, q; fig = nothing, axfig = [1, 1])
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

# Animation creation
function plot_spring_animation(frames = 0:100, output_dir = "harmonic-oscillator", create_dir = true)
    if create_dir && !isdir(output_dir)
        mkdir(output_dir)
    end

    for i in frames
        fig = plot_spring(i, zpos(i))
        save(joinpath(output_dir, "harmonic-oscillator-$(string(i, pad = 3)).png"), fig)
    end
end

# Phase space plotting
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

println("🎯 Creating demo outputs...")
if !isdir("standalone_demo")
    mkdir("standalone_demo")
end

# 1. Synthetic animation demo
println("1️⃣ Creating synthetic spring animation...")
plot_spring_animation(0:5:10, "standalone_demo/spring", true)
println("   ✅ Created spring animation frames")

# 2. Mock data phase space demo
println("2️⃣ Creating phase space plot...")
time_steps = 0:0.5:5
positions = cos.(time_steps)
velocities = -sin.(time_steps)

mock_sol = (q = positions, v = velocities, t = time_steps)
fig = plot_harmonic_oscillator(mock_sol)
save("standalone_demo/phase_space.png", fig)
println("   ✅ Created phase space plot")

# 3. Test figure
println("3️⃣ Creating basic test plot...")
test_fig = Figure()
test_ax = Axis(test_fig[1,1], title = "Basic Plot Test")
lines!(test_ax, 0..10, sin, color = :purple)
save("standalone_demo/test_plot.png", test_fig)
println("   ✅ Created test plot")

println("\n🎉 Demo completed successfully!")
println("📁 Output files:")
println("   - Spring animation frames: standalone_demo/spring/")
println("   - Phase space plot: standalone_demo/phase_space.png")
println("   - Test plot: standalone_demo/test_plot.png")
println("\n💡 To create movies:")
println("   ffmpeg -framerate 10 -i standalone_demo/spring/harmonic-oscillator-%03d.png spring.mp4")
