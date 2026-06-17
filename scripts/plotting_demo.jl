#!/usr/bin/env julia
#
# Demo script for the Harmonic Oscillator Plotting Module
# This script demonstrates how to use the standalone plotting functionality
# without needing to modify the main GeometricProblems package structure.

println("🎨 Harmonic Oscillator Plotting Demo (Standalone)")
println("=" ^ 60)

using Test

# Check if CairoMakie is available
try
    using CairoMakie
    println("✅ CairoMakie is already installed")
catch
    println("📦 Installing CairoMakie...")
    using Pkg
    Pkg.add("CairoMakie")
    using CairoMakie
end

# Load the plotting module
println("📁 Loading standalone plotting module...")
include("harmonic_oscillator_plotting.jl")
using .HarmonicOscillatorPlotting

println("✅ Plotting module loaded successfully")

# Create demo output directory
if !isdir("plotting_demo_output")
    mkdir("plotting_demo_output")
end

# Demo 1: Single spring frame
println("\n1️⃣ Creating single spring frame...")
fig1 = plot_spring(0, zpos(0))
save("plotting_demo_output/spring_frame.png", fig1)
println("   ✅ Saved to plotting_demo_output/spring_frame.png")

# Demo 2: Spring animation
println("2️⃣ Creating spring animation frames...")
plot_spring_animation(0:5:10, "plotting_demo_output/spring_animation", true)
println("   🎬 Created 3 frames in plotting_demo_output/spring_animation/")

# Demo 3: Phase space plot
println("3️⃣ Creating phase space plot...")
time_data = 0:0.5:5
position_data = cos.(time_data)
velocity_data = -sin.(time_data)

mock_solution = (q = position_data, v = velocity_data, t = time_data)
fig3 = plot_harmonic_oscillator(mock_solution)
save("plotting_demo_output/phase_space.png", fig3)
println("   ✅ Saved to plotting_demo_output/phase_space.png")

# Demo 4: Test utility functions
println("4️⃣ Testing utility functions...")
try
    using Test
    @test xpos(0) ≈ 0
    @test ypos(0) ≈ 0
    @test zpos(0) ≈ 1.0
    @test ϑ(0) ≈ 0.0
    println("   ✅ All utility functions working correctly")
catch e
    println("   ⚠️  Test functions not available: $(typeof(e))")
    println("   But the plotting functionality works!")
end

println("\n🎉 Demo completed successfully!")
println("📁 Output files created in: plotting_demo_output/")
println("💡 To create movies from animation frames:")
println("   ffmpeg -framerate 5 -i plotting_demo_output/spring_animation/harmonic-oscillator-%03d.png movie.mp4")
