# Simple Harmonic Oscillator Plotting Demo
#
# This script demonstrates the plotting functionality without requiring GeometricIntegrators

println("🎨 Simple Harmonic Oscillator Plotting Demo")
println("="^50)

# First check if CairoMakie is available
try
    using CairoMakie
    println("✅ CairoMakie is available")
catch e
    println("❌ CairoMakie not found. Please install with:")
    println("   using Pkg")
    println("   Pkg.add(\"CairoMakie\")")
    exit(1)
end

# Load the plotting module
using GeometricProblems.HarmonicOscillatorPlots

# 1. Demonstrate utility functions (always available)
println("1️⃣ Testing utility functions...")
@show xpos(0)
@show ypos(0)
@show zpos(0)
@show ϑ(0)

# Test spring geometry
spring_geom = spring(10)
println("Spring geometry at frame 10: ", size(spring_geom[1]), " points")

# 2. Create simple animation frames (synthetic data)
println("2️⃣ Creating synthetic animation...")
output_dir = "sh-la-simple-demo"
plot_spring_animation(0:4:8, output_dir, true)  # 3 frames
println("   🎬 Created frames in: $output_dir/")

# 3. Test with simulated harmonic motion data
println("3️⃣ Testing with simulated data...")

# Create mock solution data
TIME_STEPS = 50
time_data = 0:0.1:TIME_STEPS*0.1-0.1
position_data = [cos(0.7 * t) for t in time_data]  # Simulated position
velocity_data = [-0.7 * sin(0.7 * t) for t in time_data]  # Simulated velocity

# Mock solution object
mock_solution = (
    q = position_data,
    v = velocity_data,
    t = time_data
)

# Create phase space plot
println("   📊 Creating phase space plot...")
fig = plot_harmonic_oscillator(mock_solution, ntime = 25)

# Customize the plot
Axis(fig[1,1]; xlabel = "Position", ylabel = "Velocity", title = "Simulated Phase Space")

# Display or save
try
    display(fig)  # Try to display the figure
catch e
    # If display fails, save the figure
    save("simple_phase_space.png", fig)
    println("   🖼️  Saved to: simple_phase_space.png")
end

# 4. Animation from mock solution
println("4️⃣ Creating animation from mock solution...")
plot_harmonic_oscillator_animation(
    mock_solution,
    frames_to_plot = 1:10:length(time_data),
    output_dir = "solution_animation",
    create_dir = true
)
println("   🎥 Created animation frames in: solution_animation/")

println("🎉 Demo completed successfully!")
println("Files created:")
println("   - Phase space plot: simple_phase_space.png")
println("   - Synthetic animation: $output_dir/")
println("   - Solution animation: solution_animation/")
println("\nTo create movies from frames:")
println("   ffmpeg -framerate 10 -i $output_dir/harmonic-oscillator-%03d.png synthetic.mp4")
println("   ffmpeg -framerate 10 -i solution_animation/harmonic-oscillator-sol-%03d.png solution.mp4")
