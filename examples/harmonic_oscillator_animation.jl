# Harmonic Oscillator Animation Example
#
# This script demonstrates how to use the new harmonic oscillator plotting functionality
# to create animations of the spring-mass system.

println("Harmonic Oscillator Animation Example")
println("="^50)

# First, make sure we have the required packages
try
    using GeometricProblems
    using GeometricIntegrators
    using CairoMakie
catch e
    println("Installing required packages...")
    import Pkg
    Pkg.add([
        "GeometricProblems",
        "GeometricIntegrators",
        "CairoMakie"
    ])
    using GeometricProblems
    using GeometricIntegrators
    using CairoMakie
end

# Example 1: Create a visualization-friendly problem
println("Creating visualization problem...")
using GeometricProblems.HarmonicOscillatorPlots: visualization_problem

# Create a problem suitable for animation
prob = visualization_problem((
    q₀ = [1.0, 0.0],           # Start from rest at amplitude 1.0
    timespan = (0.0, 5.0),     # Simulate for 5 seconds
    timestep = 0.01            # Small timestep for smooth animation
))

println("Problem configuration:")
println("  Initial position: $(prob.q₀)")
println("  Time span: $(prob.timespan)")
println("  Time step: $(prob.Δt)")

# Example 2: Solve the problem
println("\nSolving harmonic oscillator...")
sol = integrate(prob, Gauss(2))

println("Solution info:")
println("  Number of time steps: $(length(sol.t))")
println("  Final time: $(sol.t[end])")
println("  Final position: $(sol.q[end])")

# Example 3: Basic phase space plot
println("\nCreating phase space plot...")
fig_phase = plot_harmonic_oscillator(sol)
save("examples/harmonic_oscillator_phase_space.png", fig_phase)
println("Saved phase space plot to: examples/harmonic_oscillator_phase_space.png")

# Example 4: Create animation frames
println("\nCreating animation frames...")
# Use the animation function that works with the solution
plot_harmonic_oscillator_animation(
    sol,
    frames_to_plot = 1:20:length(sol.t),  # Every 20th frame
    output_dir = "examples/harmonic_oscillator_frames",
    create_dir = true
)
println("Created $(length(1:20:length(sol.t))) animation frames in: examples/harmonic_oscillator_frames/")

# Example 5: Create synthetic animation (like original script)
println("\nCreating synthetic animation...")
plot_spring_animation(
    0:100,                          # 101 frames
    output_dir = "examples/synthetic_animation",
    create_dir = true
)
println("Created synthetic animation frames in: examples/synthetic_animation/")

println("\n" * 3)
println("✅ Animation examples completed!")
println("📁 Images and frames saved to the 'examples/' directory.")
println("⚠️  To create actual movies, you'll need to use a tool like 'ffmpeg'")

# Optional: Show how to create a movie using ffmpeg (uncomment if you have ffmpeg installed)
#println("\nCreating movie with ffmpeg...")
#run(`ffmpeg -framerate 30 -i examples/synthetic_animation/harmonic-oscillator-%03d.png \
#          -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p \
#          examples/harmonic_oscillator.mp4`)
