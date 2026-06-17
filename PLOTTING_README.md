# Harmonic Oscillator Plotting

This directory contains standalone plotting functionality for harmonic oscillators using CairoMakie. The plotting is implemented as a separate module that can be used without modifying the main GeometricProblems.jl package structure.

## 🚀 Quick Start

The plotting functionality is available as a standalone module in `scripts/harmonic_oscillator_plotting.jl`. You can use it in two ways:

### Option 1: Run the Demo Script

```bash
cd scripts
julia plotting_demo.jl
```

This will create a `plotting_demo_output/` directory with example plots and animations.

### Option 2: Use in Your Own Code

```julia
# Move to the repository root
cd /Users/benbradmin/Documents/GeometricProblems

# Include the plotting module
include("scripts/harmonic_oscillator_plotting.jl")
using .HarmonicOscillatorPlotting

# Create a spring animation
plot_spring_animation(0:100, "output_frames", true)

# Create phase space plot with mock data
sol = (q = cos.(0:0.5:10), v = -sin.(0:0.5:10), t = 0:0.5:10)
fig = plot_harmonic_oscillator(sol)
save("phase_space.png", fig)
```

## 📦 Requirements

The plotting functionality requires **CairoMakie**. If it's not installed, the demo script will automatically install it, or you can install it manually:

```julia
using Pkg
Pkg.add("CairoMakie")
```

## 🎨 Available Functions

### Core Plotting Functions

- **`plot_spring(i, q)`** - Plot a single frame of the spring-mass system
  - `i`: Frame index
  - `q`: Position of the mass
  - Returns: A Figure object

- **`plot_spring_animation(frames, output_dir, create_dir)`** - Create spring animation
  - `frames`: Range of frame indices (e.g., `0:100`)
  - `output_dir`: Directory to save PNG files
  - `create_dir`: Whether to create directory if it doesn't exist
  - Saves: Individual PNG frames for each frame

- **`plot_harmonic_oscillator(sol; ntime)`** - Plot harmonic oscillator in phase space
  - `sol`: A structure with position, velocity, and time fields
  - `ntime`: Number of time steps to plot (default: all)
  - Returns: A Figure object showing trajectory in phase space

### Utility Functions

- **`xpos(i)`, `ypos(i)`, `zpos(i)`** - Calculate position coordinates for animation
- **`ϑ(i)`** - Calculate velocity for animation
- **`spring(i)`** - Generate spring geometry for a given frame

## 📝 Example Usage

### Example 1: Basic Spring Animation

```julia
include("scripts/harmonic_oscillator_plotting.jl")
using .HarmonicOscillatorPlotting

# Create animation with 10 frames
plot_spring_animation(0:9, "my_animation", true)
```

Then create a movie:
```bash
ffmpeg -framerate 10 -i my_animation/harmonic-oscillator-%03d.png movie.mp4
```

### Example 2: Phase Space Plot

```julia
include("scripts/harmonic_oscillator_plotting.jl")
using .HarmonicOscillatorPlotting

# Create mock data representing harmonic motion
time_data = 0:0.1:10
position_data = cos.(time_data)
velocity_data = -sin.(time_data)

sol = (q = position_data, v = velocity_data, t = time_data)

# Create phase space plot
fig = plot_harmonic_oscillator(sol)
save("phase_space.png", fig)
```

### Example 3: Custom Spring Frame

```julia
include("scripts/harmonic_oscillator_plotting.jl")
using .HarmonicOscillatorPlotting

# Create a single spring frame at maximum displacement
fig = plot_spring(0, 1.0)
save("max_displacement.png", fig)

# Create spring frame at equilibrium
fig = plot_spring(10, 0.0)
save("equilibrium.png", fig)
```

## 🔧 Technical Details

### Physical Constants

The plotting uses the same physical constants as the original script:
- Spring constant `k = 0.5`
- Mass `m = 1.0`
- Angular frequency `ω = √(k/m)`
- Amplitude `A = 1.0`

### Visualization Parameters

- `n = 100` - Number of coils in the spring
- `w = 10` - Winding density
- `r = 0.1` - Radius of spring
- `rod = 0.15` - Length of rigid rods at spring ends

### File Format

- **Spring Animations**: PNG files named `harmonic-oscillator-000.png`, `harmonic-oscillator-001.png`, etc.
- **Phase Space Plots**: PNG files with customizable filenames
- **Single Frames**: PNG files that can be saved with any filename

## 📁 File Structure

```
GeometricProblems/
├── scripts/
│   ├── harmonic_oscillator_plotting.jl  # Main plotting module
│   └── plotting_demo.jl                  # Demo script
├── examples/
│   └── standalone_plotting.jl            # Alternative standalone demo
└── PLOTTING_README.md                    # This file
```

## 🎯 Why This Approach?

The plotting functionality is implemented as a standalone module rather than integrated into the main package for several reasons:

1. **No Package Structure Changes**: Doesn't modify the GeometricProblems.jl main package
2. **No Dependency Conflicts**: Avoids complex dependency resolution with the main package
3. **Flexible Usage**: Can be used independently or integrated as needed
4. **Backward Compatible**: Doesn't break any existing functionality
5. **Easy Maintenance**: Plotting code is separate and easier to maintain

## 🚀 Future Integration

While this standalone approach works well, the plotting functions could potentially be integrated into the main package structure if desired. The current implementation uses:

- Pure Julia functionality without complex module hierarchies
- Direct CairoMakie integration
- Flexible data handling for different solution formats

This makes it straightforward to eventually integrate into the main package if needed, once dependency management issues are resolved.

## 💡 Tips

1. **Frame Rate**: When creating animations, adjust the frame range and ffmpeg framerate for smooth motion
2. **Memory**: Large animations may consume memory; consider processing in batches
3. **Customization**: Modify the visualization parameters in the source file for different spring sizes
4. **Integration**: Works with both synthetic data and actual GeometricIntegrators.jl solutions

## 📧 Questions?

If you encounter any issues or have questions about the plotting functionality, please refer to the demo scripts or check the source code in `scripts/harmonic_oscillator_plotting.jl`.