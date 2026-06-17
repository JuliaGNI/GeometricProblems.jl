@doc raw"""
# Harmonic Oscillator Visualization

This module provides visualization routines for the harmonic oscillator using CairoMakie.
The visualizations include spring-mass animations and phase-space plots.

## Note

This module provides basic utility functions that don't require CairoMakie.
For full plotting functionality, use the standalone plotting scripts in `scripts/`:
```julia
include("scripts/harmonic_oscillator_plotting.jl")
using .HarmonicOscillatorPlotting
```

## Available Functions

- `xpos(i)`, `ypos(i)`, `zpos(i)`, `ϑ(i)` - Utility coordinate functions (no CairoMakie required)

## For Full Plotting Functionality

See `scripts/harmonic_oscillator_plotting.jl` for complete plotting capabilities
that work with or without CairoMakie as a weak dependency.
"""
module HarmonicOscillatorPlots

# Physical constants (always available, no CairoMakie needed)
const k = 0.5        # Spring constant
const m = 1.0        # Mass
const ω = sqrt(k/m)  # Angular frequency

# Visualization parameters
const n = 100        # Number of coils in the spring
const w = 10         # Winding density
const r = 0.1        # Radius of spring
const h = 2π / ω / n  # Time step for visualization
const A = 1.0        # Amplitude

# Utility functions (always available)
xpos(i) = sin(r*h*i)
ypos(i) = sin(r*h*i)
zpos(i) = A*cos(ω*h*i)
ϑ(i) = -ω*A*sin(ω*h*i)

# Spring dimensions
const rod = 0.15        # Length of rigid rods at spring ends
const top = zpos(0) + 1  # Top position

export xpos, ypos, zpos, ϑ

end  # module HarmonicOscillatorPlots
