@doc raw"""
# Harmonic Oscillator Visualization

This module provides visualization routines for the harmonic oscillator using CairoMakie.
The visualizations include spring-mass animations and phase-space plots.

Note: To use these plotting routines, you need to have CairoMakie installed:
```julia
using Pkg
Pkg.add("CairoMakie")
```
"""
module HarmonicOscillatorPlots

    using Requires: @require

    """
    Create a harmonic oscillator ODE problem with visualization-friendly parameters.

    This function creates a problem suitable for visualization, with parameters
    that match the animation system in the plotting routines.

    Keyword arguments:
    - `q₀::Vector=[1.0, 0.0]`: Initial condition [position, velocity]
    - `timespan::Tuple=(0.0, 10.0)`: Time span for the simulation
    - `timestep::Real=0.01`: Time step for integration

    Returns:
    - A GeometricIntegrators.ODEProblem configured for visualization
    """
    function visualization_problem(cfg::NamedTuple)
        # Import here to avoid circular dependencies
        using ..HarmonicOscillator: odeproblem

        # Extract configuration with defaults
        cfg = merge((q₀=[1.0, 0.0], timespan=(0.0, 10.0), timestep=0.01), cfg)

        # Create problem using the standard interface
        # Note: The harmonic oscillator module uses k=0.5, m=1.0 by default
        odeproblem(cfg.q₀; timespan=cfg.timespan, timestep=cfg.timestep)
    end

end  # module HarmonicOscillatorPlots

# Requires.jl functionality for optional CairoMakie dependency
function __init__()
    @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf8484" begin
        # Import plotting functions from implementation file
        include("harmonic_oscillator_plots_implementation.jl")

        # Export all the plotting functionality
        for sym in names(@__MODULE__(), all=true)
            if startswith(string(sym), "plot_") || sym in [:spring, :xpos, :ypos, :zpos, :ϑ]
                @eval export $sym
            end
        end
    end
end
