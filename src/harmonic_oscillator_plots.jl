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

    using GeometricIntegrators
    using ..HarmonicOscillator: harmonic_oscillator_ode, default_parameters
    using Requires

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
        # Extract configuration with defaults
        cfg = merge((q₀=[1.0, 0.0], timespan=(0.0, 10.0), timestep=0.01), cfg)

        # Use parameters that match animation parameters from the plot module
        # k = 0.5, m = 1.0 from the original animation script
        params = (m=1.0, k=0.5, ω=sqrt(0.5))

        harmonic_oscillator_ode(cfg.q₀; cfg.timespan, timestep=cfg.timestep, parameters=params)
    end

end  # module HarmonicOscillatorPlots

# Requires.jl functionality for optional CairoMakie dependency
function __init__()
    @require CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf8484" begin
        # Include all the plotting functions that depend on CairoMakie
        include("harmonic_oscillator_plots_implementation.jl")

        # Export the functions from the implementation
        for sym in names(@__MODULE__(), all=true)
            if startswith(string(sym), "plot_") || sym in [:spring, :xpos, :ypos, :zpos, :ϑ]
                @eval export $sym
            end
        end
    end
end
