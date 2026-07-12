module Diagnostics

    export plot_energy_error, plot_energy_drift, plot_constraint_error, plot_lagrange_multiplier

    # Shared helper used by the plotting extension to build subscripted axis labels.
    subscript(i::Integer) = i < 0 ? error("$i is negative") : join('₀' + d for d in reverse(digits(i)))

    # Generic diagnostic plots. The methods are implemented in the `DiagnosticsPlots`
    # extension, which is loaded together with `Makie`/`CairoMakie`.
    function plot_energy_error end
    function plot_energy_drift end
    function plot_constraint_error end
    function plot_lagrange_multiplier end

end
