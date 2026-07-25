module Diagnostics

    export plot_energy_error, plot_energy_drift, plot_constraint_error, plot_lagrange_multiplier
    export plot_invariant_error, plot_invariant_drift
    export plot_convergence, plot_order

    # Shared helper used by the plotting extension to build subscripted axis labels.
    subscript(i::Integer) = i < 0 ? error("$i is negative") : join('₀' + d for d in reverse(digits(i)))

    # Generic diagnostic plots. The methods are implemented in the `DiagnosticsPlots`
    # extension, which is loaded together with `Makie`/`CairoMakie`.

    # `plot_invariant_error` and `plot_invariant_drift` work for any invariant carried by a
    # problem, addressed by its key in `invariants(problem)`. The energy variants are the
    # `:h` special case, kept separate because it is by far the most common one.
    function plot_invariant_error end
    function plot_invariant_drift end
    function plot_energy_error end
    function plot_energy_drift end
    function plot_constraint_error end
    function plot_lagrange_multiplier end

    # Convergence diagnostics over a sequence of time steps.
    function plot_convergence end
    function plot_order end

end
