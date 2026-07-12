# Diagnostics

Several models provide diagnostic helpers to quantify how well a numerical integrator preserves the
geometric structure of the problem. Depending on the model these include:

* `compute_energy_error` — the (relative) drift of the Hamiltonian / energy along a solution;
* `compute_invariant_error` — the drift of a general invariant (from `GeometricSolutions`);
* `compute_casimir_error` — the drift of a Casimir of the Poisson structure (e.g. for the
  [Lotka-Volterra 3d](lotka_volterra_3d.md) model);
* `compute_momentum_error` — the drift of a momentum map / one-form.

These functions are defined per model (see the respective module documentation) and typically take
the time series, the solution, and the problem parameters, returning both the invariant time series
and its relative error.
