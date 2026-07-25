# Changelog

All notable changes from the correctness audit (see `bugs.md`) and its remediation plan
(`plan.md`). Entries are grouped by type. Bug IDs (B*/C*) and plan IDs (P*) reference those files.

Categories: **Bug fixes** = code defects (typos, wrong API calls, crashes, bad indices);
**Model fixes** = corrections to a model's mathematical/physical formulation; **New features**;
**Documentation**; **Tests**; **Repository hygiene**.

## [Unreleased]

### Bug fixes
- **Massless charged particle**: added `hamiltonian(t, q, p, params)`, which both massless modules
  were missing. All their implicit problems declare `invariants=(h=hamiltonian,)`, but
  GeometricEquations evaluates that invariant with an additional momentum/velocity slot, so
  e.g. `compute_invariant_error(t, q, p, params, invariants(iodeproblem())[:h])` raised
  `MethodError: no method matching hamiltonian(::Float64, ::Vector, ::Vector, ::NamedTuple)`.
  Mirrors the existing shim in Lotka-Volterra 2D/4D and the harmonic oscillator.
  `src/massless_charged_particle_common.jl`.
- **Lotka-Volterra 2D** (B10/P10): `lotka_volterra_2d_v_dae` now reads the velocities from
  `q[3]`/`q[4]` instead of its own uninitialised output `v[3]`/`v[4]`.
  `src/lotka_volterra_2d_common.jl`.
- **Lotka-Volterra 3D/4D plots** (B11/P11): the plot recipes' `compute_energy_error` now resolves
  `hamiltonian` (imported from the sibling model module) instead of raising `UndefVarError`.
  `src/lotka_volterra_3d_plots.jl`, `src/lotka_volterra_4d_plots.jl`.
- **Bump initial condition** (B12/P12): the vector method of `∂s` now calls `∂s` (was `s`), so it
  returns derivatives rather than values. `src/bump_initial_condition.jl`.
- **Lorenz attractor** (B13/P12): removed the export of the never-defined `plot_lorenz_attractor`.
  `src/lorenz_attractor.jl`.
- **Harmonic oscillator** (B1/P6): fixed `hamitlonian` typo in `oscillator_pdae_ϕ` and made the
  constraint energy-relative (`H − H₀`), so the `pdae`/`hdae`/`idae`/`ldae` constraints no longer
  raise `UndefVarError`. `src/harmonic_oscillator.jl`.
- **Three-body problem** (B4/P3): repaired `lodeproblem`, which referenced undefined
  `θ₀`/`p₀`/`θ̇` and used `lagrangian_variables(2)`; now uses the module ICs and
  `lagrangian_variables(6)`. `src/three_body_problem.jl`.
- **Massless charged particle** (B6/P4): fixed the argument order in the (dead) `ψ` constraint
  (`g₁(v,t,q,…)`). `src/massless_charged_particle.jl`.
- **Linear wave** (B7/P5): `lodeproblem` now calls `LODEProblem` (was a lowercase self-call →
  `MethodError`). `src/linear_wave.jl`.

### Model fixes
- **Coupled harmonic oscillator** (B2/P1): corrected the Lagrangian kinetic term
  `q̇²/(2m)` → `m·q̇²/2`, so the `LODE` and `HODE` describe the same system (previously differed by
  a factor `m²`). `src/coupled_harmonic_oscillator.jl`.
- **Three-body problem** (B3/P2): added the missing body-1↔body-3 gravitational interaction to the
  potential `V` (was only pairs (1,2) and (2,3)). `src/three_body_problem.jl`.
- **Three-body problem** (Copilot review): fixed the Lagrangian kinetic term — it reused the
  momentum-form `T(q̇)=q̇²/(2m)` instead of `m·q̇²/2`, so the LODE and HODE disagreed for non-unit
  masses (the same class of bug as the coupled oscillator); also added a `v̄` momentum→velocity map
  to the `lodeproblem`. `src/three_body_problem.jl`.
- **Massless charged particle** (B5/P4): flipped the sign-reversed `odeproblem` vector field (and
  the documented Hamiltonian-form matrix) so the ODE agrees with the variational IODE/IDAE forms;
  the one-form `ϑ = A` is retained as authoritative. `src/massless_charged_particle.jl`.
- **Harmonic oscillator** (C1/P7): threaded the mass `m` consistently through every equation form
  (ODE, PODE/HODE, SODE, IODE/LODE, degenerate IODE/LODE, DAE/PDAE constraint terms) and the
  exact/reference solution. Previously `m` was silently fixed to 1 in the dynamics while the
  Hamiltonian, Lagrangian, energy invariant and DELE derivatives already carried it, so the forms
  agreed only at `m = 1`. `src/harmonic_oscillator.jl`.
- **Lotka-Volterra 3D** (C2/P8): adopted the documented Hamiltonian `H = a·q + b·ln q` (the code
  used `− b·ln q`); flipped the `b`-term signs in the vector field `v₁,v₂,v₃` accordingly, so the
  integrated equations match the documentation (verified `v = P∇H`). `src/lotka_volterra_3d.jl`.
- **Lotka-Volterra 2D (gauge)** (B9/P9): implemented the documented gauge one-form (added the
  `d(q₁q₂)` terms `+q₂`/`+q₁` to `ϑ` and the corresponding `+1` in `dϑ₁dx₂`/`dϑ₂dx₁`); the module
  is now a genuinely distinct, gauge-equivalent variational integrator rather than a copy of the
  symmetric variant. `src/lotka_volterra_2d_gauge.jl`.

### New features
- **Generic invariant diagnostics**: `Diagnostics.plot_invariant_error` and
  `Diagnostics.plot_invariant_drift` plot the relative error (and its per-interval drift) of *any*
  invariant a problem carries, selected by its key in `invariants(problem)` (`:h` by default) or
  by passing the invariant function itself. `plot_energy_error`/`plot_energy_drift` are now the
  `:h` special case of these and keep their previous behaviour and axis labels.
  Motivated by problems with more than one conserved quantity — the point vortices below, and the
  guiding-centre toroidal momentum in `ChargedParticleDynamics` — which previously each needed
  their own bespoke plotting code. `src/diagnostics.jl`, `ext/DiagnosticsPlots.jl`.
- **Convergence diagnostics**: `Diagnostics.plot_convergence(h, ε; order)` (log-log error over
  time step, with an optional reference slope of the expected order) and
  `Diagnostics.plot_order(h, p)` (observed order over time step). `src/diagnostics.jl`,
  `ext/DiagnosticsPlots.jl`.
- **Point vortices**: added a `PointVorticesPlots` Makie extension providing
  `plot_phase_portrait`, `plot_solution` and `plot_traces` for both `PointVortices` and
  `PointVorticesLinear`, mirroring `MasslessChargedParticlePlots`. Both modules' problems now also
  declare their invariants, `invariants=(h=hamiltonian, p=angular_momentum)`, which they were
  missing entirely, so that the generic diagnostics above work on them; `angular_momentum` gained
  the four-argument `(t, q, p, params)` method the implicit formulations evaluate it with.
  `ext/PointVorticesPlots.jl`, `src/point_vortices.jl`, `src/point_vortices_linear.jl`.
- **Massless charged particle**: added `lodeproblem` and `ldaeproblem` to both
  `MasslessChargedParticle` and `MasslessChargedParticleSingular`, completing the variational
  formulations. This required the two previously missing ingredients, both derived from what the
  modules already provide: the degenerate Lagrangian `lagrangian(t, q, v, params) = ϑ⋅v − H` and
  the symplectic two-form `ω`, assembled from the per-gauge derivatives `dϑᵢdxⱼ` with the
  convention `Ωᵢⱼ = ∂ϑᵢ/∂qⱼ − ∂ϑⱼ/∂qᵢ` shared with Lotka-Volterra 2D/4D and the point vortices
  (`Ω = [0 −B; +B 0]` for both gauges, so `Ω v = −∇H`). `ω` is exported in an allocating form
  `ω(t, q, params)` mirroring `ϑ`, alongside the in-place callback
  `massless_charged_particle_ω`. Also added the `ū`/`ḡ` projection callbacks and the eight-argument
  `ψ` method that `LDAEProblem` requires, and renamed the existing `ψ` arguments to
  `(ψ, t, q, p, q̇, ṗ, params)` to match the interface it implements.
  `src/massless_charged_particle_common.jl`, `src/massless_charged_particle.jl`,
  `src/massless_charged_particle_singular.jl`.
- **Plotting → CairoMakie extensions** (⚠️ breaking): migrated the remaining package plot recipes
  from Plots.jl/`RecipesBase` to Makie extensions, matching the pattern already used for the
  harmonic oscillator and pendulum. There is now one extension per problem plus a shared diagnostics
  extension, each activated when a Makie backend (e.g. `CairoMakie`) is loaded:
  - `DiagnosticsPlots` — generic diagnostics `plot_energy_error`, `plot_energy_drift`,
    `plot_constraint_error`, `plot_lagrange_multiplier` (stubs live in the new `Diagnostics` module,
    `src/diagnostics.jl`, replacing `src/plot_recipes.jl`).
  - `LotkaVolterra2dPlots` — `plot_solution`, `plot_phase_portrait`, `plot_traces`.
  - `LotkaVolterra3dPlots` / `LotkaVolterra4dPlots` — `plot_phase_portrait`, `plot_traces`.
  - `MasslessChargedParticlePlots` — `plot_solution`, `plot_phase_portrait`, `plot_traces`.

  The old `@userplot` names (`PlotEnergyError`, `Plot_Lotka_Volterra_2d`, …) are replaced by the
  snake_case functions above, which now build and return a Makie `Figure` instead of a Plots object.
  The massless-particle module's former energy/momentum-error recipes are dropped in favour of the
  generic `Diagnostics` ones. `ext/{DiagnosticsPlots,LotkaVolterra2dPlots,LotkaVolterra3dPlots,LotkaVolterra4dPlots,MasslessChargedParticlePlots}.jl`.
- **Harmonic-oscillator / pendulum plotting helpers**: extended the `HarmonicOscillatorPlots` and
  `PendulumPlots` extensions with `plot_phase_portrait(sol)`, `plot_traces(sol)`, and
  `plot_hamiltonian(; …)` (energy-landscape contour), alongside the existing animation-frame
  `plot_solution`. `plot_phase_portrait(sol)` is now a uniform single-argument function across all
  problem extensions. `ext/{HarmonicOscillatorPlots,PendulumPlots}.jl`.
- **Documentation figures now use the extensions**: the phase-portrait / trajectory figures on the
  Lotka-Volterra 2D/3D/4D, massless-charged-particle, harmonic-oscillator and pendulum pages, the
  two Hamiltonian energy-landscape figures (harmonic oscillator, pendulum), and the pendulum
  time-trace figure now call the extension plotting functions instead of hand-rolling CairoMakie.
  `diagnostics.md` gained worked examples of `Diagnostics.plot_energy_error` and
  `plot_constraint_error`. (Ensemble/parameter-sweep overlays and the extension-less models' one-line
  trajectory plots are left as-is.)
- **Nonlinear oscillators** (B15/P13): implemented the previously empty stub modules
  **Duffing**, **Lennard-Jones**, **Morse**, and **Mathews-Lakshmanan** oscillators. Each provides
  `hamiltonian`, `lagrangian`, and EulerLagrange-generated `hodeproblem`/`lodeproblem` constructors,
  with documentation of the standard textbook forms.
  `src/{duffing,lennard_jones,morse,mathews_lakshmanan}_oscillator.jl`.

### Documentation
- **Coupled harmonic oscillator** (B2/P1): fixed docstring Hamiltonian kinetic term `q` → `p`.
- **Massless charged particle** (B5/P4): corrected the Hamiltonian-form matrix in the docstring.
- **Three-body problem**: documented the previously omitted `m₃` and `G` parameters.
- **Lotka-Volterra 4D (Lagrangian)** (C4/P12): removed the stray `b₅ log q⁵` term from the
  docstring Hamiltonian (the system has four components). `src/lotka_volterra_4d_lagrangian.jl`.
- **Linear wave** (C3/P12): reconciled the documented discrete Hamiltonian `H_h` with the
  implemented one (unweighted kinetic term, `μ²/(4Δξ²)` potential), clarified that the two boundary
  points are dynamical coordinates, and fixed LaTeX/grammar typos (`\Delta_xi`, "spaces points",
  "an completely-integrable"). `docs/src/linear_wave.md`.
- **harmonic_oscillator.md**: wrote the previously empty "Lagrangian Formulation" and "Dynamics"
  sections (the latter with a CairoMakie phase-space simulation of the default problem).
- **rigid_body.md**: explained that all trajectories lie on the Casimir sphere $\lVert x \rVert^2$,
  so the example's `[sin θ, 0, cos θ]` initial conditions (vs. the module's `[cos θ, 0, sin θ]`)
  merely pick different orientations on that sphere — not a discrepancy.
- **Default-case simulation plots**: added CairoMakie simulations of the default problem to the
  docs of the double pendulum, Lorenz attractor, Lotka-Volterra 3D/4D, point vortices, Toda
  lattice, and the four nonlinear oscillators.
- **New navigation pages**: added documentation pages (with docstrings and a CairoMakie plot) for
  the Kubo oscillator and the linear point-vortices model, and registered them in `docs/make.jl`.
- **Doc pages** (P14): filled the empty module docstrings for the Lorenz attractor and
  Lotka-Volterra 4D model, expanded the Point Vortices docstring, and added content to the
  previously-empty pages for the Kepler problem, Hénon-Heiles system, inner/outer solar system, and
  diagnostics (marked *not yet implemented* where no backing model exists).

### Tests
- **Massless charged particle two-form and Lagrangian forms**: added testsets to
  `test/massless_charged_particle_tests.jl` and `test/massless_charged_particle_singular_tests.jl`
  that compare `ω` against `ForwardDiff.jacobian` of `ϑ` (validating both the assembly and the
  hand-coded `dϑᵢdxⱼ`), check antisymmetry, `Ω₁₂ = −B`, `Ω v = −∇H`, gauge invariance of `Ω` across
  the two modules, and `∂L/∂v = ϑ`; and that integrate `lodeproblem` with
  `MidpointProjection(VPRKGauss(2))` and `ldaeproblem` with `SLRKLobattoIIIE(2)` against the
  `Gauss(8)` ODE reference. Since no integrator in GeometricIntegrators evaluates `ω` and
  `check_methods` skips it, the tests also assert the positional wiring via `functions(prob).ω` /
  `.l` — an `ω`/`l` swap would otherwise be silent (as it is in Lotka-Volterra 2D).
- **Plotting extensions**: added `test/plots_tests.jl` (wired into `test/runtests.jl`) — smoke
  tests that each migrated plot function builds a Makie `Figure` for a short integrated solution,
  covering the `Diagnostics`, Lotka-Volterra 2D/3D/4D, and massless-charged-particle plots. Added
  `CairoMakie` to `test/Project.toml`. Plots were previously untested.
- Added tests (P15): nonlinear oscillators and three-body problem (HODE↔LODE trajectory
  agreement), massless charged particle (`Ω·v = −∇H`, i.e. ODE ↔ variational sign consistency —
  the direct regression for B5), and linear wave (constructor regression for B7).
  `test/{nonlinear_oscillators,three_body,massless_charged_particle,linear_wave}_tests.jl`.
- Re-calibrated the Lotka-Volterra 3D test tolerances to the corrected dynamics (P8), and refreshed
  the stale `reference_solution`. `test/lotka_volterra_3d_tests.jl`.
- **Kubo oscillator**: migrated the SDE/PSDE/SPSDE constructors to the GeometricEquations
  noise-object API (added a `KuboNoise <: AbstractStochasticProcess` marker; the old
  `SDEProblem(m, n, v, B, …)` positional form is gone). The single-IC variants build problems and
  the multi-IC variants (`_3`) build ensembles. The Kubo test is **re-enabled** and rewritten to
  the new API (it now checks the drift/diffusion terms directly).
  `src/kubo_oscillator.jl`, `test/kubo_oscillator_tests.jl`, `test/runtests.jl`.
- Added a three-body HODE↔LODE consistency check with non-unit masses (guards the Lagrangian
  kinetic-term fix). `test/three_body_tests.jl`.

### Repository hygiene
- **Dropped plotting-only dependencies**: with the recipes gone, `RecipesBase`, `Measures`, and
  `LaTeXStrings` were removed from `[deps]`/`[compat]` (they were used only by the deleted plot
  files). Deleted `src/plot_recipes.jl` and `src/{lotka_volterra_2d,lotka_volterra_3d,lotka_volterra_4d,massless_charged_particle}_plots.jl`,
  and removed their `include`s from `src/GeometricProblems.jl`. Registered the five new extensions in
  `Project.toml`. `docs/src/lotka_volterra_2d.md`'s stale `@autodocs` on the removed
  `LotkaVolterra2dPlots` module now resolves through the problem module (docstrings live on the
  extended functions).
- **Docs plotting → CairoMakie**: migrated all documentation plotting from Plots.jl and GLMakie.jl
  to CairoMakie.jl. Rewrote the `@example`/`@eval` blocks in `abc_flow.md`, `initial_condition.md`,
  `coupled_harmonic_oscillator.md`, `massless_charged_particle.md`, `lotka_volterra_2d.md`, and
  `toda_lattice.md`, and dropped `Plots` and `GLMakie` from `docs/Project.toml`. The package's own
  plotting recipes (RecipesBase, unchanged) are unaffected — Plots.jl was only ever a docs
  dependency.
- **`.gitignore`** (P17): added `*-prev.jl` backups and the `references/` and
  `harmonic-oscillator-plots/` scratch directories. No files were deleted.
- **Lotka-Volterra 3D (Lagrangian)** (B14/P16): quarantined the unregistered, broken
  `src/lotka_volterra_3d_lagrangian.jl` — renamed its module `LotkaVolterra3d` →
  `LotkaVolterra3dLagrangian` to remove the name collision with the registered module, and added a
  warning docstring. It remains unregistered (a 3-D Lotka-Volterra admits no non-degenerate
  Lagrangian); not deleted.
