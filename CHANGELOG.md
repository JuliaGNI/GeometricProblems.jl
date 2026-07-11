# Changelog

All notable changes from the correctness audit (see `bugs.md`) and its remediation plan
(`plan.md`). Entries are grouped by type. Bug IDs (B*/C*) and plan IDs (P*) reference those files.

Categories: **Bug fixes** = code defects (typos, wrong API calls, crashes, bad indices);
**Model fixes** = corrections to a model's mathematical/physical formulation; **New features**;
**Documentation**; **Tests**; **Repository hygiene**.

## [Unreleased]

### Bug fixes
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
- **Doc pages** (P14): filled the empty module docstrings for the Lorenz attractor and
  Lotka-Volterra 4D model, expanded the Point Vortices docstring, and added content to the
  previously-empty pages for the Kepler problem, Hénon-Heiles system, inner/outer solar system, and
  diagnostics (marked *not yet implemented* where no backing model exists).

### Tests
- Added tests (P15): nonlinear oscillators and three-body problem (HODE↔LODE trajectory
  agreement), massless charged particle (`Ω·v = −∇H`, i.e. ODE ↔ variational sign consistency —
  the direct regression for B5), and linear wave (constructor regression for B7).
  `test/{nonlinear_oscillators,three_body,massless_charged_particle,linear_wave}_tests.jl`.
- Re-calibrated the Lotka-Volterra 3D test tolerances to the corrected dynamics (P8), and refreshed
  the stale `reference_solution`. `test/lotka_volterra_3d_tests.jl`.
- The Kubo oscillator test remains **disabled**: re-enabling it surfaced a pre-existing
  GeometricEquations SDE-API incompatibility (old `SDEProblem(m, n, v, B, …)` vs the new
  `noise`-object API). The Kubo model itself is correct; the SDE-API migration is tracked as a
  follow-up. `test/runtests.jl`.

### Repository hygiene
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
