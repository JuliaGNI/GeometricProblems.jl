# Changelog

All notable changes to GeometricProblems.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Categories: **Bug fixes** = code defects (typos, wrong API calls, crashes, bad indices);
**Model fixes** = corrections to a model's mathematical/physical formulation; **Added**;
**Changed**; **Documentation**; **Tests**; **Repository hygiene**.

> Development notes for the 0.7.0 correctness audit — the original findings report, its
> remediation plan and the execution log — are archived under [`docs/dev/`](docs/dev/).

## [Unreleased]

### Model fixes
- **Symplectic two-form `ω` of a regular Lagrangian**: `ω` is now the full ``2n × 2n`` Lagrange
  two-form on `(q, q̇)`, not the ``n × n`` `q`–`q` block. A regular Lagrangian is a second-order
  system of `n` equations, equivalently first order in `2n`; a degenerate one is already first order
  in `n` and keeps its ``n × n`` form. This makes the convention consistent within this package and
  with EulerLagrange, which had generated an identically-zero `ω` (it antisymmetrised `∂L/∂z`, giving
  `d(dL) ≡ 0`) until JuliaGNI/EulerLagrange.jl#24.

  **This supersedes the `ω!` entry under Bug fixes below**, which recorded `Ω = 0` for the regular
  harmonic oscillator. That was right for the `q`–`q` block but is not the two-form the LODE wants:
  `HarmonicOscillator.ω!` now returns the ``2 × 2`` canonical `[0 -m; m 0]`. `degenerate_ω!` is
  unchanged. `src/harmonic_oscillator.jl`.

### Added
- **Outer solar system**: hand-written vector fields. `hodeproblem`/`lodeproblem` no longer generate
  the equations of motion symbolically by default; they use in-place `v`, `f`, `ϑ`, `g` and `ω!`
  built on a shared `∇V!` kernel that computes each of the `N(N-1)/2 = 15` pairwise distances once
  per force evaluation. The generated force was 2.7× slower. `hamiltonian_system` and
  `lagrangian_system` remain, and `symbolic = true` routes the problems through them for
  cross-checking; the two agree to one ulp. `src/outer_solar_system.jl`.
- `benchmark/outer_solar_system.jl` and `benchmark/simplify_evaluation.jl`, backing the two changes
  above with construction time, generated-code size and per-call evaluation cost.
- **All symbolically generated problems now pass `nanmath = true`** to
  `LagrangianSystem`/`HamiltonianSystem`/`DegenerateLagrangianSystem` (25 call sites across 13
  modules). The generated code uses the `NaNMath` variants of `log`, `sqrt`, `^` and friends, so an
  inadmissible state returns `NaN` instead of throwing a `DomainError`. Implicit integrators solve a
  nonlinear system per step and can probe states outside a model's domain — a negative population in
  Lotka-Volterra, a non-positive radius in Lennard-Jones or Morse — where an exception aborts the run
  while a `NaN` lets the solver reject the trial step and continue. The keyword is new in
  EulerLagrange 0.5; before that the behaviour was hardcoded off.

### Bug fixes
- **Lotka-Volterra 2D**: `ldaeproblem` and `ldaeproblem_slrk` passed the symplectic matrix `ω` and
  the Lagrangian `l` to `LDAEProblem` in the wrong positional order, so `functions(prob).ω` *was*
  the Lagrangian and vice versa. 0.7.2 fixed this in `lodeproblem` but not in the two `LDAE`
  constructors. `src/lotka_volterra_2d_equations.jl`.
- **Harmonic oscillator**: `ω!` had no `(Ω, t, q, v, params)` method, so evaluating the symplectic
  matrix of `lodeproblem`, `degenerate_lodeproblem` or `ldaeproblem` raised a `MethodError`; it also
  wrote a 2×2 matrix into the 1×1 one the single-degree-of-freedom forms allocate, and its
  hard-coded value `[0 -1; +1 0]` was neither the mass-aware nor the sign-correct two-form. It is
  now split in two: `ω!` for the regular Lagrangian (`Ω = 0`, since `ϑ = m v` does not depend on
  `q`, and size-agnostic) and `degenerate_ω!` for the degenerate formulation
  (`Ω = [0 m; -m 0]`, so that `Ω v = -∇H` reproduces the degenerate velocity field). Both gained
  the five-argument method. `degenerate_lodeproblem` now also passes `degenerate_lagrangian`
  instead of the regular `lagrangian`, which was dead code. `src/harmonic_oscillator.jl`.
- **Lotka-Volterra 4D**: added the missing `(Ω, t, q, v, params)` method of
  `lotka_volterra_4d_ω`, which `lodeproblem`/`ldaeproblem`/`ldaeproblem_secondary` evaluate.
  `src/lotka_volterra_4d.jl`.
- **Harmonic oscillator**: the energy constraint of the differential-algebraic forms
  (`daeproblem`, `pdaeproblem`, `hdaeproblem`, `idaeproblem`, `ldaeproblem`) referenced the
  *module-level default* initial condition, so `ϕ = H - H₀` did not vanish at `t₀` for any other
  initial data. `H₀` is now captured per problem from the problem's own initial condition.
  `src/harmonic_oscillator.jl`.
- **Harmonic oscillator**: the default initial momentum of `degenerate_iodeproblem` /
  `degenerate_lodeproblem` was built from a mass-free `ϑ`, while the problems' own `ϑ` returns
  `m·q₂` since 0.7.0. The two disagreed by a factor `m`, masked only because the default initial
  condition has zero velocity. `ϑ`/`ϑ₁`/`ϑ₂` now take the parameters, and the default is resolved
  against the parameters actually passed. `src/harmonic_oscillator.jl`.

### Changed

The three entries below alter a public signature, but each one repairs behaviour that never worked
as documented, so they are bug fixes rather than deliberate API changes: nothing could have
depended on the old form doing what it claimed.

- **Requires EulerLagrange 0.5.** That release stops calling `Symbolics.simplify` on the
  Lagrangian/Hamiltonian by default, which was never a win: measured over all 13 EulerLagrange-based
  problems here (88 generated functions), `simplify = true` was faster **zero** times, slower up to
  15× (`OuterSolarSystem`'s Hamiltonian force), and cost 17.5 s against 0.75 s to build in
  aggregate. Three local workarounds for the old default are therefore gone: the explicit
  `simplify = false` in `src/linear_wave.jl`, the `simplify = N ≤ 10` size heuristic in
  `src/toda_lattice.jl`, and the explicit setting in `src/outer_solar_system.jl`. EulerLagrange 0.5
  also generates code with common subexpression elimination, so results may differ from 0.4 in the
  last bit.

- **Linear wave**: the number of interior points `N` was carried in `default_parameters()` but
  **silently ignored** — `hamiltonian`, `lagrangian` and the symbolic system all sized themselves
  from the module constant `Ñ = 256`, so passing a different `N` changed nothing at all. It cannot
  be a system parameter: it fixes the number of degrees of freedom and the summation bounds, so it
  does not survive `symbolize`. Following `TodaLattice`, it is now a leading positional argument of
  `hodeproblem`/`lodeproblem`/`hodeensemble`/`lodeensemble` (defaulting to `Ñ`), `hamiltonian` and
  `lagrangian` take it as a trailing argument, and `default_parameters()` is `(μ = 0.6,)`. The
  two-argument `hodeproblem(q₀, p₀)` form still works and infers `N` from `length(q₀) - 2`.
  `src/linear_wave.jl`.
- **Plotting**: `LotkaVolterra3d.plot_traces` and `LotkaVolterra4d.plot_traces` took the parameter
  tuple as their second argument, while the Lotka-Volterra 2D, massless-charged-particle and
  point-vortex versions took the problem — an inconsistency across five extensions that are
  otherwise identical in shape. All five now take the problem, which lets them read both the `:h`
  invariant and the parameters from it and allows dispatch on the problem type. Call
  `plot_traces(sol, ode)` instead of `plot_traces(sol, parameters(ode))`.
  `ext/LotkaVolterra3dPlots.jl`, `ext/LotkaVolterra4dPlots.jl`.
- **Perturbed pendulum**: the derived perturbation coefficient
  `A = 0.3ϵ sin(2ϕ) + 0.7ϵ sin(3ϕ)` was stored in `default_parameters()` next to `ϵ` and `ϕ`, so
  overriding either silently kept a stale `A` and produced an inconsistent system. It is computed
  from `ϵ` and `ϕ` by the new `perturbation(params)`, which `symbolize` handles, so the symbolic
  Hamiltonian and Lagrangian now track the dependence correctly. `src/perturbed_pendulum.jl`.

### Documentation
- **Outer solar system**: the module described "the five outer planets (Jupiter, Saturn, Uranus,
  Neptune) and Pluto" — four planets and Pluto. Also documented that building the problem is very
  slow (see *Known follow-ups*). `src/outer_solar_system.jl`.
- Moved the `plot_*` docstrings of the point-vortex and massless-charged-particle extensions from
  the private `_plot_*` helpers onto the public `Module.plot_*` methods, so that `?` and
  `@autodocs` find them. `ext/PointVorticesPlots.jl`, `ext/MasslessChargedParticlePlots.jl`.
- Noted in `src/kubo_oscillator.jl` why the ensembles assemble `EnsembleProblem` by hand:
  GeometricEquations exports `SDEEnsemble`/`PSDEEnsemble`/`SPSDEEnsemble` only as type aliases,
  without the convenience constructors that `ODEEnsemble` and friends have.
- Restructured this changelog into released versions with a migration guide (below).

### Tests
- **`test/outer_solar_system_tests.jl` is now part of `runtests.jl`**, which the note at the end of
  that file had deferred until the problem stopped taking minutes to construct. It runs in ~12 s and
  gained a second testset checking the hand-written vector fields against the generated ones and
  `∇V!` against an independent body-by-body double sum.
- **`test/lode_wiring_tests.jl`**: each problem now carries whether its Lagrangian is regular, and
  the expected `ω` is sized accordingly (``2n × 2n`` regular, ``n × n`` degenerate). Added an
  `any(!iszero, Ω)` assertion — antisymmetry alone is satisfied vacuously by a zero matrix, which is
  exactly what EulerLagrange used to return.
- **`test/lode_wiring_tests.jl`** (new, wired into `runtests.jl`): asserts for *every*
  `lodeproblem`/`ldaeproblem` in the package that `functions(prob).ω` accepts the velocity slot the
  LODE/LDAE evaluate it with and yields an antisymmetric matrix of the right size, and that
  `functions(prob).l` is a finite scalar-valued Lagrangian. No integrator in GeometricIntegrators
  evaluates `ω`, and `check_methods` skips it, so both an `ω`/`l` swap and a missing five-argument
  method are otherwise completely silent — this is the direct regression test for the three `ω`
  bugs above.
- **Harmonic oscillator**: added a non-unit-mass testset (exact solution in both the velocity and
  momentum conventions, degenerate `ϑ`/`p₀` consistency, `Ω = [0 m; -m 0]`, `Ω v = -∇H`) and a
  DAE-energy-constraint testset (`ϕ(t₀) = 0` for a non-default initial condition across all five
  DAE forms). `test/harmonic_oscillator_tests.jl`.
- **`test/default_timespan_timestep_tests.jl`**: completed the module list, which was missing the
  seven modules added after the #83/#82 fix (`HenonHeilesPotential`, `OuterSolarSystem`,
  `PerturbedPendulum`, `MasslessChargedParticleSingular` and the three `LotkaVolterra2d` gauge
  variants) — precisely the newest code the guard was meant to cover.
- Explained the Lotka-Volterra 3D tolerance relaxation of 0.7.0 in place.
  `test/lotka_volterra_3d_tests.jl`.

### Repository hygiene
- Archived the 0.7.0 audit working files (`bugs.md`, `plan.md`, `memory.md`) under `docs/dev/`,
  with a `README.md` explaining what they are. They were agent working logs sitting in the package
  root and shipping in the released tarball.

### Known follow-ups
- `LotkaVolterra3d` declares only `h` as an invariant, although `casimir` is conserved too;
  declaring it would make `Diagnostics.plot_invariant_error(sol; invariant = :c)` work out of the
  box.
- No documentation navigation pages yet for some models; see `docs/make.jl`.

## [0.7.3] — 2026-07-25

### Added
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

### Tests
- Smoke tests for the new diagnostics and point-vortex plots, and documentation examples of
  `plot_invariant_error` on the point vortices' second invariant. `test/plots_tests.jl`,
  `docs/src/diagnostics.md`.

## [0.7.2] — 2026-07-24

### Bug fixes
- **Lotka-Volterra 2D**: `lodeproblem` passed `lagrangian` and `lotka_volterra_2d_ω` to
  `LODEProblem` in the wrong positional order. (The same bug in `ldaeproblem`/`ldaeproblem_slrk`
  was missed; it is fixed under *Unreleased* above.) `src/lotka_volterra_2d_equations.jl`.

### Added
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

### Tests
- **Massless charged particle two-form and Lagrangian forms**: added testsets to
  `test/massless_charged_particle_tests.jl` and `test/massless_charged_particle_singular_tests.jl`
  that compare `ω` against `ForwardDiff.jacobian` of `ϑ` (validating both the assembly and the
  hand-coded `dϑᵢdxⱼ`), check antisymmetry, `Ω₁₂ = −B`, `Ω v = −∇H`, gauge invariance of `Ω` across
  the two modules, and `∂L/∂v = ϑ`; and that integrate `lodeproblem` with
  `MidpointProjection(VPRKGauss(2))` and `ldaeproblem` with `SLRKLobattoIIIE(2)` against the
  `Gauss(8)` ODE reference. Since no integrator in GeometricIntegrators evaluates `ω` and
  `check_methods` skips it, the tests also assert the positional wiring via `functions(prob).ω` /
  `.l`. (Generalised to the whole package under *Unreleased*.)

## [0.7.1] — 2026-07-24

### Added
- **`MasslessChargedParticleSingular`**: the massless charged particle in a gauge whose vector
  potential has a single non-vanishing component (`A₂ = 0`), as required by the degenerate
  variational integrator (DVI) method. Same magnetic and electric fields, and hence the same
  dynamics, as `MasslessChargedParticle`; the shared equations moved to the new
  `src/massless_charged_particle_common.jl`, included by both modules.
  `src/massless_charged_particle_singular.jl`.
- **`HenonHeilesPotential`**: the Hénon-Heiles system as an EulerLagrange-generated
  `hodeproblem`/`lodeproblem`. `src/henon_heiles_potential.jl`.
- **`OuterSolarSystem`**: gravitational N-body model of the Sun, the four outer planets and Pluto
  (18 degrees of freedom), with initial data and constants from Hairer, Lubich & Wanner,
  *Geometric Numerical Integration*, Chapter I.2.4. `src/outer_solar_system.jl`.
- **`PerturbedPendulum`**: a mathematical pendulum with a non-separable perturbation.
  `src/perturbed_pendulum.jl`.

### Bug fixes
- Pinned plot time-axis limits to the downsampled time range, and fixed the Lotka-Volterra 2D and
  massless-charged-particle invariant-error helpers. Added `xlims` to the harmonic-oscillator and
  pendulum traces.

### Changed
- Requires GeometricSolutions v0.6.5.

## [0.7.0] — 2026-07-13

The bulk of this release is a correctness audit of every model in the package, together with a
migration of all plotting from Plots.jl/`RecipesBase` recipes to Makie package extensions and a
standardisation of the problem-construction interface. The audit's findings report, remediation
plan and execution log are archived under [`docs/dev/`](docs/dev/); the bug IDs (`B*`/`C*`) and
plan IDs (`P*`) below refer to them.

### ⚠️ Migration guide

| Before | After |
| --- | --- |
| `default_parameters` (a constant) | `default_parameters()` (a function; `default_parameters(T)` for a given float type) |
| `MyProblem.timespan`, `MyProblem.timestep` (constants) | `MyProblem.DEFAULT_TIMESPAN`, `MyProblem.DEFAULT_TIMESTEP` |
| `odeproblem(Float32)`, `hodeproblem(q₀, p₀, Float32)` | `odeproblem(Float32.(x₀); parameters = default_parameters(Float32))` |
| `Pendulum.hamiltonian(t, q, p)` | `Pendulum.hamiltonian(t, q, p, params)` |
| `LorenzAttractor.lorenz_attractor_ode()` | `LorenzAttractor.odeproblem()` |
| `KuboOscillator.kubo_oscillator_sde_1()` / `_2()` | `KuboOscillator.sdeproblem()` |
| `KuboOscillator.kubo_oscillator_sde_3()` | `KuboOscillator.sdeensemble()` |
| `kubo_oscillator_psde_*` / `_spsde_*` | `psdeproblem`/`psdeensemble`, `spsdeproblem`/`spsdeensemble` |
| `KuboOscillator` `SDEProblem(m, n, v, B, …)` | `SDEProblem(v, B, KuboNoise(), …)` (noise-object API) |
| `PointVortices.ϑ1` … `ϑ4` | `ϑ₁` … `ϑ₄` (likewise `PointVorticesLinear`) |
| `PointVortices.idoeproblem_dg` | `PointVortices.iodeproblem_dg` (the export was always spelled this way) |
| `LotkaVolterra3d` parameters `A1, A2, A3, B1, B2, B3` | `A₁, A₂, A₃, B₁, B₂, B₃` |
| `hamiltonian_iode`, `hamiltonian_pode` | `hamiltonian` (now has the `(t, q, p, params)` method) |
| `@userplot` recipes (`PlotEnergyError`, `Plot_Lotka_Volterra_2d`, …) | `plot_energy_error`, `plot_phase_portrait`, … returning a Makie `Figure` |

`LotkaVolterra3d`'s vector field and Hamiltonian changed sign convention in the `b`-terms (see
*Model fixes*), so its `reference_solution` and any stored results change accordingly.

### Bug fixes
- **Lotka-Volterra 2D** (B10/P10): `lotka_volterra_2d_v_dae` now reads the velocities from
  `q[3]`/`q[4]` instead of its own uninitialised output `v[3]`/`v[4]`.
  `src/lotka_volterra_2d_common.jl`.
- **Lotka-Volterra 3D/4D plots** (B11/P11): the plot recipes' `compute_energy_error` now resolves
  `hamiltonian` (imported from the sibling model module) instead of raising `UndefVarError`.
  (These recipe files were subsequently replaced by the Makie extensions, below.)
- **Bump initial condition** (B12/P12): the vector method of `∂s` now calls `∂s` (was `s`), so it
  returns derivatives rather than values. `src/bump_initial_condition.jl`.
- **Lorenz attractor** (B13/P12): removed the export of the never-defined `plot_lorenz_attractor`.
  `src/lorenz_attractor.jl`.
- **Harmonic oscillator** (B1/P6): fixed the `hamitlonian` typo in `oscillator_pdae_ϕ` and made the
  constraint energy-relative (`H − H₀`), so the `pdae`/`hdae`/`idae`/`ldae` constraints no longer
  raise `UndefVarError`. `src/harmonic_oscillator.jl`.
- **Three-body problem** (B4/P3): repaired `lodeproblem`, which referenced undefined
  `θ₀`/`p₀`/`θ̇` and used `lagrangian_variables(2)`; now uses the module ICs and
  `lagrangian_variables(6)`. `src/three_body_problem.jl`.
- **Massless charged particle** (B6/P4): fixed the argument order in the (dead) `ψ` constraint
  (`g₁(v,t,q,…)`). `src/massless_charged_particle.jl`.
- **Massless charged particle**: added `hamiltonian(t, q, p, params)`, which both massless modules
  were missing. All their implicit problems declare `invariants=(h=hamiltonian,)`, but
  GeometricEquations evaluates that invariant with an additional momentum/velocity slot, so
  e.g. `compute_invariant_error(t, q, p, params, invariants(iodeproblem())[:h])` raised
  `MethodError: no method matching hamiltonian(::Float64, ::Vector, ::Vector, ::NamedTuple)`.
  Mirrors the existing shim in Lotka-Volterra 2D/4D and the harmonic oscillator.
  `src/massless_charged_particle_common.jl`.
- **Linear wave** (B7/P5): `lodeproblem` now calls `LODEProblem` (was a lowercase self-call →
  `MethodError`). `src/linear_wave.jl`.
- **`timespan`/`timestep` shadowing** (#83, #82): PR #81 had renamed the module-level default
  constants `tspan`/`tstep` to `timespan`/`timestep`, shadowing the `GeometricBase.timespan` /
  `GeometricBase.timestep` functions inside every problem module. Renamed to
  `DEFAULT_TIMESPAN`/`DEFAULT_TIMESTEP`.

### Model fixes
- **Coupled harmonic oscillator** (B2/P1): corrected the Lagrangian kinetic term
  `q̇²/(2m)` → `m·q̇²/2`, so the `LODE` and `HODE` describe the same system (previously differed by
  a factor `m²`). `src/coupled_harmonic_oscillator.jl`.
- **Three-body problem** (B3/P2): added the missing body-1↔body-3 gravitational interaction to the
  potential `V` (was only pairs (1,2) and (2,3)). `src/three_body_problem.jl`.
- **Three-body problem**: fixed the Lagrangian kinetic term — it reused the momentum-form
  `T(q̇)=q̇²/(2m)` instead of `m·q̇²/2`, so the LODE and HODE disagreed for non-unit masses (the same
  class of bug as the coupled oscillator); also added a `v̄` momentum→velocity map to the
  `lodeproblem`. `src/three_body_problem.jl`.
- **Massless charged particle** (B5/P4): flipped the sign-reversed `odeproblem` vector field (and
  the documented Hamiltonian-form matrix) so the ODE agrees with the variational IODE/IDAE forms;
  the one-form `ϑ = A` is retained as authoritative. The corrected field is the physical
  **E**×**B** drift `v = (E₂/B, −E₁/B)`; the previous one was the time-reversed flow.
  `src/massless_charged_particle.jl`.
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

### Added
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
  generic `Diagnostics` ones.
- **Harmonic-oscillator / pendulum plotting helpers**: extended the `HarmonicOscillatorPlots` and
  `PendulumPlots` extensions with `plot_phase_portrait(sol)`, `plot_traces(sol)`, and
  `plot_hamiltonian(; …)` (energy-landscape contour), alongside the existing animation-frame
  `plot_solution`. `plot_phase_portrait(sol)` is now a uniform single-argument function across all
  problem extensions. `ext/{HarmonicOscillatorPlots,PendulumPlots}.jl`.
- **Nonlinear oscillators** (B15/P13): implemented the previously empty stub modules
  **Duffing**, **Lennard-Jones**, **Morse**, and **Mathews-Lakshmanan** oscillators. Each provides
  `hamiltonian`, `lagrangian`, and EulerLagrange-generated `hodeproblem`/`lodeproblem` constructors,
  with documentation of the standard textbook forms.
  `src/{duffing,lennard_jones,morse,mathews_lakshmanan}_oscillator.jl`.
- **EulerLagrange ensembles** (#64): `hodeensemble`/`lodeensemble` for the remaining
  EulerLagrange-based problems (double pendulum, three-body, linear wave and the four nonlinear
  oscillators), accepting either a vector of initial conditions or a vector of parameter sets.

### Changed
- **Standardised problem interface**: every problem module now exposes `DEFAULT_TIMESPAN`,
  `DEFAULT_TIMESTEP` and a `default_parameters(::Type{T} = Float64)` function, and every
  `hamiltonian`/`lagrangian` method takes `params`. The `::Type{T}` positional constructor
  overloads of `HarmonicOscillator` and `Pendulum` and the parameter-free `Pendulum.hamiltonian`
  convenience methods are gone; see the migration guide above.

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
- **Documentation figures now use the extensions**: the phase-portrait / trajectory figures on the
  Lotka-Volterra 2D/3D/4D, massless-charged-particle, harmonic-oscillator and pendulum pages, the
  two Hamiltonian energy-landscape figures, and the pendulum time-trace figure now call the
  extension plotting functions instead of hand-rolling CairoMakie. `diagnostics.md` gained worked
  examples of `Diagnostics.plot_energy_error` and `plot_constraint_error`.
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
- **Plotting extensions**: added `test/plots_tests.jl` (wired into `test/runtests.jl`) — smoke
  tests that each migrated plot function builds a Makie `Figure` for a short integrated solution,
  covering the `Diagnostics`, Lotka-Volterra 2D/3D/4D, and massless-charged-particle plots. Added
  `CairoMakie` to `test/Project.toml`. Plots were previously untested.
- Added tests (P15): nonlinear oscillators and three-body problem (HODE↔LODE trajectory
  agreement), massless charged particle (`Ω·v = −∇H`, i.e. ODE ↔ variational sign consistency —
  the direct regression for B5), and linear wave (constructor regression for B7).
  `test/{nonlinear_oscillators,three_body,massless_charged_particle,linear_wave}_tests.jl`.
- **Kubo oscillator**: migrated the SDE/PSDE/SPSDE constructors to the GeometricEquations
  noise-object API (added a `KuboNoise <: AbstractStochasticProcess` marker; the old
  `SDEProblem(m, n, v, B, …)` positional form is gone). The single-IC variants build problems and
  the multi-IC variants build ensembles. The Kubo test is **re-enabled** and rewritten to the new
  API (it now checks the drift/diffusion terms directly).
  `src/kubo_oscillator.jl`, `test/kubo_oscillator_tests.jl`, `test/runtests.jl`.
- Re-calibrated the Lotka-Volterra 3D test tolerances to the corrected dynamics (P8), and refreshed
  the stale `reference_solution`. `test/lotka_volterra_3d_tests.jl`.
- Added regression tests for the `timespan`/`timestep` shadowing (#83, #82) and for the
  EulerLagrange ensembles (#64). `test/default_timespan_timestep_tests.jl`,
  `test/eulerlagrange_ensembles_tests.jl`.

### Repository hygiene
- **Dropped plotting-only dependencies**: with the recipes gone, `RecipesBase`, `Measures`, and
  `LaTeXStrings` were removed from `[deps]`/`[compat]` (they were used only by the deleted plot
  files). Deleted `src/plot_recipes.jl` and
  `src/{lotka_volterra_2d,lotka_volterra_3d,lotka_volterra_4d,massless_charged_particle}_plots.jl`,
  and removed their `include`s from `src/GeometricProblems.jl`. Registered the five new extensions in
  `Project.toml`. `docs/src/lotka_volterra_2d.md`'s stale `@autodocs` on the removed
  `LotkaVolterra2dPlots` module now resolves through the problem module (docstrings live on the
  extended functions).
- **Docs plotting → CairoMakie**: migrated all documentation plotting from Plots.jl and GLMakie.jl
  to CairoMakie.jl, and dropped `Plots` and `GLMakie` from `docs/Project.toml`.
- **`.gitignore`** (P17): added `*-prev.jl` backups and the `notebooks/`, `obsolete/`,
  `prototyping/` and `references/` scratch directories (none of which were tracked). No files were
  deleted.
- **Lotka-Volterra 3D (Lagrangian)** (B14/P16): quarantined the unregistered, broken
  `src/lotka_volterra_3d_lagrangian.jl` — renamed its module `LotkaVolterra3d` →
  `LotkaVolterra3dLagrangian` to remove the name collision with the registered module, and added a
  warning docstring. It remains unregistered (a 3-D Lotka-Volterra admits no non-degenerate
  Lagrangian); not deleted.

[Unreleased]: https://github.com/JuliaGNI/GeometricProblems.jl/compare/v0.7.3...HEAD
[0.7.3]: https://github.com/JuliaGNI/GeometricProblems.jl/compare/v0.7.2...v0.7.3
[0.7.2]: https://github.com/JuliaGNI/GeometricProblems.jl/compare/v0.7.1...v0.7.2
[0.7.1]: https://github.com/JuliaGNI/GeometricProblems.jl/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/JuliaGNI/GeometricProblems.jl/compare/v0.6.25...v0.7.0
