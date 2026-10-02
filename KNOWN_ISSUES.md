# Known issues

Defects found during review and recorded, not fixed.

## KI-1 · Comments outside `test/` name test files that were renamed

- **Kind:** docs (stale comment)
- **Evidence:** the test files were renamed in the unified test layout. These lines still name
  the old names: `benchmark/linear_wave.jl:373`, `examples/point_vortices.jl:8`,
  `examples/massless_charged_particle_2d.jl:8`, `examples/lotka_volterra_2d.jl:10`,
  `examples/lotka_volterra_4d.jl:8`. Comments under `src/` name them too:
  `src/toda_lattice.jl:94`, `src/linear_wave.jl:93`, `src/lotka_volterra_2d_equations.jl:36`
  and `:59`, `src/massless_charged_particle_common.jl:23` and `:40`.
  `grep -rn '_tests\.jl' src benchmark examples` lists them.
- **Fix:** change each name to the new path under `test/`.

## KI-2 · The Euler-Lagrange ensemble tests only check that the members differ

- **Kind:** weak test
- **Evidence:** `test/integration/eulerlagrange_ensembles.jl` asserts only that the ensemble has
  two members whose solutions differ. `mutate.jl` mutants on the Duffing Hamiltonian
  (`src/duffing_oscillator.jl:43`, `β * q[1]^4 / 4` → `/ 3`) and Lagrangian (`:48`, the same
  change) both SURVIVED.
- **Fix:** compare each ensemble member with the matching single problem, or with a reference
  value.

## KI-4 · The test, Doctests and Documentation jobs fail because the test and docs environments do not resolve

- **Kind:** upstream
- **Evidence:** `test/Project.toml` and `docs/Project.toml` require `GeometricIntegrators = "0.18"`,
  and the package requires `GeometricBase = "0.15.0"`. In the General registry every
  GeometricIntegrators release from 0.18.0 to 0.18.5 bounds GeometricBase to 0.14.8–0.14.12. On
  Julia 1.11.9, a temporary environment that develops this repository and adds
  GeometricIntegrators 0.18 fails with `Unsatisfiable requirements detected`. In one run the
  resolver reported it for package GeometricBase [9a0b12b7]: `restricted to versions 0.15 by
  GeometricProblems`, `restricted by compatibility requirements with GeometricIntegrators
  [dcce2d33] to versions: 0.14.8 - 0.14.12 — no versions left`. In another run it reported it for
  package GeometricIntegrators; the package named depends on the resolve order, the cause is the
  same. Every test job of `CI.yml`, its `Doctests - ubuntu-latest` job and the
  `Documentation` job of `Documenter.yml` instantiate one of these environments and fail there.
  The jobs heal when GeometricIntegrators registers a release for GeometricBase 0.15 (0.18.6).
  Then a `workflow_dispatch` of `CI.yml` and of `Documenter.yml` on `main` must be green, and a
  later pull request deletes this entry.
- **Fix:** none in this repository; wait for the GeometricIntegrators release.

## KI-5 · The examples environment bounds a GeometricIntegrators release that is not registered, and two examples call names that 0.18.0 renamed

- **Kind:** not verified
- **Evidence:** `examples/Project.toml` bounds `GeometricIntegrators = "0.18.6"`, the release for
  GeometricBase 0.15, and the environment does not resolve until it registers. The bound moves from
  0.17 to 0.18.6 and skips 0.18.0–0.18.5. GeometricIntegrators 0.18.0
  renames `SymplecticEulerA` and `SymplecticEulerB`, which `examples/harmonic_oscillator.jl:82`
  and `:88` and `examples/pendulum.jl:67` and `:72` call. Through GeometricIntegrators the
  environment also moves from SimpleSolvers 0.10 to 0.14: 0.11 removes `Backtracking`'s `α₀`,
  0.12 stops line-search warnings inside `solver_step!`, and 0.13 changes the default linear
  solver for LAPACK element types. No example calls SimpleSolvers directly, so only the last
  change can reach them, through the implicit integrators. No example was run against the new
  bounds; CI runs none.
- **Fix:** after GeometricIntegrators 0.18.6 registers, resolve the examples environment, run each
  example, and fix the callers in a later pull request.
