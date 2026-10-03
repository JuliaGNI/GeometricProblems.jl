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

## KI-5 · The examples environment bounds a GeometricIntegrators release that is not registered, and no example has run against it

- **Kind:** not verified
- **Evidence:** `examples/Project.toml` bounds `GeometricIntegrators = "0.18.6"`, the release for
  GeometricBase 0.15, and the environment does not resolve until it registers. The bound skips
  0.18.0–0.18.5, and the examples are written for GeometricIntegrators 0.17 and SimpleSolvers 0.10.
  GeometricIntegrators 0.18.0 renames its `SymplecticEulerA` and `SymplecticEulerB` to
  `SymplecticEulerARK` and `SymplecticEulerBRK`, and re-exports the GeometricIntegratorsBase
  integrators of the unsuffixed names. So the calls in `examples/harmonic_oscillator.jl:82` and
  `:88` and `examples/pendulum.jl:67` and `:72` resolve, but run the GeometricIntegratorsBase
  integrators instead of the Runge-Kutta ones. Through GeometricIntegrators the environment
  resolves SimpleSolvers 0.14: 0.11 removes `Backtracking`'s `α₀`, 0.12 stops line-search warnings
  inside `solver_step!`, and 0.13 changes the default linear solver for LAPACK element types. No
  example calls SimpleSolvers directly, so only the last change can reach them, through the
  implicit integrators. No example has run against these bounds; CI runs none.
- **Fix:** after GeometricIntegrators 0.18.6 registers, resolve the examples environment, run each
  example, and fix what fails in a later pull request.
