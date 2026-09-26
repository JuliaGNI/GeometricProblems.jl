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

## KI-3 · `test/quality/aqua.jl` has run only on Julia 1.13

- **Kind:** not verified
- **Evidence:** a local run gave 10 Pass and 1 Broken. The ambiguity and persistent-task checks
  can differ on the 1.10 floor. The CI `min` job is the first run there.
- **Fix:** read the `min` job of the pull request's CI.
