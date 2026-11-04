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

## KI-5 · Seven of the nine examples have no recorded run against the bounds of the examples environment

- **Kind:** not verified
- **Evidence:** `examples/Project.toml` bounds `GeometricIntegrators = "0.19"` and resolves with
  GeometricIntegrators 0.19.0. `examples/harmonic_oscillator.jl` runs to the end and writes its five
  figures. `examples/pendulum.jl` fails (KI-6). The other seven were written for
  GeometricIntegrators 0.17 and SimpleSolvers 0.10, and no run against these bounds is recorded.
  The environment resolves SimpleSolvers 0.14 through GeometricIntegrators: 0.11 removes
  `Backtracking`'s `α₀`, 0.12 stops line-search warnings inside `solver_step!`, and 0.13 changes
  the default linear solver for LAPACK element types. No example calls SimpleSolvers directly, so
  only the last change can reach them, through the implicit integrators. CI runs no example.
- **Fix:** run each of the seven examples, and fix what fails in a later pull request.

## KI-6 · `examples/pendulum.jl` calls `hamiltonian` without the parameters

- **Kind:** defect
- **Evidence:** `examples/pendulum.jl:18` and `:27` call `hamiltonian.(t, q, p)`. The method in
  `src/pendulum.jl` takes `(t, q, p, params)`, so the script stops with
  `MethodError: no method matching hamiltonian(::Float64, ::Float64, ::Float64)` at line 58, in the
  first section, before any symplectic Euler call. `examples/harmonic_oscillator.jl` passes
  `Ref(params)`.
- **Fix:** pass `default_parameters()` as `harmonic_oscillator.jl` does, and run the script.
