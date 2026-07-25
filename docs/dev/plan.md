# GeometricProblems.jl — Remediation Plan

Companion to [`bugs.md`](bugs.md). Bug IDs (B1–B15, C1–C4) refer to that report.

Author decisions incorporated:
- **Massless charged particle (B5):** the one-form `ϑ = A` is correct, so the variational
  IODE/IDAE forms are authoritative. (This supersedes the earlier "documented Hamiltonian form is
  correct" reading: the `odeproblem` vector field and the documented Hamiltonian-form sign are the
  ones to fix.)
- **Lotka-Volterra 3D (C2):** `H = a·q + b·ln q` (the docstring) is authoritative.
- **Harmonic oscillator (C1):** the mass `m` must be used consistently.
- **No deletions** of orphan doc pages or empty stub modules — fill/implement instead.

Recommended order: Tier 1 (wrong numerics / crashes) → Tier 2 (form & doc consistency) →
Tier 3 (completeness: stubs, docs, tests, hygiene). Each item lists the fix and how to verify it.

---

## Tier 1 — correctness bugs

### P1. Coupled harmonic oscillator kinetic term (B2)
`src/coupled_harmonic_oscillator.jl:64` — change `q̇[i]^2 / (2*mᵢ)` to `mᵢ * q̇[i]^2 / 2` for
both DOF, so `L = Σ mᵢ q̇ᵢ²/2 − V`. This makes the Legendre transform give `p = m q̇` and
`H = p²/(2m) + V`, matching the coded `hamiltonian` (`:57`) and `v̄` (`:68`).
Also fix the docstring (`:12`): kinetic terms are `p₁²/(2m₁) + p₂²/(2m₂)`.
**Verify:** integrate `hodeproblem()` and `lodeproblem()` from the same IC; assert trajectory
agreement to ~1e-10 (currently they diverge by `m²`).

### P2. Three-body missing interaction (B3)
`src/three_body_problem.jl:55` — add the (1,3) pair to `V`:
`− params.G * params.m₁ * params.m₃ / √((q[1] − q[5])^2 + (q[2] − q[6])^2)`.
**Verify:** with three equal masses on a symmetric configuration, the net force on each body is
nonzero and the total momentum/energy are conserved along an integration.

### P3. Three-body `lodeproblem` (B4)
`src/three_body_problem.jl:106-110` — repair by mirroring the working `hodeproblem`: use the
module's own `q₀`/`p₀`, `v̄` = the module's velocity function, and `lagrangian_variables(6)`.
**Verify:** `lodeproblem()` constructs without error and its trajectory matches `hodeproblem()`.

### P4. Massless charged particle sign (B5)
Keep the one-form `ϑ = A` and the derived `dϑ` (`:98-101`), `ω!`, and `f` (`:104-105, :142-146`)
**unchanged** — the `IODE`/`IDAE` variational forms are correct (they give the physical E×B drift
`v = (E₂/B, −E₁/B)`). Correct the sign-reversed forms to match:
- `odeproblem` vector field (`:78-79`): flip both signs → `v₁ = +E₂/B`, `v₂ = −E₁/B`.
- Documentation (`massless_charged_particle.md` and the module docstring): the Hamiltonian form
  should read `ẋ = (1/B)[[0,-1],[1,0]]∇φ` (i.e. flip the current sign), so it equals
  `(E₂/B, −E₁/B)`, consistent with `ω = dϑ` for `ϑ = A`.
Separately fix the dead `ψ` argument order (B6): call `g₁(v,t,q,params)`, `g₂(v,t,q,params)`.
**Verify:** integrate `odeproblem()` and `iodeproblem()` from the same IC; assert trajectory
agreement to ~1e-10. (This is the direct regression test for the original defect.)

### P5. Linear wave `lodeproblem` (B7)
`src/linear_wave.jl:69` — call `LODEProblem(lag_sys, timespan, timestep, q₀, p₀; parameters)`
(capital, from EulerLagrange), matching the `HODEProblem` path (`:59`).
**Verify:** `lodeproblem()` constructs without error and its trajectory matches `hodeproblem()`.

### P6. Harmonic oscillator `hamitlonian` typo + constraint (B1)
`src/harmonic_oscillator.jl:412` — `hamitlonian` → `hamiltonian`, and return
`hamiltonian(t,q,p,params) − hamiltonian(t₀,q₀,p₀,params)` to match `oscillator_dae_ϕ` (`:368`).
**Verify:** add a smoke test that evaluates the `pdae`/`hdae`/`idae`/`ldae` constraint arrays
(they should be finite and ≈0 at the initial condition).

---

## Tier 2 — form / documentation consistency

### P7. Harmonic oscillator — thread mass `m` (C1)
Make every vector field and the exact/reference solution consistent with `H = p²/(2m) + kq²/2`
(the DELE derivatives already carry `m`, so leave them). Concretely:
- `oscillator_ode_v` (`:143`): `v[2] = −(k/m) x[1]`.
- `oscillator_pode_v` / `oscillator_pdae_v` / `oscillator_iode_v` (`:171,:379,:288-291`): `q̇ = p/m`.
- `oscillator_iode_ϑ` (`:272-275`): `p = m v`; degenerate forms (`:314-338`) likewise carry `m`.
- `oscillator_dae_u` (`:364`) and `oscillator_pdae_ψ` (`:418`): include the `1/m` / `m` factors.
- `oscillator_sode_*` (`:237-261`): `−(k/m)` in the momentum flow.
- `A`/`ϕ`/`exact_solution` (`:88-92`): use `m ω²` (i.e. `1/(m·ω²)` factors), and `p = m q̇`.
**Verify:** run every form (ODE/PODE/HODE/IODE/LODE/SODE/DAE/DELE) with a non-unit mass
(e.g. `m = 2`) from a common IC and assert mutual trajectory agreement and energy conservation;
compare against `exact_solution` to ~1e-10.

### P8. Lotka-Volterra 3D — adopt the documented Hamiltonian (C2)
`src/lotka_volterra_3d.jl:86-89` — change `hamiltonian` to `H = a·q + b·log q`
(`+ B1*log(q[1]) + B2*log(q[2]) + B3*log(q[3])`), and flip the `b`-term signs in the vector field
`v₁,v₂,v₃` (`:63-76`) to match the docstring equations (`:9-11`).
**Verify:** assert `v == P∇H` component-wise (symbolically or numerically) with the new `H`, and
that energy + the Casimir `q₁q₂q₃` are conserved along an integration. Update the tight
`reference_solution` if the tests rely on it.

### P9. Lotka-Volterra 2D gauge variant (B9)
`src/lotka_volterra_2d_gauge.jl` — implement the documented gauge one-form: add `+q[2]` to `ϑ₁`
and `+q[1]` to `ϑ₂`, and add `+1` to `dϑ₁dx₂` and `dϑ₂dx₁`. This makes it a genuinely distinct,
gauge-equivalent Lagrangian.
**Verify:** the ODE trajectory is unchanged (gauge terms are an exact one-form, so `ω` is
unchanged, `= −1/(q₁q₂)`), but the variational integrator now differs from the symmetric variant
(the gauge vs symmetric test tolerances should no longer be identical).

### P10. Lotka-Volterra 2D DAE velocity (B10)
`src/lotka_volterra_2d_equations.jl:126-132` — `v[1] = q[3]`, `v[2] = q[4]`.
**Verify:** construct and integrate `daeproblem()`; energy/Casimir conserved (add a small test).

### P11. Plot-recipe energy error (B11)
`src/lotka_volterra_3d_plots.jl` / `_4d_plots.jl` — route `compute_energy_error` through
`invariants(equ)[:h]` (as the 2D plots do), or import `hamiltonian` from the model module.
**Verify:** call the plot recipe on a solution without `UndefVarError`.

### P12. Doc/formula fixes (C3, C4, B12, B13)
- `src/lotka_volterra_4d_lagrangian.jl:9` — remove the stray `b₅ log q⁵` term.
- `src/bump_initial_condition.jl:39` — vector `∂s` should call `∂s`, not `s`.
- `src/lorenz_attractor.jl:9` — remove the undefined `plot_lorenz_attractor` export
  (or define the recipe).
- `docs/src/linear_wave.md` — reconcile the stated `H_h` with the coded `H`, correct the BC
  description (boundary nodes are dynamical), fix the LaTeX typos.

---

## Tier 3 — completeness, documentation, tests, hygiene

### P13. Implement the nonlinear-oscillator stubs (B15) — do NOT delete
Implement `duffing_oscillator.jl`, `lennard_jones_oscillator.jl`,
`mathews_lakshmanan_oscillator.jl`, `morse_oscillator.jl` with their standard forms, e.g.:
- Duffing: `H = p²/2 + α q²/2 + β q⁴/4` (add damping/forcing only if a non-autonomous form is wanted).
- Lennard-Jones: `V(r) = 4ε[(σ/r)¹² − (σ/r)⁶]`, `H = p²/2 + V(r)`.
- Mathews–Lakshmanan: `L = ½(q̇² − ω²q²)/(1 + λ q²)`.
- Morse: `V(r) = D(1 − e^{−a(r−r₀)})²`, `H = p²/(2m) + V(r)`.
Follow the harmonic-oscillator module as the structural template; export only names that are
defined. Until implemented, at minimum the premature `export hamiltonian` should not reference an
undefined symbol (define `hamiltonian` first). Add matching content to `nonlinear_oscillators.md`.

### P14. Fill the orphan / stub doc pages (§4 of bugs.md) — do NOT delete
Add content (and, where sensible, backing models) for `kepler_problem.md`, `henon_heiles.md`,
`inner_solar_system.md`, `outer_solar_system.md`, `diagnostics.md`; and flesh out the thin pages
(`lotka_volterra_4d.md`, `point_vortices.md`, `lorenz_attractor.md`, a page for `kubo_oscillator`
and `point_vortices_linear`, and the empty sections in `harmonic_oscillator.md`). Populate empty
module docstrings (`lotka_volterra_4d.jl`, `lorenz_attractor.jl`) so the autodoc pages render.
Fix the `rigid_body.md` example ICs and note the Casimir; add `m₃`/`G` to the three-body docstring.

### P15. Tests (§5 of bugs.md)
Add test files for `linear_wave`, `massless_charged_particle`, `three_body_problem`; re-enable the
Kubo test (`test/runtests.jl:18`). New tests should include the cross-form trajectory checks used
in P1–P8 (ODE vs HODE/IODE/LODE agreement, energy/Casimir conservation) plus a DAE-constraint and
plot-recipe smoke test (to catch B1/B10/B11).

### P16. `lotka_volterra_3d_lagrangian.jl` (B14)
Do not delete. Either repair (resolve the `module LotkaVolterra3d` name collision, fix the
`ϑ(t,q)` method, and document that the 3-D LV admits no nondegenerate Lagrangian — so this is an
experimental/degenerate construction) or leave it unregistered with a header comment marking it
WIP. Do not `include` it until it is consistent.

### P17. Hygiene (§6 of bugs.md)
Add the `-prev.jl` backups and the `examples/`, `prototyping/`, `notebooks/`, `references/`,
`harmonic-oscillator-plots/` scratch directories to `.gitignore` (do not delete).

---

## Global verification

- Baseline first: `julia --project=. -e 'using Pkg; Pkg.test()'` (confirm current suite passes,
  and see which newly added tests fail before fixes).
- After each fix, re-run the relevant test; after all fixes, run the full suite.
- For every model with multiple forms, assert ODE ↔ HODE/PODE/IODE/LODE trajectory agreement from
  a common IC to a tight tolerance, and `∇H · ẋ ≈ 0` (energy) plus any Casimir conservation.
- Rebuild the docs: `julia --project=docs docs/make.jl` — confirm no page errors and that
  autodocs render.
