# GeometricProblems.jl — Correctness Audit / Bug Report

Audit of every model included in `src/GeometricProblems.jl`, checking (1) consistency
between the different equation forms of each model, (2) agreement with the literature,
(3) other bugs, plus a record of deliberate simplifications and a documentation overview.

All **HIGH** findings were verified by direct inspection of the source. Line numbers are
indicative and may drift as files change.

Legend — Severity:
- **HIGH** — wrong numerics or a crash on a normal call path.
- **LATENT** — crash / wrong result only on a code path not exercised by the current tests.
- **INCONSISTENCY** — forms disagree, or code disagrees with its own docs (often correct only at default parameters).
- **DOC** — documentation defect.
- **SIMPLIFICATION** — deliberate omission of physical constants (recorded, not necessarily a bug).
- **HYGIENE** — repository cleanliness.

---

## 1. Confirmed bugs

### B1 — Harmonic oscillator: `hamitlonian` typo in DAE constraint — HIGH (latent)
`src/harmonic_oscillator.jl:412` — `oscillator_pdae_ϕ` calls `hamitlonian(t, q, p, params)`
(misspelled) → `UndefVarError` whenever the constraint is evaluated. It is wired into
`pdaeproblem`, `hdaeproblem`, `idaeproblem`, and `ldaeproblem`.
Secondary: the constraint returns `H` rather than `H − H₀`; compare the correct
`oscillator_dae_ϕ` (`:368`) which returns `hamiltonian(t,x,params) − hamiltonian(t₀,x₀,params)`.
Uncaught because the tests only construct the problems, never evaluate `ϕ`.

### B2 — Coupled harmonic oscillator: Lagrangian kinetic term inverted in `m` — HIGH
`src/coupled_harmonic_oscillator.jl:64` — kinetic energy is written `q̇²/(2m)`; it must be
`m·q̇²/2`. Legendre-transforming the coded `L` gives `p = q̇/m ⇒ q̇ = m·p` and
`H = m·p²/2 + V`, which is **not** the coded `H = p²/(2m) + V` (`:57`). Consequently the
`LODE` integrates `q̈ = −m ∂V/∂q` while the `HODE` and `v̄` (`:68`, `q̇ = p/m`) integrate
`q̈ = −(1/m) ∂V/∂q`. With the defaults `m₁=2, m₂=1` these are physically different systems
(differ by `m²`).
Related DOC bug: docstring Hamiltonian (`:12`) writes the kinetic term with `q` instead of `p`.

### B3 — Three-body problem: missing 1↔3 gravitational interaction — HIGH
`src/three_body_problem.jl:55` — the potential contains only the (1,2) and (2,3) pairs:
`V = −G m₁m₂/√((q₁−q₃)²+(q₂−q₄)²) − G m₂m₃/√((q₃−q₅)²+(q₄−q₆)²)`.
The (1,3) term `−G m₁m₃/√((q₁−q₅)²+(q₂−q₆)²)` is missing, so bodies 1 and 3 do not interact —
this is a "chain," not a gravitational three-body problem.

### B4 — Three-body problem: broken `lodeproblem` — HIGH (uncaught, no test)
`src/three_body_problem.jl:106-110` — copy-paste from the double pendulum. Defaults
`q₀ = θ₀` and `v̄ = θ̇` reference names that do not exist in the module, and
`lagrangian_variables(2)` should be `(6)`. Calling `lodeproblem()` throws `UndefVarError`.
(`hodeproblem` is fine.)

### B5 — Massless charged particle: `odeproblem` sign-reversed vs the variational (IODE/IDAE) form — HIGH (uncaught, no test)
`src/massless_charged_particle.jl` — the `IODE`/`IDAE` Euler–Lagrange forms are built from the
**correct one-form `ϑ = A`** (`:82-85`, with the derived `dϑ`/`ω!`/`f`) and yield the physical
E×B guiding-centre drift `v = (E₂/B, −E₁/B)`. These are the authoritative dynamics. The
`odeproblem` vector field (`:78-79`, `v₁ = −E₂/B, v₂ = +E₁/B`) is the **exact negation**
(time-reversed), so `odeproblem` and `iodeproblem` integrate mirror-image trajectories. The
documented Hamiltonian form (`ẋ = (1/B)[[0,1],[-1,0]]∇φ`) carries the same reversed sign as the
ODE. Energy `H = φ` is conserved either way, so the invariant does not detect it. Fix is on the
ODE (and documentation) side; keep `ϑ = A` (see plan P4).

### B6 — Massless charged particle: wrong argument order in `ψ` — LATENT
`src/massless_charged_particle.jl:176-177` — `massless_charged_particle_ψ` calls
`g₁(t,q,v,…)`/`g₂(t,q,v,…)` but `g₁`,`g₂` are defined as `g₁(v,t,q,…)` (`:107-108`) — arguments
shuffled. Dead code (`ψ` is never wired into a problem constructor), but incorrect.

### B7 — Linear wave: `lodeproblem` calls a lowercase constructor → `MethodError` — HIGH (uncaught, no test)
`src/linear_wave.jl:69` — recursively calls lowercase `lodeproblem(lag_sys, …)` instead of
EulerLagrange's `LODEProblem` (capital). It resolves to the module's own
`lodeproblem(q₀, p₀; …)`, finds no matching method for `(::LagrangianSystem, ::Tuple, …)`, and
throws. The `HODE` path correctly uses `HODEProblem` (`:59`).

### B8 — Point vortices: exported name unbound (`idoeproblem_dg` typo) — MEDIUM
`src/point_vortices.jl` — the export list (`:11`) names `iodeproblem_dg`, but the function is
defined `idoeproblem_dg` ("ido" transposed, `:244`). The exported name is unbound and the DG
constructor is unreachable. It also passes `v=point_vortices_v` where `v̄=` is intended.
(`point_vortices_linear.jl` spells it correctly.)

### B9 — Lotka-Volterra 2D "gauge" is a byte-identical copy of "symmetric" — MEDIUM (inconsistency)
`src/lotka_volterra_2d_gauge.jl:18-25` — the module body is identical to
`lotka_volterra_2d_symmetric.jl` except for the module name (verified by diff). The documented
gauge terms (`+q₂`, `+q₁` added to `ϑ`, and the matching `+1` in `dϑ₁dx₂`/`dϑ₂dx₁`) were never
applied, so `LotkaVolterra2dGauge` yields the *same* variational integrator as the symmetric
variant — defeating the module's stated purpose. ODE dynamics are still correct (shared
`equations.jl`), so the ODE test passes and the gauge/symmetric test tolerances are identical.

### B10 — Lotka-Volterra 2D: DAE velocity reads its own output — LATENT
`src/lotka_volterra_2d_equations.jl:126-132` — `lotka_volterra_2d_v_dae` sets `v[1]=v[3]`,
`v[2]=v[4]`, reading the uninitialized output array; it should read `q[3]`, `q[4]` (the
velocities, per `ϕ_dae` at `:192-196`). Latent: the DAE is constructed but never integrated.

### B11 — Lotka-Volterra 3D/4D plot recipes reference undefined `hamiltonian` — LATENT
`src/lotka_volterra_3d_plots.jl:11`, `src/lotka_volterra_4d_plots.jl:8` — `compute_energy_error`
closes over a bare `hamiltonian` that is never imported from the model module → `UndefVarError`
when the recipe runs. (The 2D plots correctly route through `invariants(equ)[:h]`.)

### B12 — `bump_initial_condition`: vector `∂s` returns values, not derivatives — LATENT
`src/bump_initial_condition.jl:37-40` — the vector method of `∂s` calls `s` instead of `∂s`.
Latent: the Toda / linear-wave default IC path uses the (correct) scalar `∂s`, but any vector
caller (or a non-zero-`p₀` route) gets wrong results.

### B13 — Lorenz attractor: exported `plot_lorenz_attractor` is never defined — LOW
`src/lorenz_attractor.jl:9` — exports a symbol that has no definition anywhere.

### B14 — `lotka_volterra_3d_lagrangian.jl`: unregistered and broken — DEAD CODE
`src/lotka_volterra_3d_lagrangian.jl` is **not** `include`d in `src/GeometricProblems.jl`, and:
- `module LotkaVolterra3d` (`:4`) name-collides with `src/lotka_volterra_3d.jl` (`:40`).
- `ϑ(t, q)` (`:92-96`) calls a method signature that does not exist → `MethodError`.
- The coded `v` does not satisfy `Ω v = −∇H` for the active `ϑ`.
- Mathematically, the 3-D Lotka-Volterra is an odd-dimensional Poisson system and does **not**
  admit this nondegenerate Lagrangian formulation (a 3×3 antisymmetric `Ω` is singular).
Uses stale API (`h=` keyword instead of `invariants=(h=…,)`). Empty docstring.

### B15 — Four nonlinear-oscillator modules are empty stubs exporting an undefined name — INCOMPLETE
`src/duffing_oscillator.jl`, `src/lennard_jones_oscillator.jl`,
`src/mathews_lakshmanan_oscillator.jl`, `src/morse_oscillator.jl` — each is a 9-line module that
`export hamiltonian` but never defines it (any access → `UndefVarError`) and implements no
dynamics. `docs/src/nonlinear_oscillators.md` autodocs against them and renders empty.

---

## 2. Consistency issues (code correct at defaults, but inconsistent across forms / with docs)

### C1 — Harmonic oscillator: mass `m` dropped from the dynamics — INCONSISTENCY / SIMPLIFICATION
`H`, `L`, the energy invariant, and the DELE derivatives all carry `m`, but every vector field
and the exact/reference solution silently assume `m = 1`:
`oscillator_ode_v` (`:143`, `−k x₁` should be `−(k/m) x₁`), `oscillator_pode_v` (`:171`,
`q̇ = p` should be `p/m`), `oscillator_iode_v`/`oscillator_iode_ϑ` (`:288-291`, `:272-275`,
`p = v` should be `m v`), `oscillator_pdae_v` (`:379`), `oscillator_dae_u` (`:364`),
the degenerate forms (`:314-338`), `oscillator_pdae_ψ` (`:418`), and `A`/`ϕ`/`exact_solution`
(`:88-92`, which use `k` where `m ω²` is needed). All forms agree only at the default `m = 1`.

### C2 — Lotka-Volterra 3D: code Hamiltonian sign differs from the docstring — INCONSISTENCY
Code `hamiltonian` (`src/lotka_volterra_3d.jl:86-89`) is `H = a·q − b·log q`, with `+b` terms in
the vector field; the docstring (`:9-11`, `:28`) states `H = a·q + b·ln q` with `−b` terms. The
code is internally self-consistent (`v = P∇H` verified, energy conserved), but the **documented
equations are not the ones integrated** — they differ for asymmetric `b`, including at the
defaults `B1=0, B2=1, B3=1`. Per the author, the **docstring is authoritative**.

### C3 — Linear wave: coded `H` ≠ documented `H_h`; boundary conditions not enforced — INCONSISTENCY / DOC
Code uses unweighted kinetic `Σ p²/2` and potential coefficient `μ²/(4Δx²)`;
`docs/src/linear_wave.md:13` states a `Δx`-weighted energy. The equations of motion are the
correct semidiscrete wave equation, but the stated Hamiltonian does not equal the coded one.
The doc also claims Dirichlet BC (`q = 0` on `∂Ω`), while the boundary nodes are actually
dynamical (free-end-like: `q̈₁ = μ²(q₂−q₁)/(2Δx²)`).

### C4 — Doc-only formula errors — DOC
- `src/coupled_harmonic_oscillator.jl:12` — kinetic energy written with `q` instead of `p` (see B2).
- `src/lotka_volterra_4d_lagrangian.jl:9` — stray `b₅ log q⁵` term (only 4 components exist).

---

## 3. Deliberate simplifications (recorded)

| Model | Simplification |
|---|---|
| Harmonic oscillator | Mass `m` effectively fixed to 1 in the dynamics & exact solution (see C1; to be fixed). |
| Kubo oscillator | Frequency and mass fixed to 1; only noise intensity `ν` is a parameter. |
| Linear wave | Advertised `N` parameter **ignored** — grid hardcoded `Ñ=256` (`:27,30,36,39,59,69`); wave speed `= μ`. |
| Point vortices (linear) | Separation `d = 1` hardcoded, not a parameter. |
| Three-body | `m₁ = m₂ = m₃ = G = 1` (per `jin2020sympnets`). |
| ABC flow | Default `A=0.5, B=1, C=1` (not the classic `A=B=C=1` / `A=√3,B=√2,C=1`). |
| Rigid body | Only the ODE (Euler equations) provided; Lie–Poisson/Hamiltonian structure not exposed. |
| Coupled harmonic oscillator | Sigmoid-coupled potential is a custom (non-literature) construction, acknowledged in docstring. |
| Pendulum | Defaults `l=m=g=1`; nonstandard but **documented** sign convention (θ from upward vertical). Not a bug. |
| Lotka-Volterra (all) | `NaNMath.log` guard; large blocks of commented-out alternative one-forms/matrices. |

---

## 4. Documentation overview

**Adequate** (states the equations / real prose): `pendulum.md` (excellent),
`double_pendulum.md`, `toda_lattice.md`, `abc_flow.md`, `lotka_volterra_2d.md`,
`massless_charged_particle.md` (but carries the sign issue B5), `initial_condition.md`.

**Incomprehensive** (thin, partial, or contains an error): `harmonic_oscillator.md`
("Lagrangian Formulation" and "Dynamics" are empty headers), `coupled_harmonic_oscillator.md`
(short + `q`/`p` typo), `rigid_body.md` (equations only; example ICs swap sin/cos vs the module),
`three_body_problem.md` (no equations; docstring omits `m₃`, `G`), `lotka_volterra_3d.md`
(math present but sign mismatch C2), `linear_wave.md` (`H` mismatch, wrong BC claim, LaTeX typos:
`\Delta_xi`, "equidistantly spaces points", "an completely-integrable"),
`nonlinear_oscillators.md` (autodocs against four empty modules).

**None / stub**: `lotka_volterra_4d.md` (empty module docstring → near-empty page),
`point_vortices.md` (title only; deformed-`S` variant unexplained), `lorenz_attractor.md`
(empty docstring + autodocs only), `point_vortices_linear` (no page), `kubo_oscillator`
(no page, no docstrings), Duffing / Lennard-Jones / Mathews-Lakshmanan / Morse (empty).

**Orphan doc pages** wired into `docs/make.jl` with **no backing source model** and empty content
(title only): `kepler_problem.md`, `henon_heiles.md`, `inner_solar_system.md`,
`outer_solar_system.md`, `diagnostics.md`.

---

## 5. Test-coverage gaps (these hide several bugs above)

- Kubo oscillator test is **commented out** (`test/runtests.jl:18`).
- **No test at all** for `linear_wave`, `massless_charged_particle`, `three_body_problem`
  (this is exactly why B3–B5 and B7 went unnoticed), nor for the four oscillator stubs.
- Existing tests construct problems and cross-compare forms / energy error, but never evaluate
  DAE constraint functions (hides B1, B10) nor invoke plot recipes (hides B11).

---

## 6. Repository hygiene

- Untracked backup/scratch files: `src/*-prev.jl` (8), `test/*-prev.jl`, and untracked
  `examples/`, `prototyping/`, `notebooks/`, `references/`, `harmonic-oscillator-plots/`. None are
  `include`d; recommend git-ignoring rather than deleting.
- `src/lotka_volterra_3d_lagrangian.jl` (B14) is untracked, unregistered, and broken.

---

## 7. Verified correct (no change needed)

Pendulum (all forms), double pendulum (incl. mass-matrix inversion & Hamiltonian), rigid-body
Euler equations + Casimir, Lorenz equations, Toda lattice (periodic BC), ABC flow, Kubo
SDE/PSDE/SPSDE, Lotka-Volterra 2D core dynamics, Lotka-Volterra 4D and 4D-Lagrangian (mutually
consistent to ~`eps()`), harmonic-oscillator DELE derivatives `D1Ld`/`D2Ld` (midpoint & trapezoidal).
