
# GeometricProblems.jl

*Collection of example problems with interesting geometric structure for
[GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl).*

[![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliagni.github.io/GeometricProblems.jl/stable)
[![Latest Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://juliagni.github.io/GeometricProblems.jl/latest)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE.md)
[![Build Status](https://github.com/JuliaGNI/GeometricProblems.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaGNI/GeometricProblems.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/JuliaGNI/GeometricProblems.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaGNI/GeometricProblems.jl)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.3740036.svg)](https://doi.org/10.5281/zenodo.3740036)


#### Example Problems

- [x] ABC Flow,
- [x] Exponential Growth,
- [ ] Fermi-Pasta-Ulam Problem,
- [ ] Hénon-Heiles System,
- [ ] Kepler Problem,
- [x] Lorenz Attractor in 3D,
- [x] Lotka-Volterra in 2D,
- [x] Lotka-Volterra in 3D,
- [x] Lotka-Volterra in 4D,
- [x] Massless Charged Particle,
- [x] Harmonic Oscillator,
- [x] Coupled Harmonic Oscillator,
- [ ] Nonlinear Oscillators,
    - [ ] Duffing Oscillator,
    - [ ] Lennard-Jones Oscillator,
    - [ ] Mathews-Lakshmanan Oscillator,
    - [ ] Morse Oscillator,
- [x] Pendulum,
    - [x] Mathematical Pendulum,
    - [x] Double Pendulum,
- [x] Planar Point Vortices,
- [ ] Rigid Body,
- [ ] Chaplygin Sleigh,
- [ ] Inner Solar System,
- [ ] Outer Solar System,
- [ ] Heavy Top,
- [x] Toda lattice.


See [ChargedParticleDynamics.jl](https://github.com/JuliaPlasma/ChargedParticleDynamics.jl) for

- Charged Particle Motion in various electromagnetic Fields,
- Pauli Particle Dynamics in various electromagnetic Fields,
- Guiding Center Dynamics in various magnetic fields,
- Gyrokinetic Dynamics in various magnetic fields.


See [GeometricExamples.jl](https://github.com/JuliaGNI/GeometricExamples.jl) for
example scripts that run these problems with the integrators implemented in
[GeometricIntegrators.jl](https://github.com/JuliaGNI/GeometricIntegrators.jl).


## Development

### Git hooks

Two hooks live in `.githooks`. They are **not active in a fresh clone** — `core.hooksPath` is local
configuration and does not travel with a push — so enable them once per clone:

```sh
git config core.hooksPath .githooks
```

**`pre-commit`** acts on **staged `.jl` files only**, and exits immediately when a commit stages
none, so a documentation- or workflow-only commit is not slowed down by it:

- **JuliaFormatter `--check`**, honouring this repository's own `.JuliaFormatter.toml` — **blocks**
  the commit. Formatting is mechanical and always fixable.
- **`fatou lint`**, when `fatou` is installed — **advisory only**, and deliberately so: its
  `unused-import` rule does not follow `include`, so it flags the load-bearing imports of every
  module file.
- **`using <Package>`**, which catches a syntax error or a broken `include` — **blocks**.

**`pre-push`** runs the full test suite with `--check-bounds=auto`, but **only when pushing to
`main` or `master`**; a topic branch is left to CI. It prints nothing for **10–30 minutes**, which
looks exactly like a network hang and is not one. If you do interrupt it, check for an orphaned
Julia process that the killed hook left behind.

Either hook can be bypassed for a single command with `--no-verify`, for a change you know it does
not apply to:

```sh
git commit --no-verify
git push --no-verify
```

The hooks are generated from one shared copy and are byte-identical across the related
repositories, so edit them there rather than here — a local edit is silently undone by the next
install.
