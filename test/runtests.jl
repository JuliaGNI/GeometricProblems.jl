using SafeTestsets

const GROUPS = isempty(ARGS) ? ["core", "slow"] : ARGS

if "core" in GROUPS
    @safetestset "Aqua" include("quality/aqua.jl")
    @safetestset "Unicode normalization (NFC)" include("integration/unicode_normalization.jl")
    @safetestset "Default timespan/timestep constants" include("integration/default_timespan_timestep.jl")
    @safetestset "Bump initial condition: test derivative." include("bump_initial_condition.jl")
    @safetestset "ABC Flow" include("abc_flow.jl")
    @safetestset "Coupled Harmonic Oscillator" include("coupled_harmonic_oscillator.jl")
    @safetestset "Double Pendulum" include("double_pendulum.jl")
    @safetestset "Harmonic Oscillator" include("harmonic_oscillator.jl")
    @safetestset "Kubo Oscillator" include("kubo_oscillator.jl")
    @safetestset "Nonlinear Oscillators" include("nonlinear_oscillators.jl")
    @safetestset "Linear Wave" include("linear_wave.jl")
    @safetestset "Massless Charged Particle" include("massless_charged_particle.jl")
    @safetestset "Massless Charged Particle (singular)" include("massless_charged_particle_singular.jl")
    @safetestset "Three-Body Problem" include("three_body_problem.jl")
    @safetestset "Lorenz Attractor" include("lorenz_attractor.jl")
    @safetestset "Lotka-Volterra 2D" include("lotka_volterra_2d.jl")
    @safetestset "Lotka-Volterra 2D with singular Lagrangian" include("lotka_volterra_2d_singular.jl")
    @safetestset "Lotka-Volterra 2D with symmetric Lagrangian" include("lotka_volterra_2d_symmetric.jl")
    @safetestset "Lotka-Volterra 2D with symmetric Lagrangian with gauge terms" include("lotka_volterra_2d_gauge.jl")
    @safetestset "Lotka-Volterra 3D" include("lotka_volterra_3d.jl")
    @safetestset "Outer Solar System" include("outer_solar_system.jl")
    @safetestset "Pendulum" include("pendulum.jl")
    @safetestset "Perturbed Pendulum" include("perturbed_pendulum.jl")
    @safetestset "Point Vortices" include("point_vortices.jl")
    @safetestset "Point Vortices (linear)" include("point_vortices_linear.jl")
    @safetestset "Rigid Body" include("rigid_body.jl")
    @safetestset "Toda Lattice" include("toda_lattice.jl")
    @safetestset "Henon Heiles Potential" include("henon_heiles_potential.jl")
    @safetestset "LODE/LDAE ω and l wiring" include("integration/lode_wiring.jl")
    @safetestset "Poincaré invariants extensions" include("integration/poincare_invariants.jl")
end
if "slow" in GROUPS
    @safetestset "Lotka-Volterra 4D" include("lotka_volterra_4d.jl")
    @safetestset "Lotka-Volterra 4D (Lagrangian)" include("lotka_volterra_4d_lagrangian.jl")
    @safetestset "Euler-Lagrange ensembles" include("integration/eulerlagrange_ensembles.jl")
    @safetestset "Plotting extensions" include("integration/plots.jl")
end
