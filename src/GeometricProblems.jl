"""
GeometricProblems.jl is a collection of ODEs and DAEs with interesting geometric structures.
"""
module GeometricProblems

    include("bump_initial_condition.jl")

    include("diagnostics.jl")
    include("plot_recipes.jl")

    include("abc_flow.jl")
    include("double_pendulum.jl")
    include("duffing_oscillator.jl")
    include("harmonic_oscillator.jl")
    include("kubo_oscillator.jl")
    include("lennard_jones_oscillator.jl")
    include("linear_wave.jl")
    include("lorenz_attractor.jl")
    include("lotka_volterra_2d.jl")
    include("lotka_volterra_2d_gauge.jl")
    include("lotka_volterra_2d_singular.jl")
    include("lotka_volterra_2d_symmetric.jl")
    include("lotka_volterra_2d_plots.jl")
    include("lotka_volterra_3d.jl")
    include("lotka_volterra_3d_plots.jl")
    include("lotka_volterra_4d.jl")
    include("lotka_volterra_4d_lagrangian.jl")
    include("lotka_volterra_4d_plots.jl")
    include("massless_charged_particle.jl")
    include("massless_charged_particle_plots.jl")
    include("coupled_harmonic_oscillator.jl")
    include("mathews_lakshmanan_oscillator.jl")
    include("morse_oscillator.jl")
    include("pendulum.jl")
    include("point_vortices.jl")
    include("point_vortices_linear.jl")
    include("rigid_body.jl")
    include("three_body_problem.jl")
    include("toda_lattice.jl")
    include("outer_solar_system.jl")
    include("perturbed_pendulum.jl")
    include("henon_heiles_potential.jl")
end
