"""
GeometricProblems.jl is a collection of ODEs and DAEs with interesting geometric structures.
"""
module GeometricProblems

    include("diagnostics.jl")
    include("plot_recipes.jl")

    include("double_pendulum.jl")
    include("duffing_oscillator.jl")
    include("harmonic_oscillator.jl")
    include("kubo_oscillator.jl")
    include("lennard_jones_oscillator.jl")
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
    include("toda_lattice.jl")

end
