using SafeTestsets

@safetestset "Unicode normalization (NFC)                                                     " begin
    include("unicode_normalization_tests.jl")
end
@safetestset "Default timespan/timestep constants                                             " begin
    include("default_timespan_timestep_tests.jl")
end
@safetestset "Bump initial condition: test derivative.                                        " begin
    include("bump_initial_condition_test_derivatives.jl")
end
@safetestset "ABC Flow                                                                        " begin
    include("abc_flow_tests.jl")
end
@safetestset "Coupled Harmonic Oscillator                                                     " begin
    include("coupled_harmonic_oscillator_tests.jl")
end
@safetestset "Double Pendulum                                                                 " begin
    include("double_pendulum_tests.jl")
end
@safetestset "Harmonic Oscillator                                                             " begin
    include("harmonic_oscillator_tests.jl")
end
@safetestset "Kubo Oscillator                                                                 " begin
    include("kubo_oscillator_tests.jl")
end
@safetestset "Nonlinear Oscillators                                                           " begin
    include("nonlinear_oscillators_tests.jl")
end
@safetestset "Linear Wave                                                                     " begin
    include("linear_wave_tests.jl")
end
@safetestset "Massless Charged Particle                                                       " begin
    include("massless_charged_particle_tests.jl")
end
@safetestset "Massless Charged Particle (singular)                                            " begin
    include("massless_charged_particle_singular_tests.jl")
end
@safetestset "Three-Body Problem                                                              " begin
    include("three_body_tests.jl")
end
@safetestset "Lorenz Attractor                                                                " begin
    include("lorenz_attractor_tests.jl")
end
@safetestset "Lotka-Volterra 2D                                                               " begin
    include("lotka_volterra_2d_tests.jl")
end
@safetestset "Lotka-Volterra 2D with singular Lagrangian                                      " begin
    include("lotka_volterra_2d_singular_tests.jl")
end
@safetestset "Lotka-Volterra 2D with symmetric Lagrangian                                     " begin
    include("lotka_volterra_2d_symmetric_tests.jl")
end
@safetestset "Lotka-Volterra 2D with symmetric Lagrangian with gauge terms                    " begin
    include("lotka_volterra_2d_gauge_tests.jl")
end
@safetestset "Lotka-Volterra 3D                                                               " begin
    include("lotka_volterra_3d_tests.jl")
end
@safetestset "Lotka-Volterra 4D                                                               " begin
    include("lotka_volterra_4d_tests.jl")
end
@safetestset "Lotka-Volterra 4D (Lagrangian)                                                  " begin
    include("lotka_volterra_4d_lagrangian_tests.jl")
end
@safetestset "Outer Solar System                                                              " begin
    include("outer_solar_system_tests.jl")
end
@safetestset "Pendulum                                                                        " begin
    include("pendulum_tests.jl")
end
@safetestset "Perturbed Pendulum                                                              " begin
    include("perturbed_pendulum_tests.jl")
end
@safetestset "Point Vortices                                                                  " begin
    include("point_vortices_tests.jl")
end
@safetestset "Point Vortices (linear)                                                         " begin
    include("point_vortices_linear_tests.jl")
end
@safetestset "Rigid Body                                                                      " begin
    include("rigid_body_test.jl")
end
@safetestset "Toda Lattice                                                                    " begin
    include("toda_lattice_tests.jl")
end
@safetestset "Henon Heiles Potential                                                          " begin
    include("henon_heiles_potential_tests.jl")
end
@safetestset "Euler-Lagrange ensembles                                                        " begin
    include("eulerlagrange_ensembles_tests.jl")
end
@safetestset "LODE/LDAE ω and l wiring                                                        " begin
    include("lode_wiring_tests.jl")
end
@safetestset "Plotting extensions                                                             " begin
    include("plots_tests.jl")
end
