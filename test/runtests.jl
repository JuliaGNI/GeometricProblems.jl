using SafeTestsets

@safetestset "Bump initial condition: test derivative.                                        " begin include("bump_initial_condition_test_derivatives.jl") end
@safetestset "ABC Flow                                                                        " begin include("abc_flow_tests.jl") end
@safetestset "Double Pendulum                                                                 " begin include("double_pendulum_tests.jl") end
@safetestset "Harmonic Oscillator                                                             " begin include("harmonic_oscillator_tests.jl") end
# @safetestset "Kubo Oscillator                                                                 " begin include("kubo_oscillator_tests.jl") end
@safetestset "Lorenz Attractor                                                                " begin include("lorenz_attractor_tests.jl") end
@safetestset "Lotka-Volterra 2D                                                               " begin include("lotka_volterra_2d_tests.jl") end
@safetestset "Lotka-Volterra 2D with singular Lagrangian                                      " begin include("lotka_volterra_2d_singular_tests.jl") end
@safetestset "Lotka-Volterra 2D with symmetric Lagrangian                                     " begin include("lotka_volterra_2d_symmetric_tests.jl") end
@safetestset "Lotka-Volterra 2D with symmetric Lagrangian with gauge terms                    " begin include("lotka_volterra_2d_gauge_tests.jl") end
@safetestset "Lotka-Volterra 3D                                                               " begin include("lotka_volterra_3d_tests.jl") end
@safetestset "Lotka-Volterra 4D                                                               " begin include("lotka_volterra_4d_tests.jl") end
@safetestset "Lotka-Volterra 4D (Lagrangian)                                                  " begin include("lotka_volterra_4d_lagrangian_tests.jl") end
@safetestset "Point Vortices                                                                  " begin include("point_vortices_tests.jl") end
@safetestset "Point Vortices (linear)                                                         " begin include("point_vortices_linear_tests.jl") end
@safetestset "Rigid Body                                                                      " begin include("rigid_body_test.jl") end