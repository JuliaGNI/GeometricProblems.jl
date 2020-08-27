
using SafeTestsets

@safetestset "Harmonic Oscillator                                                             " begin include("harmonic_oscillator_tests.jl") end
@safetestset "Lorenz Attractor                                                                " begin include("lorenz_attractor_tests.jl") end
@safetestset "Lotka-Volterra 2D                                                               " begin include("lotka_volterra_2d_tests.jl") end
@safetestset "Lotka-Volterra 2D with singular Lagrangian                                      " begin include("lotka_volterra_2d_singular_tests.jl") end
@safetestset "Lotka-Volterra 2D with symmetric Lagrangian                                     " begin include("lotka_volterra_2d_symmetric_tests.jl") end
@safetestset "Lotka-Volterra 2D with symmetric Lagrangian with gauge terms                    " begin include("lotka_volterra_2d_gauge_tests.jl") end
@safetestset "Lotka-Volterra 3D                                                               " begin include("lotka_volterra_3d_tests.jl") end
@safetestset "Point Vortices                                                                  " begin include("point_vortices_tests.jl") end
