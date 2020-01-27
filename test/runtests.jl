
using SafeTestsets

@safetestset "Lotka-Volterra 2D Tests                                                         " begin include("lotka_volterra_2d_tests.jl") end
@safetestset "Point Vortices                                                                  " begin include("point_vortices_tests.jl") end
