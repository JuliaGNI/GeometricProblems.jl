__precompile__()

module GeometricProblems

    export ExponentialGrowth
    export LorenzAttractor
    export LotkaVolterra2d
    export Pendulum
    export PointVortices


    include("exponential_growth.jl")
    include("lorenz_attractor.jl")
    include("lotka_volterra_2d.jl")
    include("pendulum.jl")
    include("point_vortices.jl")

end
