__precompile__()

module GeometricProblems

    using GeometricIntegrators

    export ExponentialGrowth, LotkaVolterra2d, Pendulum, PointVortices

    include("exponential_growth.jl")
    include("lotka_volterra_2d.jl")
    include("pendulum.jl")
    include("point_vortices.jl")

end
