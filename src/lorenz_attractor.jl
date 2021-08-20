@doc raw"""
# Lorenz Attractor

"""
module LorenzAttractor

    using GeometricEquations

    export lorenz_attractor_ode, plot_lorenz_attractor


    Δt = 0.01
    nt = 1000

    const q₀ = [1., 1., 1.]
    const default_params = (σ = 10., ρ = 28., β = 8/3)

    function lorenz_attractor_v(t, x, v, params=default_params)
        σ, ρ, β = params
        v[1] = σ * (x[2] - x[1])
        v[2] = x[1] * (ρ - x[3]) - x[2]
        v[3] = x[1] * x[2] - β * x[3]
        nothing
    end


    function lorenz_attractor_ode(q₀=q₀)
        ODE(lorenz_attractor_v, q₀)
    end

end
