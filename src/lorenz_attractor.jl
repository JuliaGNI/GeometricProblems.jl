@doc raw"""
# Lorenz Attractor

"""
module LorenzAttractor

    using GeometricEquations

    export lorenz_attractor_ode, plot_lorenz_attractor


    const Δt = 0.01
    const nt = 1000
    const tspan = (0.0, Δt*nt)

    const q₀ = [1., 1., 1.]

    const default_parameters = (σ = 10., ρ = 28., β = 8/3)
    const reference_solution = [-4.902687541134471, -3.743872921802973, 24.690858102790042]

    function lorenz_attractor_v(v, t, x, params)
        σ, ρ, β = params
        v[1] = σ * (x[2] - x[1])
        v[2] = x[1] * (ρ - x[3]) - x[2]
        v[3] = x[1] * x[2] - β * x[3]
        nothing
    end


    function lorenz_attractor_ode(q₀=q₀; tspan = tspan, tstep = Δt, parameters = default_parameters)
        ODEProblem(lorenz_attractor_v, tspan, tstep, q₀; parameters = parameters)
    end

end
