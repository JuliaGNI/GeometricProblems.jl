@doc raw"""
# Lorenz Attractor

The Lorenz system is a three-dimensional, non-Hamiltonian ODE exhibiting chaotic dynamics
(deterministic, aperiodic flow on a strange attractor):
```math
\begin{aligned}
\dot{x} &= \sigma (y - x) , \\
\dot{y} &= x (\rho - z) - y , \\
\dot{z} &= x y - \beta z .
\end{aligned}
```

System parameters:
* `σ`: Prandtl number
* `ρ`: (reduced) Rayleigh number
* `β`: geometric factor

The default parameters ``(\sigma, \rho, \beta) = (10, 28, 8/3)`` are the classic values of
Lorenz (1963), for which the system is chaotic. Only the plain `ODEProblem` form is provided:
the flow is dissipative (volume-contracting) and carries no Hamiltonian or symplectic structure.
"""
module LorenzAttractor

    using GeometricEquations

    export lorenz_attractor_ode


    const Δt = 0.01
    const nt = 1000
    const DEFAULT_TIMESPAN = (0.0, Δt*nt)

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


    function lorenz_attractor_ode(q₀=q₀; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        ODEProblem(lorenz_attractor_v, timespan, timestep, q₀; parameters = parameters)
    end

end
