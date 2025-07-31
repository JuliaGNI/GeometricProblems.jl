@doc raw"""
undampled Duffing oscillator
The `DuffingOscillator` module provides a Hamiltonian problem to be solved
in the GeometricIntegrators.jl ecosystem.
The Duffing oscillator is a nonlinear second-order differential equation
describing the motion of a damped and driven oscillator with a nonlinear
restoring force. It is often used to model systems with a cubic nonlinearity.
The equation of motion is given by:
```
    dot(x) = y
    dot(y) =  αx - βx^3 - δ * y 
```
where `x` is the displacement, `y` is the velocity, and `δ` is the damping coefficient.
here we only consider the undamped case, i.e., `δ = 0`, which can be described by Hamiltonian and Lagrangian systems.
"""
module DuffingOscillator

    using GeometricEquations
    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const default_parameters= (
        δ = 0.0,  # damping coefficient
        m = 2.0,  # mass
        α = -1.0,  # linear stiffness
        β = 1.0,  # cubic stiffness
    )
    
    const p₀ = [1.2]
    const q₀ = [1.2]
    const x₀ = vcat(q₀, p₀)

    const timestep = 0.1
    const timespan = (0.0, 10.0)

    function hamiltonian(t, q, p, params)
        @unpack δ, m, α, β = params
        kinetic = p[1]^2 / (2 * m)
        potential = 0.5 * α * q[1]^2 + 0.25 * β * q[1]^4
        return kinetic + potential
    end

    function lagrangian(t, q, v, params)
        @unpack δ, m, α, β = params
        kinetic = 0.5 * m * v[1]^2
        potential = 0.5 * α * q[1]^2 + 0.25 * β * q[1]^4
        return kinetic - potential
    end

    function θ̇(t, q, p, params)
        @unpack δ, m, α, β = params
        p[1]/m
    end

    function θ̇(v, t, q, p, params)
        v[1] = θ̇(t, q, p, params)
        nothing
    end


    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(1)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, x, v = lagrangian_variables(1)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, timespan, timestep, q₀, p₀; v̄ = θ̇, parameters = parameters)
    end

end
