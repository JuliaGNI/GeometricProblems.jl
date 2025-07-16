module PerturbedPendulum

    using GeometricEquations
    using EulerLagrange
    using LinearAlgebra
    using Parameters


    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const default_parameters= (
        ω = 0.5,
        ϵ = 0.5,
        m = 1.0,
        ϕ = pi/3,
    )
    
    const p₀ = [0.0]
    const q₀ = [0.5]
    const x₀ = vcat(q₀, p₀)

    const timestep = 0.1
    const timespan = (0.0, 10.0)

    function hamiltonian(t,q,p,params)
        @unpack ω,ϵ,ϕ = params
        (p[1]^2)/2 - ω^2 * cos(q[1]) - q[1]* p[1] * (0.3 * ϵ * sin(2*ϕ) + 0.7 * ϵ * sin(3*ϕ))
    end

    function lagrangian(t,q,v,params)
        @unpack ω,ϵ,ϕ = params
        (v[1]^2)/2 + ω^2 * cos(q[1]) + q[1]* v[1] * (0.3 * ϵ * sin(2*ϕ) + 0.7 * ϵ * sin(3*ϕ)) + 0.5 * q[1]^2 * (0.3 * ϵ * sin(2*ϕ) + 0.7 * ϵ * sin(3*ϕ))^2  
    end

    function θ̇(t, q, p, params)
        @unpack ω,ϵ,ϕ = params
        p[1] - q[1]* (0.3 * ϵ * sin(2*ϕ) + 0.7 * ϵ * sin(3*ϕ))
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