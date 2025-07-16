module HenonHeilesPotential

    using GeometricEquations
    using EulerLagrange
    using LinearAlgebra
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem

    const default_parameters= (
        λ = 1.0,
    )
    
    const p₀ = [0.3, 0.3]
    const q₀ = [0.2, 0.2]
    const x₀ = vcat(q₀, p₀)

    const timestep = 0.1
    const timespan = (0.0, 10.0)

    function hamiltonian(t,q,p,params)
        @unpack λ = params
        0.5 * (p[1]^2 + p[2]^2) + 0.5 * (q[1]^2 + q[2]^2) + λ * (q[1]^2 * q[2] - q[2]^3 / 3) 
    end

    function lagrangian(t,q,v,params)
        @unpack λ = params
        0.5 * (v[1]^2 + v[2]^2) - 0.5 * (q[1]^2 + q[2]^2) - λ * (q[1]^2 * q[2] - q[2]^3 / 3) 
    end        

    θ̇₁(t, q, p, params) = p[1]
    θ̇₂(t, q, p, params) = p[2]

    function θ̇(v, t, q, p, params)
        v[1] = θ̇₁(t, q, p, params)
        v[2] = θ̇₂(t, q, p, params)
        nothing
    end


    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, q, p = hamiltonian_variables(2)
        sparams = symbolize(parameters)
        ham_sys = HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
        HODEProblem(ham_sys, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        t, x, v = lagrangian_variables(2)
        sparams = symbolize(parameters)
        lag_sys = LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
        LODEProblem(lag_sys, timespan, timestep, q₀, p₀; v̄ = θ̇, parameters = parameters)
    end

end