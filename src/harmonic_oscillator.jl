@doc raw"""
# Harmonic Oscillator

"""
module HarmonicOscillator

    using GeometricEquations
    using GeometricSolutions
    using Parameters

    export harmonic_oscillator_ode, harmonic_oscillator_iode, harmonic_oscillator_pode, harmonic_oscillator_hode, harmonic_oscillator_sode,
           harmonic_oscillator_dae, harmonic_oscillator_idae, harmonic_oscillator_pdae, harmonic_oscillator_hdae

    export hamiltonian, compute_energy_error


    const t₀ = 0.0
    const Δt = 0.1
    const nt = 10
    const tspan = (t₀, Δt*nt)

    const k = 0.5
    const ω = √k

    const default_parameters = (k=k, ω=ω)
    
    ϑ₁(t,q) = q[2]
    ϑ₂(t,q) = zero(eltype(q))

    function ϑ(q)
        p = zero(q)
        p[1] = ϑ₁(0,q)
        p[2] = ϑ₂(0,q)
        return p
    end

    function hamiltonian(t, q, params)
        @unpack k = params
        q[2]^2 / 2 + k * q[1]^2 / 2
    end

    function hamiltonian(p, t, q, params)
        @unpack k = params
        p[1]^2 / 2 + k * q[1]^2 / 2
    end

    function lagrangian(v, t, q, params)
        ϑ₁(t,q) * v[1] + ϑ₂(t,q) * v[2] - hamiltonian(t, q, params)
    end


    const q₀ = [0.5, 0.0]
    const z₀ = [0.5, 0.0, 0.5]
    const p₀ = ϑ(q₀)

    const A = sqrt(q₀[2]^2 / k + q₀[1]^2)
    const ϕ = asin(q₀[1] / A)

    const reference_solution_q = A * sin(ω * Δt * nt + ϕ)
    const reference_solution_p = ω * A * cos(ω * Δt * nt + ϕ)

    const reference_solution = [reference_solution_q, reference_solution_p]
    

    function oscillator_ode_v(v, t, q, params)
        @unpack k = params
        v[1] = q[2]
        v[2] = -k * q[1]
        nothing
    end

    function harmonic_oscillator_ode(q₀=q₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert size(x₀,1) == 2
        # ODE(oscillator_ode_v, q₀; invariants=(h=hamiltonian,), parameters=params)
        ODEProblem(oscillator_ode_v, tspan, tstep, q₀; invariants = (h=hamiltonian,), parameters = parameters)
    end


    function oscillator_pode_v(v, t, q, p, params)
        v[1] = p[1]
        nothing
    end

    function oscillator_pode_f(f, t, q, p, params)
        @unpack k = params
        f[1] = -k*q[1]
        nothing
    end

    function harmonic_oscillator_pode(q₀=[q₀[1]], p₀=[p₀[1]]; parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert length(q₀) == length(p₀)
        # @assert all([length(q) == length(p) == 1 for (q,p) in zip(q₀,p₀)])
        # @assert size(q₀,1) == size(p₀,1) == 1
        PODEProblem(oscillator_pode_v, oscillator_pode_f, tspan, tstep, q₀, p₀; invariants = (h=hamiltonian,), parameters = parameters)
    end


    function harmonic_oscillator_hode(q₀=[q₀[1]], p₀=[p₀[1]]; parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert length(q₀) == length(p₀)
        # @assert all([length(q) == length(p) == 1 for (q,p) in zip(q₀,p₀)])
        # @assert size(q₀,1) == size(p₀,1) == 1
        HODEProblem(oscillator_pode_v, oscillator_pode_f, hamiltonian, tspan, tstep, q₀, p₀; parameters = parameters)
    end


    function oscillator_sode_v_1(v, t, q, params)
        v[1] = q[2]
        v[2] = 0
        nothing
    end

    function oscillator_sode_v_2(v, t, q, params)
        @unpack k = params
        v[1] = 0
        v[2] = -k * q[1]
        nothing
    end

    function oscillator_sode_q_1(q₁, t₁, q₀, t₀, params)
        q₁[1] = q₀[1] + (t₁ - t₀) * q₀[2]
        q₁[2] = q₀[2]
        nothing
    end

    function oscillator_sode_q_2(q₁, t₁, q₀, t₀, params)
        @unpack k = params
        q₁[1] = q₀[1]
        q₁[2] = q₀[2] - (t₁ - t₀) * k * q₀[1]
        nothing
    end

    function harmonic_oscillator_sode(q₀=q₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        SODEProblem((oscillator_sode_v_1, oscillator_sode_v_2),
                    (oscillator_sode_q_1, oscillator_sode_q_2),
                    tspan, tstep, q₀; parameters = parameters)
    end


    function oscillator_iode_ϑ(p, t, q, params)
        p[1] = q[2]
        p[2] = 0
        nothing
    end

    function oscillator_iode_ϑ(p, t, q, v, params)
        oscillator_iode_ϑ(p, t, q, params)
    end

    function oscillator_iode_f(f, t, q, v, params)
        @unpack k = params
        f[1] = -k*q[1]
        f[2] = v[1] - q[2]
        nothing
    end

    function oscillator_iode_g(g, t, q, λ, params)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function oscillator_iode_v(v, t, q, params)
        @unpack k = params
        v[1] = q[2]
        v[2] = -k*q[1]
        nothing
    end

    function harmonic_oscillator_iode(q₀=q₀, p₀=ϑ(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert size(q₀,1) == size(p₀,1) == 2
        IODEProblem(oscillator_iode_ϑ, oscillator_iode_f,
             oscillator_iode_g, tspan, tstep, q₀, p₀;
             invariants = (h=hamiltonian,), parameters = parameters,
             v̄ = oscillator_iode_v)
    end


    function oscillator_dae_v(v, t, q, params)
        @unpack k = params
        v[1] = q[2]
        v[2] = -k * q[1]
        v[3] = q[2] - k * q[1]
        nothing
    end

    function oscillator_dae_u(u, t, q, λ, params)
        u[1] = -λ[1]
        u[2] = -λ[1]
        u[3] = +λ[1]
    end

    function oscillator_dae_ϕ(ϕ, t, q, params)
        ϕ[1] = q[3] - q[1] - q[2]
    end

    function harmonic_oscillator_dae(q₀=z₀, λ₀=[zero(eltype(q₀))]; parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert size(z₀,1) == 3
        # @assert size(λ₀,1) == 1
        # @assert all([length(z) == 3 for z in z₀])
        # @assert all([length(λ) == 1 for λ in λ₀])
        # DAE(oscillator_dae_v, oscillator_dae_u, oscillator_dae_ϕ,
        #     z₀, λ₀; invariants=(h=hamiltonian,), parameters=params, v̄=oscillator_ode_v)
        DAEProblem(oscillator_ode_v, oscillator_dae_u, oscillator_dae_ϕ, tspan, tstep, q₀, λ₀;
                    v̄ = oscillator_ode_v, invariants = (h=hamiltonian,), parameters = parameters)
    end


    function oscillator_idae_u(u, t, q, p, λ, params)
        u[1] = λ[1]
        u[2] = λ[2]
        nothing
    end

    function oscillator_idae_g(g, t, q, p, λ, params)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function oscillator_hdae_ū(u, t, q, p, λ, params)
        u[1] = λ[1]
        u[2] = λ[2]
        nothing
    end

    function oscillator_hdae_ḡ(g, t, q, p, λ, params)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function oscillator_idae_ϕ(ϕ, t, q, p, params)
        ϕ[1] = p[1] - q[2]
        ϕ[2] = p[2]
        nothing
    end

    function oscillator_hdae_ψ(ψ, t, q, p, v, f, params)
        ψ[1] = f[1] - v[2]
        ψ[2] = f[2]
        nothing
    end

    function harmonic_oscillator_idae(q₀=q₀, p₀=ϑ(q₀), λ₀=zero(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert size(q₀,1) == size(p₀,1) == size(λ₀,1) == 2
        IDAEProblem(oscillator_iode_ϑ, oscillator_iode_f,
                    oscillator_idae_u, oscillator_idae_g, oscillator_idae_ϕ,
                    tspan, tstep, q₀, p₀, λ₀; v̄ = oscillator_iode_v, invariants = (h=hamiltonian,), parameters = parameters)
    end

    function oscillator_pdae_v(v, t, q, p, params)
        @unpack k = params
        v[1] = q[2]
        v[2] = -k*q[1]
        nothing
    end

    function oscillator_pdae_f(f, t, q, p, params)
        @unpack k = params
        f[1] = -k*q[1]
        f[2] = p[1] - q[2]
        nothing
    end

    function harmonic_oscillator_pdae(q₀=q₀, p₀=ϑ(q₀), λ₀=zero(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert size(q₀,1) == size(p₀,1) == 2
        PDAEProblem(oscillator_pdae_v, oscillator_pdae_f,
                    oscillator_idae_u,  oscillator_idae_g, oscillator_idae_ϕ,
                    tspan, tstep, q₀, p₀, λ₀; invariants=(h=hamiltonian,), parameters = parameters)
    end

    function harmonic_oscillator_hdae(q₀=q₀, p₀=ϑ(q₀), λ₀=zero(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        # @assert size(q₀,1) == size(p₀,1) == 2
        HDAEProblem(oscillator_pdae_v, oscillator_pdae_f, 
                    oscillator_idae_u, oscillator_idae_g, oscillator_idae_ϕ,
                    oscillator_idae_ū, oscillator_idae_ḡ, oscillator_idae_ψ,
                    hamiltonian, tspan, tstep, q₀, p₀, λ₀; parameters = parameters)
    end


    function compute_energy_error(t, q::DataSeries{T}, params) where {T}
        h = DataSeries(T, q.nt)
        e = DataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i], params)
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

end
