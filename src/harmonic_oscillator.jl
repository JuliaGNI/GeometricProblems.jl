@doc raw"""
# Harmonic Oscillator

"""
module HarmonicOscillator

    using GeometricEquations
    using GeometricSolutions
    using Parameters

    export odeproblem, podeproblem, hodeproblem, iodeproblem, lodeproblem, sodeproblem,
           daeproblem, pdaeproblem, hdaeproblem, idaeproblem, ldaeproblem,
           degenerate_iodeproblem, degenerate_lodeproblem

    export odeensemble, podeensemble, hodeensemble

    export hamiltonian, lagrangian

    export compute_energy_error, exact_solution


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

    function ω!(ω, t, q, params)
        ω[1,1] = 0
        ω[1,2] = -1
        ω[2,1] = +1
        ω[2,2] = 0
        nothing
    end

    function hamiltonian(t, q, params)
        @unpack k = params
        q[2]^2 / 2 + k * q[1]^2 / 2
    end

    function hamiltonian(t, q, p, params)
        @unpack k = params
        p[1]^2 / 2 + k * q[1]^2 / 2
    end

    function lagrangian(t, q, v, params)
        @unpack k = params
        v[1]^2 / 2 - k * q[1]^2 / 2
    end

    function degenerate_lagrangian(t, q, v, params)
        ϑ₁(t,q) * v[1] + ϑ₂(t,q) * v[2] - hamiltonian(t, q, params)
    end

    function discrete_lagrangian(h,q̄,q,params)
        h*lagrangian(0, (q̄ + q)/2, (q - q̄)/h, params)
    end

    A(q, p, params) = q * sqrt(1 + p^2 / q^2 / params.k)
    ϕ(q, p, params) = atan(p / q / params.ω)

    exact_solution_q(t, q₀, p₀, t₀, params) = A(q₀, p₀, params) * cos(params.ω * (t-t₀) - ϕ(q₀, p₀, params))
    exact_solution_p(t, q₀, p₀, t₀, params) = - params.ω * A(q₀, p₀, params) * sin(params.ω * (t-t₀) - ϕ(q₀, p₀, params))

    exact_solution_q(t, q₀::AbstractVector, p₀::AbstractVector, t₀, params) = exact_solution_q(t, q₀[1], p₀[1], t₀, params)
    exact_solution_p(t, q₀::AbstractVector, p₀::AbstractVector, t₀, params) = exact_solution_p(t, q₀[1], p₀[1], t₀, params)
    
    exact_solution_q(t, x₀::AbstractVector, t₀, params) = exact_solution_q(t, x₀[1], x₀[2], t₀, params)
    exact_solution_p(t, x₀::AbstractVector, t₀, params) = exact_solution_p(t, x₀[1], x₀[2], t₀, params)
    exact_solution(t, x₀::AbstractVector, t₀, params) = [exact_solution_q(t, x₀, t₀, params), exact_solution_p(t, x₀, t₀, params)]

    
    const q₀ = [0.5]
    const p₀ = [0.0]
    const x₀ = vcat(q₀, p₀)

    const xmin = [-2., -2.]
    const xmax = [+2., +2.]
    const nsamples = [10, 10]

    const reference_solution_q = exact_solution_q(Δt * nt, q₀[1], p₀[1], t₀, default_parameters)
    const reference_solution_p = exact_solution_p(Δt * nt, q₀[1], p₀[1], t₀, default_parameters)

    const reference_solution = [reference_solution_q, reference_solution_p]
    

    function _ode_samples(qmin, qmax, nsamples)
        qs = [range(qmin[i], qmax[i]; length = nsamples[i]) for i in eachindex(qmin, qmax, nsamples)]

        samples = vec(collect.(collect(Base.Iterators.product(qs...))))

        (q = samples,)
    end

    function _pode_samples(qmin, qmax, pmin, pmax, qsamples, psamples)
        qs = [range(qmin[i], qmax[i]; length = qsamples[i]) for i in eachindex(qmin, qmax, qsamples)]
        ps = [range(pmin[i], pmax[i]; length = psamples[i]) for i in eachindex(pmin, pmax, psamples)]

        qsamples = vec(collect.(collect(Base.Iterators.product(qs...))))
        psamples = vec(collect.(collect(Base.Iterators.product(ps...))))
        zsamples = Base.Iterators.product(qsamples, psamples)

        return (
            q = vec([zs[1] for zs in zsamples]),
            p = vec([zs[2] for zs in zsamples]),
        )
    end


    function oscillator_ode_v(v, t, x, params)
        @unpack k = params
        v[1] = x[2]
        v[2] = -k * x[1]
        nothing
    end

    function odeproblem(x₀ = x₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(x₀) == 2
        ODEProblem(oscillator_ode_v, tspan, tstep, x₀; invariants = (h=hamiltonian,), parameters = parameters)
    end

    function odeensemble(qmin = xmin, qmax = xmax, nsamples = nsamples; parameters = default_parameters, tspan = tspan, tstep = Δt)
        samples = _ode_samples(qmin, qmax, nsamples)
        ODEEnsemble(oscillator_ode_v, tspan, tstep, samples...; invariants = (h=hamiltonian,), parameters = parameters)
    end

    function exact_solution!(sol::GeometricSolution, prob::ODEProblem)
        for n in eachtimestep(sol)
            sol.q[n] .= exact_solution(sol.t[n], sol.q[0], sol.t[0], parameters(prob))
        end
        return sol
    end

    function exact_solution(prob::ODEProblem)
        exact_solution!(GeometricSolution(prob), prob)
    end


    function oscillator_pode_v(v, t, q, p, params)
        v[1] = p[1]
        nothing
    end

    function oscillator_pode_f(f, t, q, p, params)
        @unpack k = params
        f[1] = -k * q[1]
        nothing
    end

    function podeproblem(q₀ = q₀, p₀ = p₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 1
        PODEProblem(oscillator_pode_v, oscillator_pode_f, tspan, tstep, q₀, p₀; invariants = (h=hamiltonian,), parameters = parameters)
    end

    function hodeproblem(q₀ = q₀, p₀ = p₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 1
        HODEProblem(oscillator_pode_v, oscillator_pode_f, hamiltonian, tspan, tstep, q₀, p₀; parameters = parameters)
    end

    function podeensemble(qmin = [xmin[1]], qmax = [xmax[1]], pmin = [xmin[2]], pmax = [xmax[2]], qsamples = [nsamples[1]], psamples = [nsamples[2]]; parameters = default_parameters, tspan = tspan, tstep = Δt)
        samples = _pode_samples(qmin, qmax, pmin, pmax, qsamples, psamples)     
        PODEEnsemble(oscillator_pode_v, oscillator_pode_f, tspan, tstep, samples...; invariants = (h=hamiltonian,), parameters = parameters)
    end

    function hodeensemble(qmin = [xmin[1]], qmax = [xmax[1]], pmin = [xmin[2]], pmax = [xmax[2]], qsamples = [nsamples[1]], psamples = [nsamples[2]]; parameters = default_parameters, tspan = tspan, tstep = Δt)
        samples = _pode_samples(qmin, qmax, pmin, pmax, qsamples, psamples)     
        HODEEnsemble(oscillator_pode_v, oscillator_pode_f, hamiltonian, tspan, tstep, samples...; parameters = parameters)
    end

    function exact_solution!(sol::GeometricSolution, prob::Union{PODEProblem,HODEProblem})
        for n in eachtimestep(sol)
            sol.q[n] = [exact_solution_q(sol.t[n], sol.q[0], sol.p[0], sol.t[0], parameters(prob))]
            sol.p[n] = [exact_solution_p(sol.t[n], sol.q[0], sol.p[0], sol.t[0], parameters(prob))]
        end
        return sol
    end

    function exact_solution(prob::Union{PODEProblem,HODEProblem})
        exact_solution!(GeometricSolution(prob), prob)
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

    function sodeproblem(x₀ = x₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        SODEProblem((oscillator_sode_v_1, oscillator_sode_v_2),
                    (oscillator_sode_q_1, oscillator_sode_q_2),
                    tspan, tstep, x₀; v̄ = oscillator_ode_v, parameters = parameters)
    end

    function oscillator_iode_ϑ(p, t, q, v, params)
        p[1] = v[1]
        nothing
    end

    function oscillator_iode_f(f, t, q, v, params)
        @unpack k = params
        f[1] = -k * q[1]
        nothing
    end

    function oscillator_iode_g(g, t, q, v, λ, params)
        g[1] = λ[1]
        nothing
    end

    function oscillator_iode_v(v, t, q, p, params)
        v[1] = p[1]
        nothing
    end

    function iodeproblem(q₀=q₀, p₀=p₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 1
        IODEProblem(oscillator_iode_ϑ, oscillator_iode_f,
             oscillator_iode_g, tspan, tstep, q₀, p₀;
             invariants = (h=hamiltonian,), parameters = parameters,
             v̄ = oscillator_iode_v)
    end

    function lodeproblem(q₀=q₀, p₀=p₀; parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 1
        LODEProblem(oscillator_iode_ϑ, oscillator_iode_f,
             oscillator_iode_g, ω!, lagrangian,
             tspan, tstep, q₀, p₀;
             invariants = (h=hamiltonian,), parameters = parameters,
             v̄ = oscillator_iode_v)
    end


    function degenerate_oscillator_iode_ϑ(p, t, q, v, params)
        p[1] = q[2]
        p[2] = 0
        nothing
    end

    function degenerate_oscillator_iode_f(f, t, q, v, params)
        @unpack k = params
        f[1] = -k * q[1]
        f[2] = v[1] - q[2]
        nothing
    end

    function degenerate_oscillator_iode_g(g, t, q, v, λ, params)
        g[1] = 0
        g[2] = λ[1]
        nothing
    end

    function degenerate_oscillator_iode_v(v, t, q, p, params)
        @unpack k = params
        v[1] = q[2]
        v[2] = -k * q[1]
        nothing
    end

    function degenerate_iodeproblem(q₀ = x₀, p₀ = ϑ(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 2
        IODEProblem(degenerate_oscillator_iode_ϑ, degenerate_oscillator_iode_f,
        degenerate_oscillator_iode_g, tspan, tstep, q₀, p₀;
             invariants = (h=hamiltonian,), parameters = parameters,
             v̄ = degenerate_oscillator_iode_v)
    end

    function degenerate_lodeproblem(q₀ = x₀, p₀ = ϑ(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 2
        LODEProblem(degenerate_oscillator_iode_ϑ, degenerate_oscillator_iode_f,
             degenerate_oscillator_iode_g, ω!, lagrangian,
             tspan, tstep, q₀, p₀;
             invariants = (h=hamiltonian,), parameters = parameters,
             v̄ = degenerate_oscillator_iode_v)
    end


    function oscillator_dae_u(u, t, x, λ, params)
        @unpack k = params
        u[1] = k * x[1] * λ[1]
        u[2] = x[2] * λ[1]
    end

    function oscillator_dae_ϕ(ϕ, t, x, params)
        ϕ[1] = hamiltonian(t, x, params) - hamiltonian(t₀, x₀, params)
    end

    function daeproblem(x₀=x₀, λ₀=[zero(eltype(x₀))]; parameters = default_parameters, tspan = tspan, tstep = Δt)
        DAEProblem(oscillator_ode_v, oscillator_dae_u, oscillator_dae_ϕ, tspan, tstep, x₀, λ₀;
                    v̄ = oscillator_ode_v, invariants = (h=hamiltonian,), parameters = parameters)
    end


    function oscillator_pdae_v(v, t, q, p, params)
        @unpack k = params
        v[1] = p[1]
        nothing
    end

    function oscillator_pdae_f(f, t, q, p, params)
        @unpack k = params
        f[1] = -k * q[1]
        nothing
    end

    function oscillator_pdae_u(u, t, q, p, λ, params)
        @unpack k = params
        u[1] = k * q[1] * λ[1]
        nothing
    end

    function oscillator_pdae_g(g, t, q, p, λ, params)
        g[1] = p[1] * λ[1]
        nothing
    end

    function oscillator_pdae_ū(u, t, q, p, λ, params)
        @unpack k = params
        u[1] = k * q[1] * λ[1]
        nothing
    end

    function oscillator_pdae_ḡ(g, t, q, p, λ, params)
        g[1] = p[1] * λ[1]
        nothing
    end

    function oscillator_pdae_ϕ(ϕ, t, q, p, params)
        ϕ[1] = hamitlonian(t, q, p, params)
        nothing
    end

    function oscillator_pdae_ψ(ψ, t, q, p, q̇, ṗ, params)
        @unpack k = params
        ψ[1] = p[1] * ṗ[1] + k * q[1] * q̇[1]
        nothing
    end

    function pdaeproblem(q₀ = q₀, p₀ = p₀, λ₀ = zero(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 1
        PDAEProblem(oscillator_pdae_v, oscillator_pdae_f,
                    oscillator_pdae_u, oscillator_pdae_g, oscillator_pdae_ϕ,
                    tspan, tstep, q₀, p₀, λ₀; invariants=(h=hamiltonian,), parameters = parameters)
    end

    function hdaeproblem(q₀ = q₀, p₀ = p₀, λ₀ = zero(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == 1
        HDAEProblem(oscillator_pdae_v, oscillator_pdae_f, 
                    oscillator_pdae_u, oscillator_pdae_g, oscillator_pdae_ϕ,
                    oscillator_pdae_ū, oscillator_pdae_ḡ, oscillator_pdae_ψ,
                    hamiltonian, tspan, tstep, q₀, p₀, λ₀; parameters = parameters)
    end


    oscillator_idae_u(u, t, q, v, p, λ, params) = oscillator_pdae_u(u, t, q, p, λ, params)
    oscillator_idae_g(g, t, q, v, p, λ, params) = oscillator_pdae_g(g, t, q, p, λ, params)
    oscillator_idae_ū(u, t, q, v, p, λ, params) = oscillator_pdae_ū(u, t, q, p, λ, params)
    oscillator_idae_ḡ(g, t, q, v, p, λ, params) = oscillator_pdae_ḡ(g, t, q, p, λ, params)
    oscillator_idae_ϕ(ϕ, t, q, v, p, params) = oscillator_pdae_ϕ(ϕ, t, q, p, params)
    oscillator_idae_ψ(ψ, t, q, v, p, q̇, ṗ, params) = oscillator_pdae_ψ(ψ, t, q, p, q̇, ṗ, params)

    function idaeproblem(q₀ = q₀, p₀ = p₀, λ₀ = zero(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == length(λ₀) == 1
        IDAEProblem(oscillator_iode_ϑ, oscillator_iode_f,
                    oscillator_idae_u, oscillator_idae_g, oscillator_idae_ϕ,
                    tspan, tstep, q₀, p₀, λ₀; v̄ = oscillator_iode_v, invariants = (h=hamiltonian,), parameters = parameters)
    end

    function ldaeproblem(q₀ = q₀, p₀ = p₀, λ₀ = zero(q₀); parameters = default_parameters, tspan = tspan, tstep = Δt)
        @assert length(q₀) == length(p₀) == length(λ₀) == 1
        LDAEProblem(oscillator_iode_ϑ, oscillator_iode_f,
                    oscillator_idae_u, oscillator_idae_g, oscillator_idae_ϕ, ω!, lagrangian,
                    tspan, tstep, q₀, p₀, λ₀; v̄ = oscillator_iode_v, invariants = (h=hamiltonian,), parameters = parameters)
    end


    function exact_solution(probs::Union{ODEEnsemble,PODEEnsemble,HODEEnsemble})
        sols = EnsembleSolution(probs)
        for (sol, prob) in zip(sols, probs)
            exact_solution!(sol, prob)
        end
        return sols
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
