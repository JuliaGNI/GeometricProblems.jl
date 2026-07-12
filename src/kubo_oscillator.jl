@doc raw"""
# Kubo Oscillator

The Kubo oscillator is a unit-frequency harmonic oscillator, $\dot{q}_1 = q_2$, $\dot{q}_2 = -q_1$,
driven by multiplicative (Stratonovich) noise. It is a standard test problem for stochastic
geometric integrators. The module provides it as a stochastic differential equation
(`kubo_oscillator_sde_*`), a partitioned SDE (`kubo_oscillator_psde_*`), a split partitioned SDE
(`kubo_oscillator_spsde_*`), and the underlying deterministic ODE (`kubo_oscillator_ode`). The
`_1`/`_2` variants build single problems and the `_3` variants build ensembles of initial
conditions. The noise process is represented by [`KuboNoise`](@ref).

System parameter: `ν` — the noise intensity.
"""
module KuboOscillator

    using GeometricEquations
    using Parameters

    import GeometricBase: AbstractStochasticProcess

    export kubo_oscillator_sde_1, kubo_oscillator_psde_1, kubo_oscillator_spsde_1
    export kubo_oscillator_sde_2, kubo_oscillator_psde_2, kubo_oscillator_spsde_2
    export kubo_oscillator_sde_3, kubo_oscillator_psde_3, kubo_oscillator_spsde_3
    export kubo_oscillator_ode
    export KuboNoise

    """
        KuboNoise <: AbstractStochasticProcess

    Marker for the (one-dimensional, unit-intensity Wiener) noise process driving the Kubo
    oscillator. The GeometricEquations SDE/PSDE/SPSDE API expects a noise object of type
    `AbstractStochasticProcess`; the concrete realisation of the increments is supplied by the
    stochastic integrator.
    """
    struct KuboNoise <: AbstractStochasticProcess end

    q_init_A=[0.5, 0.0]
    q_init_B=[[ 0.5, 0.0],
              [ 0.0, 0.5],
              [-0.5, 0.0]]

    const noise_intensity = 0.1

    const Δt = 0.01
    const nt = 10
    const DEFAULT_TIMESPAN = (0.0, Δt*nt)

    const default_parameters = (ν = noise_intensity,)


    function kubo_oscillator_sde_v(v, t, q, params)
        v[1]=  q[2]
        v[2]= -q[1]
    end


    function kubo_oscillator_sde_B(B::AbstractVector, t, q, params)
        @unpack ν = params
        B[1] = +ν * q[2]
        B[2] = -ν * q[1]
    end

    function kubo_oscillator_sde_B(B::AbstractMatrix, t, q, params)
        @unpack ν = params
        for j in axes(B, 2)
            B[1,j] = +ν * q[2]
            B[2,j] = -ν * q[1]
        end
    end


    function kubo_oscillator_sde_1(q₀=q_init_A; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_A - interpreted as one random initial conditions with one sample path
        # 1-dimensional noise
        SDEProblem(kubo_oscillator_sde_v, kubo_oscillator_sde_B, KuboNoise(), timespan, timestep, q₀; parameters = parameters)
    end

    function kubo_oscillator_sde_2(q₀=q_init_A; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_A - single deterministic initial condition
        # Generating 3 sample paths
        # 1-dimensional noise
        SDEProblem(kubo_oscillator_sde_v, kubo_oscillator_sde_B, KuboNoise(), timespan, timestep, q₀; parameters = parameters)
    end

    function kubo_oscillator_sde_3(q₀=q_init_B; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_B holds several initial conditions -> ensemble problem
        equ = SDE(kubo_oscillator_sde_v, kubo_oscillator_sde_B, KuboNoise();
                  parameters = GeometricEquations.parameter_types(parameters))
        ics = [(q = StateVariable(x),) for x in q₀]
        GeometricEquations.EnsembleProblem(equ, timespan, timestep, ics, parameters)
    end


    # ODE

    function kubo_oscillator_ode(q₀=q_init_A; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        ODEProblem(kubo_oscillator_sde_v, timespan, timestep, q₀; parameters = parameters)
    end


    # PSDE

    q_init_C=[0.5]
    p_init_C=[0.0]

    q_init_D=[[0.5], [0.0], [-0.5]]
    p_init_D=[[0.0], [0.5], [ 0.0]]


    function kubo_oscillator_psde_v(v, t, q, p, params)
        v[1] =  p[1]
    end

    function kubo_oscillator_psde_f(f, t, q, p, params)
        f[1] = -q[1]
    end

    function kubo_oscillator_psde_B(B, t, q, p, params)
        @unpack ν = params
        B[1,1] = +ν * p[1]
    end

    function kubo_oscillator_psde_G(G, t, q, p, params)
        @unpack ν = params
        G[1,1] = -ν * q[1]
    end


    function kubo_oscillator_psde_1(q₀=q_init_C, p₀=p_init_C; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_C - interpreted as a single random initial condition with one sample path
        # 1-dimensional noise
        PSDEProblem(kubo_oscillator_psde_v, kubo_oscillator_psde_f,
                    kubo_oscillator_psde_B, kubo_oscillator_psde_G, KuboNoise(),
                    timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function kubo_oscillator_psde_2(q₀=q_init_C, p₀=p_init_C; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_C - single deterministic initial condition
        # Generating 3 sample paths
        # 1-dimensional noise
        PSDEProblem(kubo_oscillator_psde_v, kubo_oscillator_psde_f,
                    kubo_oscillator_psde_B, kubo_oscillator_psde_G, KuboNoise(),
                    timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function kubo_oscillator_psde_3(q₀=q_init_D, p₀=p_init_D; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_D / p_init_D hold several initial conditions -> ensemble problem
        equ = PSDE(kubo_oscillator_psde_v, kubo_oscillator_psde_f,
                   kubo_oscillator_psde_B, kubo_oscillator_psde_G, KuboNoise();
                   parameters = GeometricEquations.parameter_types(parameters))
        ics = [(q = StateVariable(x), p = StateVariable(y)) for (x, y) in zip(q₀, p₀)]
        GeometricEquations.EnsembleProblem(equ, timespan, timestep, ics, parameters)
    end


    # SPSDE

    function kubo_oscillator_spsde_v(v, t, q, p, params)
        v[1] =  p[1]
    end

    function kubo_oscillator_spsde_f1(f, t, q, p, params)
        f[1] = -q[1]
    end

    function kubo_oscillator_spsde_f2(f, t, q, p, params)
        f[1] = 0
    end

    function kubo_oscillator_spsde_B(B, t, q, p, params)
        @unpack ν = params
        B[1,1] = +ν * p[1]
    end

    function kubo_oscillator_spsde_G1(G, t, q, p, params)
        @unpack ν = params
        G[1,1] = -ν * q[1]
    end

    function kubo_oscillator_spsde_G2(G, t, q, p, params)
        G[1,1] = 0
    end


    function kubo_oscillator_spsde_1(q₀ = q_init_C, p₀ = p_init_C; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_C - interpreted as a single random initial condition with one sample path
        # 1-dimensional noise
        SPSDEProblem(kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2,
                     kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, KuboNoise(),
                     timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function kubo_oscillator_spsde_2(q₀ = q_init_C, p₀ = p_init_C; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_C - single deterministic initial condition
        # Generating 3 sample paths
        # 1-dimensional noise
        SPSDEProblem(kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2,
                     kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, KuboNoise(),
                     timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function kubo_oscillator_spsde_3(q₀ = q_init_D, p₀ = p_init_D; timespan = DEFAULT_TIMESPAN, timestep = Δt, parameters = default_parameters)
        # q_init_D / p_init_D hold several initial conditions -> ensemble problem
        equ = SPSDE(kubo_oscillator_spsde_v, kubo_oscillator_spsde_f1, kubo_oscillator_spsde_f2,
                    kubo_oscillator_spsde_B, kubo_oscillator_spsde_G1, kubo_oscillator_spsde_G2, KuboNoise();
                    parameters = GeometricEquations.parameter_types(parameters))
        ics = [(q = StateVariable(x), p = StateVariable(y)) for (x, y) in zip(q₀, p₀)]
        GeometricEquations.EnsembleProblem(equ, timespan, timestep, ics, parameters)
    end

end
