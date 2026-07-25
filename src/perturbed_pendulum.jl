@doc raw"""
    PerturbedPendulum

A mathematical pendulum with a non-separable perturbation. The Hamiltonian is
```math
H(q, p) = \frac{p^2}{2} - \omega^2 \cos(q) - q \, p \, A(\epsilon, \phi) ,
```
where ``q`` is the angle, ``p`` is the conjugate momentum, and
```math
A(\epsilon, \phi) = 0.3 \, \epsilon \sin(2\phi) + 0.7 \, \epsilon \sin(3\phi) .
```
The corresponding Lagrangian obtained via Legendre transform is
```math
L(q, \dot{q}) = \frac{\dot{q}^2}{2} + \omega^2 \cos(q)
              + q \, \dot{q} \, A(\epsilon, \phi)
              + \frac{1}{2} q^2 A(\epsilon, \phi)^2 .
```

System parameters:
* `ω`: natural frequency of the unperturbed pendulum
* `ϵ`: perturbation amplitude
* `ϕ`: perturbation phase
"""
module PerturbedPendulum

    using EulerLagrange
    using Parameters

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hamiltonian_system, lagrangian_system
    export default_parameters

    const DEFAULT_TIMESTEP = 0.1
    const DEFAULT_TIMESPAN = (0.0, 10.0)

    const q₀ = [0.5]
    const p₀ = [0.0]

    default_parameters(::Type{T} = Float64) where {T} = (ω = T(0.5), ϵ = T(0.5), ϕ = T(π/3))

    # The perturbation coefficient A(ϵ, ϕ) = 0.3 ϵ sin(2ϕ) + 0.7 ϵ sin(3ϕ). It is *derived* from ϵ
    # and ϕ, so it must not be stored alongside them in the parameter tuple: a user overriding ϵ or
    # ϕ would otherwise silently keep a stale A. `symbolize` handles the expression, so the
    # symbolic Hamiltonian/Lagrangian pick up the dependence on ϵ and ϕ correctly.
    perturbation(params) = 0.3 * params.ϵ * sin(2 * params.ϕ) + 0.7 * params.ϵ * sin(3 * params.ϕ)

    # velocity from momentum: q̇ = ∂H/∂p = p - q*A
    function θ̇(t, q, p, params)
        p[1] - q[1] * perturbation(params)
    end

    function θ̇(v, t, q, p, params)
        v[1] = θ̇(t, q, p, params)
        nothing
    end

    function hamiltonian(t, q, p, params)
        @unpack ω = params
        A = perturbation(params)
        p[1]^2 / 2 - ω^2 * cos(q[1]) - q[1] * p[1] * A
    end

    function lagrangian(t, q, v, params)
        @unpack ω = params
        A = perturbation(params)
        v[1]^2 / 2 + ω^2 * cos(q[1]) + q[1] * v[1] * A + q[1]^2 * A^2 / 2
    end

    function hamiltonian_system(parameters::NamedTuple)
        t, q, p = hamiltonian_variables(1)
        sparams  = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(1)
        sparams  = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
    end

    """
        hodeproblem(q₀, p₀; timespan, timestep, parameters)

    Hamiltonian formulation of the perturbed pendulum.

    Constructor with default arguments:
    ```
    hodeproblem(
        q₀ = $(q₀),
        p₀ = $(p₀);
        timespan   = $(DEFAULT_TIMESPAN),
        timestep   = $(DEFAULT_TIMESTEP),
        parameters = $(default_parameters())
    )
    ```
    """
    function hodeproblem(q₀ = q₀, p₀ = p₀;
                         timespan   = DEFAULT_TIMESPAN,
                         timestep   = DEFAULT_TIMESTEP,
                         parameters = default_parameters())
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀;
                    parameters = parameters)
    end

    """
        lodeproblem(q₀, p₀; timespan, timestep, parameters)

    Lagrangian formulation of the perturbed pendulum.

    Constructor with default arguments:
    ```
    lodeproblem(
        q₀ = $(q₀),
        p₀ = $(p₀);
        timespan   = $(DEFAULT_TIMESPAN),
        timestep   = $(DEFAULT_TIMESTEP),
        parameters = $(default_parameters())
    )
    ```
    """
    function lodeproblem(q₀ = q₀, p₀ = p₀;
                         timespan   = DEFAULT_TIMESPAN,
                         timestep   = DEFAULT_TIMESTEP,
                         parameters = default_parameters())
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀;
                    v̄ = θ̇, parameters = parameters)
    end

end
