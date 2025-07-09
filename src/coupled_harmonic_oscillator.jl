@doc raw"""
    CoupledHarmonicOscillator

The `CoupledHarmonicOscillator` module provides functions `hodeproblem` and `lodeproblem`
each returning a Hamiltonian or Lagrangian problem, respectively, to be solved
in the GeometricIntegrators.jl ecosystem.
The actual code is generated with EulerLagrange.jl.

The coupled harmonic oscillator is a collection of two point masses that are connected to a fixed wall with spring constants ``k_1`` and ``k_2`` and are furthermore coupled nonlinearly resulting in the Hamiltonian: 

```math 
H(q_1, q_2, p_1, p_2) = \frac{q_1^2}{2m_1} + \frac{q_2^2}{2m_2} + k_1\frac{q_1^2}{2} + k_2\frac{q_2^2}{2} +  k\sigma(q_1)\frac{(q_2 - q_1)^2}{2},
```
where ``\sigma(x) = 1 / (1 + e^{-x})`` is the sigmoid activation function. 


System parameters:
* `k₁`: spring constant of mass 1
* `k₂`: spring constant of mass 2
* `m₁`: mass 1
* `m₂`: mass 2
* `k`: coupling strength between the two masses. 
"""
module CoupledHarmonicOscillator

    using EulerLagrange
    using LinearAlgebra
    using Parameters
    using GeometricEquations: HODEEnsemble

    export hamiltonian, lagrangian
    export hodeproblem, lodeproblem
    export hodeensemble
    export hamiltonian_system, lagrangian_system

    const timespan = (0.0, 100.0)
    const timestep = 0.4

    const default_parameters = (
        m₁ = 2.,
        m₂ = 1.,
        k₁ = 1.5,
        k₂ = 0.3,
        k = 1.0
    )

    const q₀ = [1., 0.]
    const p₀ = [2., 0.]

    function σ(x::T) where {T<:Real}
        T(1) / (T(1) + exp(-x))
    end

    function hamiltonian(t, q, p, parameters)
        @unpack k₁, k₂, m₁, m₂, k = parameters

        p[1] ^ 2 / (2 * m₁) + p[2] ^ 2 / (2 * m₂) + k₁ * q[1] ^ 2 / 2 + k₂ * q[2] ^ 2 / 2  + k * σ(q[1]) * (q[2] - q[1]) ^2 / 2
    end


    function lagrangian(t, q, q̇, parameters)
        @unpack k₁, k₂, m₁, m₂, k = parameters

        q̇[1] ^ 2 / (2 * m₁) + q̇[2] ^ 2 / (2 * m₂) - k₁ * q[1] ^ 2 / 2 - k₂ * q[2] ^ 2 / 2  - k * σ(q[1]) * (q[2] - q[1]) ^2 / 2
    end

    function v̄(v, t, q, p, parameters)
        v[1] = p[1] / parameters.m₁
        v[2] = p[2] / parameters.m₂
        nothing
    end

    function hamiltonian_system(parameters::NamedTuple)
        t, q, p = hamiltonian_variables(2)
        sparams = symbolize(parameters)
        HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
    end

    function lagrangian_system(parameters::NamedTuple)
        t, x, v = lagrangian_variables(2)
        sparams = symbolize(parameters)
        LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
    end

    _parameters(p::NamedTuple) = p
    _parameters(p::AbstractVector) = p[begin]

    
    """
        Hamiltonian problem for coupled oscillator

    Constructor with default arguments:
    ```
    hodeproblem(
        q₀ = $(q₀),
        p₀ = $(p₀);
        timespan = $(timespan),
        timestep = $(timestep),
        parameters = $(default_parameters)
    )
    ```
    """
    function hodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        HODEProblem(hamiltonian_system(parameters), timespan, timestep, q₀, p₀; parameters = parameters)
    end

    function hodeensemble(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        eqs = functions(hamiltonian_system(_parameters(parameters)))
        HODEEnsemble(eqs.v, eqs.f, eqs.H, timespan, timestep, q₀, p₀; parameters = parameters)
    end

    """
        Lagrangian problem for the coupled oscillator

    Constructor with default arguments:
    ```
    lodeproblem(
        q₀ = $(q₀),
        p₀ = $(p₀);
        timespan = $(timespan),
        timestep = $(timestep),
        parameters = $(default_parameters)
    )
    ```
    """
    function lodeproblem(q₀ = q₀, p₀ = p₀; timespan = timespan, timestep = timestep, parameters = default_parameters)
        LODEProblem(lagrangian_system(parameters), timespan, timestep, q₀, p₀; v̄ = v̄, parameters = parameters)
    end

end
