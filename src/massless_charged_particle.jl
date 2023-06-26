@doc raw"""
# Massless charged particle in 2D

The Lagrangian is given by
```math
L(x, \dot{x}) = A(x) \cdot \dot{x} - \phi (x) ,
```
with magnetic vector potential
```math
A(x) = \frac{A_0}{2} \big( 1 + x_1^2 + x_2^2 \big) \begin{pmatrix}
- x_2 \\
+ x_1 \\
\end{pmatrix} ,
```
electrostatic potential
```math
\phi(x) =  E_0 \, \big( \cos (x_1) + \sin(x_2) \big) ,
```
and magnetic and electric fields
```math
\begin{aligned}
B(x) &= \nabla \times A(x) = A_0 \, (1 + 2 x_1^2 + 2 x_2^2) , \\
E(x) &= - \nabla \phi(x) = E_0 \, \big( \sin x_1, \, - \cos x_2 \big)^T .
\end{aligned}
```

The Hamiltonian form of the equations of motion reads
```math
\dot{x} = \frac{1}{B(x)} \begin{pmatrix}
\hphantom{-} 0 & + 1 \\
- 1 & \hphantom{+} 0 \\
\end{pmatrix} \nabla \phi (x) .
```

"""
module MasslessChargedParticle

    using GeometricEquations
    
    using ..Diagnostics
    
    import ..Diagnostics: compute_invariant_error, compute_momentum_error

    export ϑ, A, B, ϕ, E, hamiltonian
    export massless_charged_particle_ode, massless_charged_particle_iode,
           massless_charged_particle_idae, massless_charged_particle_idae_spark
    export compute_energy_error, compute_momentum_error


    # default simulation parameters
    const Δt = 0.2
    const nt = 5000
    const tspan = (0.0, Δt*nt)

    # default initial conditions and parameters
    q₀ = [1.0, 1.0]
    
    default_parameters = (A₀ = 1.0, E₀ = 1.0)

    # components of the vector potential
    A₁(q, params) = - params[:A₀] * q[2] * (1 + q[1]^2 + q[2]^2) / 2
    A₂(q, params) = + params[:A₀] * q[1] * (1 + q[1]^2 + q[2]^2) / 2

    A(q, params) = [A₁(q, params), A₂(q, params)]

    # z-componend of the magnetic field
    B(q, params) = params[:A₀] * (1 + 2 * q[1]^2 + 2 * q[2]^2)

    # electrostatic potential
    ϕ(q, params) = params[:E₀] * (cos(q[1]) + sin(q[2]))
    # ϕ(q, params) = E₀ * (q[1]^2 + q[2]^2)

    # components of the electric field
    E₁(q, params) = + params[:E₀] * sin(q[1])
    E₂(q, params) = - params[:E₀] * cos(q[2])
    # E₁(q, params) = + params[:E₀] * q[1]
    # E₂(q, params) = + params[:E₀] * q[2]

    E(q, params) = [E₁(q, params), E₂(q, params)]

    # components of the velocity
    v₁(t, q, params) = - E₂(q, params) / B(q, params)
    v₂(t, q, params) = + E₁(q, params) / B(q, params)

    # components of the one-form (symplectic potential)
    ϑ₁(t, q, params) = A₁(q, params)
    ϑ₂(t, q, params) = A₂(q, params)

    ϑ(t, q, params) = [ϑ₁(t, q, params), ϑ₂(t, q, params)]

    function ϑ(t::Number, q::AbstractVector, params::NamedTuple, k::Int)
        if k == 1
            ϑ₁(t, q, params)
        elseif k == 2
            ϑ₂(t, q, params)
        else
            throw(BoundsError(ϑ,k))
        end
    end

    # derivatives of the one-form components
    dϑ₁dx₁(t, q, params) = - params[:A₀] * q[1] * q[2]
    dϑ₁dx₂(t, q, params) = - params[:A₀] * (1 + q[1]^2 + 3 * q[2]^2) / 2
    dϑ₂dx₁(t, q, params) = + params[:A₀] * (1 + 3 * q[1]^2 + q[2]^2) / 2
    dϑ₂dx₂(t, q, params) = + params[:A₀] * q[1] * q[2]

    # components of the force
    f₁(v, t, q, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₂dx₁(t, q, params) * v[2]
    f₂(v, t, q, params) = dϑ₁dx₂(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

    g₁(v, t, q, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₁dx₂(t, q, params) * v[2]
    g₂(v, t, q, params) = dϑ₂dx₁(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

    # Hamiltonian (total energy)
    hamiltonian(t, q, params) = ϕ(q, params)

    # components of the gradient of the Hamiltonian
    dHd₁(t, q, params) = - E₁(q, params)
    dHd₂(t, q, params) = - E₂(q, params)


    function massless_charged_particle_dH(dH, t, q, params)
        dH[1] = dHd₁(t, q, params)
        dH[2] = dHd₂(t, q, params)
        nothing
    end 

    function massless_charged_particle_v(v, t, q, params)
        v[1] = v₁(t, q, params)
        v[2] = v₂(t, q, params)
        nothing
    end

    function massless_charged_particle_v(v, t, q, p, params)
        massless_charged_particle_v(v, t, q, params)
    end

    function massless_charged_particle_ϑ(Θ, t, q, params)
        Θ[1] = ϑ₁(t, q, params)
        Θ[2] = ϑ₂(t, q, params)
        nothing 
    end

    massless_charged_particle_ϑ(Θ, t, q, v, params) = massless_charged_particle_ϑ(Θ, t, q, params)

    function massless_charged_particle_f(f, t, q, v, params)
        f[1] = f₁(v, t, q, params) - dHd₁(t, q, params)
        f[2] = f₂(v, t, q, params) - dHd₂(t, q, params)
        nothing
    end

    function massless_charged_particle_f̄(f, t, q, v, params)
        f[1] = - dHd₁(t, q, params)
        f[2] = - dHd₂(t, q, params)
        nothing
    end

    function massless_charged_particle_g(g, t, q, v, params)
        g[1] = f₁(v, t, q, params)
        g[2] = f₂(v, t, q, params)
        nothing
    end

    massless_charged_particle_g(g, t, q, p, v, params) = massless_charged_particle_g(g, t, q, v, params)

    function massless_charged_particle_u(u, t, q, v, params)
        u .= v
        nothing
    end

    massless_charged_particle_u(u, t, q, p, v, params) = massless_charged_particle_u(u, t, q, v, params)

    function massless_charged_particle_ϕ(ϕ, t, q, p, params)
        ϕ[1] = p[1] - ϑ₁(t,q,params)
        ϕ[2] = p[2] - ϑ₂(t,q,params)
        nothing
    end

    function massless_charged_particle_ψ(ψ, t, q, p, v, f, params)
        ψ[1] = f[1] - g₁(t,q,v,params)
        ψ[2] = f[2] - g₂(t,q,v,params)
        nothing
    end



    "Creates an ODE object for the massless charged particle in 2D."
    function massless_charged_particle_ode(q₀=q₀; params=default_parameters)
        ODE(massless_charged_particle_v, q₀; invariants=(h=hamiltonian,), parameters=params)
    end

    "Creates an implicit ODE object for the massless charged particle in 2D."
    function massless_charged_particle_iode(q₀=q₀; params=default_parameters)
        IODE(massless_charged_particle_ϑ, massless_charged_particle_f,
                massless_charged_particle_g, q₀, ϑ(0., q₀, params);
                v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
                invariants=(h=hamiltonian,), parameters=params)
    end

    "Creates an implicit DAE object for the massless charged particle in 2D."
    function massless_charged_particle_idae(q₀=q₀; params=default_parameters)
        IDAE(massless_charged_particle_ϑ, massless_charged_particle_f,
                massless_charged_particle_u, massless_charged_particle_g,
                massless_charged_particle_ϕ, q₀, ϑ(0., q₀, params), zero(q₀);
                v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
                invariants=(h=hamiltonian,), parameters=params)
    end

    "Creates an implicit DAE object for the massless charged particle in 2D."
    function massless_charged_particle_idae_spark(q₀=q₀; params=default_parameters)
        IDAE(massless_charged_particle_ϑ, massless_charged_particle_f̄,
                massless_charged_particle_u, massless_charged_particle_g,
                massless_charged_particle_ϕ, q₀, ϑ(0., q₀, params), zero(q₀);
                v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
                invariants=(h=hamiltonian,), parameters=params)
    end


    compute_energy_error(t,q,params) = compute_invariant_error(t,q, (t,q) -> hamiltonian(t,q,params))
    compute_momentum_error(t,q,p,params::NamedTuple) = compute_momentum_error(t, q, p, (t,q,k) -> ϑ(t,q,params,k))

end
