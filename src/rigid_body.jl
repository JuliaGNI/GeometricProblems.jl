@doc raw"""
# Rigid body

```math
\begin{aligned}
    \dot{x} = Ayz \\
    \dot{y} = Bxz \\ 
    \dot{z} = Cxy
\end{aligned},
```

where ``A = (I_2 - I_3)/(I_2I_3)``, ``B = (I_3 - I_1)/(I_3I_1)`` and ``C = (I_1 - I_2)/(I_1I_2)``; ``I_{\cdot}`` are the *principal components of inertia*. 

The initial condition and the default parameters are taken from [bajars2023locally](@cite).
"""
module RigidBody

    using GeometricEquations 
    using GeometricSolutions
    using Parameters 

    export odeproblem, q₀, default_parameters, tspan, tstep

    const tspan = (0.0, 100.0)
    const tstep = 0.1

    const default_parameters = (
        I₁ = 2.,
        I₂ = 1.,
        I₃ = 2. / 3.
    )

    const q₀ = (q = [cos(1.1), 0., sin(1.1)], )

    function rigid_body_v(v, t, q, params)
        @unpack I₁, I₂, I₃ = params
        A = (I₂ - I₃) / (I₂ * I₃)
        B = (I₃ - I₁) / (I₃ * I₁)
        C = (I₁ - I₂) / (I₁ * I₂)
        v[1] = A * q[2] * q[3]
        v[2] = B * q[1] * q[3]
        v[3] = C * q[1] * q[2]
        
        nothing
    end

    function odeproblem(q₀ = q₀; tspan = tspan, tstep = tstep, parameters = default_parameters)
        ODEProblem(rigid_body_v, tspan, tstep, q₀; parameters = parameters)
    end

end