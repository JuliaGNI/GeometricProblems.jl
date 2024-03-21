@doc raw"""
# ABC Flow

```math
\begin{aligned}
    \dot{x} = A\sin(z) + C\cos(y) \\
    \dot{y} = B\sin(x) + A\cos(z) \\ 
    \dot{z} = C\sin(y) + B\cos(x)
\end{aligned}
```
"""
module ABCFlow

    using GeometricEquations 
    using GeometricSolutions
    using Parameters 

    export odeproblem, q₀, default_parameters, tspan, tstep

    const tspan = (0.0, 100.0)
    const tstep = 0.1

    const default_parameters = (
        A = 0.5,
        B = 1.,
        C = 1.
    )

    const q₀ = [0., 0., 0.]

    function abc_flow_v(v, t, q, params)
        @unpack A, B, C = params
        v[1] = A * sin(q[3]) + C * cos(q[2])
        v[2] = B * sin(q[1]) + A * cos(q[3])
        v[3] = C * sin(q[2]) + B * cos(q[1])
        
        nothing
    end

    function odeproblem(q₀ = (q = q₀); tspan = tspan, tstep = tstep, parameters = default_parameters)
        ODEProblem(abc_flow_v, tspan, tstep, q₀; parameters = parameters)
    end

end
