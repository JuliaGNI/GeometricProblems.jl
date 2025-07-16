@doc raw"""
# Lotka-Volterra Model in 3D

The Lotka–Volterra model in 3D is an example of a Hamiltonian system with degenerate Poisson structure.

The equations read
```math
\begin{aligned}
\dot{q}_{1} &= q_{1} (           -  a_{2} q_{2} + a_{3} q_{3} - b_{2} + b_{3} ) , \\
\dot{q}_{2} &= q_{2} ( \hphantom{-} a_{1} q_{1} - a_{3} q_{3} + b_{1} - b_{3} ) , \\
\dot{q}_{3} &= q_{3} (           -  a_{1} q_{1} + a_{2} q_{2} - b_{1} + b_{2} ) , \\
\end{aligned}
```
which can be written in Poisson-form as
```math
\dot{q} = P(q) \nabla H(q) ,
```
with Poisson matrix
```math
P(q) = \begin{pmatrix}
            0 & - q_{1} q_{2} & \hphantom{-} q_{1} q_{3} \\
\hphantom{-} q_{1} q_{2} &             0 & - q_{2} q_{3} \\
- q_{1} q_{3} & \hphantom{-} q_{2} q_{3} & 0             \\
\end{pmatrix} ,
```
and Hamiltonian
```math
H(q) = a_{1} q_{1} + a_{2} q_{2} + a_{3} q_{3} + b_{1} \ln q_{1} + b_{2} \ln q_{2} + b_{3} \ln q_{3} .
```

References:
* A. M. Perelomov. Selected topics on classical integrable systems,
  Troisième cycle de la physique, expanded version of lectures
  delivered in May 1995.

* Yuri B. Suris. Integrable discretizations for lattice systems: local
  equations of motion and their Hamiltonian properties,
  Rev. Math. Phys. 11, pp. 727–822, 1999.
"""
module LotkaVolterra3d

    using GeometricBase
    using GeometricEquations
    using GeometricSolutions
    using Parameters

    export odeproblem

    export hamiltonian, casimir
    export compute_energy_error, compute_casimir_error


    const Δt = 0.01
    const nt = 1000
    const timespan = (0.0, Δt*nt)

    const default_parameters = (A1 = 1.0, A2 = 1.0, A3 = 1.0, B1 = 0.0, B2 = 1.0, B3 = 1.0)
    const reference_solution = [0.39947308320241187, 1.9479527336244262, 2.570183075433086]


    function v₁(t, q, params)
        @unpack A1, A2, A3, B1, B2, B3 = params
        q[1] * ( - A2 * q[2] + A3 * q[3] + B2 - B3)
    end

    function v₂(t, q, params)
        @unpack A1, A2, A3, B1, B2, B3 = params
        q[2] * ( + A1 * q[1] - A3 * q[3] - B1 + B3)
    end

    function v₃(t, q, params)
        @unpack A1, A2, A3, B1, B2, B3 = params
        q[3] * ( - A1 * q[1] + A2 * q[2] + B1 - B2)
    end


    const X₀ = 1.0
    const Y₀ = 1.0
    const Z₀ = 2.0
    const q₀ = [X₀, Y₀, Z₀]
    const v₀ = [v₁(0, q₀, default_parameters), v₂(0, q₀, default_parameters), v₃(0, q₀, default_parameters)]


    function hamiltonian(t, q, params)
        @unpack A1, A2, A3, B1, B2, B3 = params
        A1*q[1] + A2*q[2] + A3*q[3] - B1*log(q[1]) - B2*log(q[2]) - B3*log(q[3])
    end

    hamiltonian_iode(v, t, q, params) = hamiltonian(t, q, params)

    function casimir(t, q, params)
        log(q[1]) + log(q[2]) + log(q[3])
    end


    function lotka_volterra_3d_v(v, t, q, params)
        v[1] = v₁(t, q, params)
        v[2] = v₂(t, q, params)
        v[3] = v₃(t, q, params)
        nothing
    end


    function odeproblem(q₀=q₀; timespan=timespan, timestep=Δt, parameters=default_parameters)
        ODEProblem(lotka_volterra_3d_v, timespan, timestep, q₀; parameters=parameters, invariants=(h=hamiltonian,))
    end


    function compute_energy_error(t::TimeSeries, q::DataSeries{T}, params) where {T}
        h = DataSeries(T, ntime(q))
        e = DataSeries(T, ntime(q))

        for i in axes(q,1)
            h[i] = hamiltonian(t[i], q[i], params)
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_casimir_error(t::TimeSeries, q::DataSeries{T}, params) where {T}
        c = DataSeries(T, ntime(q))
        e = DataSeries(T, ntime(q))

        for i in axes(q,1)
            c[i] = casimir(t[i], q[i], params)
            e[i] = (c[i] - c[0]) / c[0]
        end

        (c, e)
    end

end
