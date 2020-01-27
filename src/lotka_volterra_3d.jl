@doc raw"""
    Lotka-Volterra model in 3D

The equations read
```math
\begin{align*}
\dot{q}_{1} &= q_{1} ( - A_{2} q_{2} + A_{3} q_{3} + B_{2} - B_{3} ) , &
\dot{q}_{2} &= q_{2} (   A_{1} q_{1} - A_{3} q_{3} - B_{1} + B_{3} ) , &
\dot{q}_{3} &= q_{3} ( - A_{1} q_{1} + A_{2} q_{2} + B_{1} - B_{2} ) ,
\end{align*}
```
which can be written in Poisson-form as
```math
\dot{q} = P(q) \nabla H(q) ,
```
with Poisson matrix
```math
P(q) = \begin{pmatrix}
            0 & - q_{1} q_{2} &   q_{1} q_{3} \\
  q_{1} q_{2} &             0 & - q_{2} q_{3} \\
- q_{1} q_{3} &   q_{2} q_{3} & 0             \\
\end{pmatrix} ,
```
and Hamiltonian
```math
H(q) = A_{1} q_{1} + A_{2} q_{2} + A_{3} q_{3} - B_{1} \ln q_{1} - B_{2} \ln q_{2} - B_{3} \ln q_{3} .
```

References:
    A. M. Perelomov. Selected topics on classical integrable systems,
                Troisième cycle de la physique, expanded version of lectures
                delivered in May 1995.

    Yuri B. Suris. Integrable discretizations for lattice systems: local
                equations of motion and their Hamiltonian properties,
                Rev. Math. Phys. 11, pp. 727–822, 1999.
"""
module LotkaVolterra3d

    using GeometricIntegrators.Equations
    using GeometricIntegrators.Solutions

    export lotka_volterra_3d_ode

    export hamiltonian, casimir
    export compute_energy_error, compute_casimir_error


    Δt = 0.01
    nt = 1000

    const A1=+1.0
    const A2=+1.0
    const A3=+1.0

    const B1= 0.0
    const B2=+1.0
    const B3=+1.0

    const X0=1.0
    const Y0=1.0
    const Z0=2.0


    function v₁(t, q)
        q[1] * ( - A2 * q[2] + A3 * q[3] + B2 - B3)
    end

    function v₂(t, q)
        q[2] * ( + A1 * q[1] - A3 * q[3] - B1 + B3)
    end

    function v₃(t, q)
        q[3] * ( - A1 * q[1] + A2 * q[2] + B1 - B2)
    end


    const q₀=[X0, Y0, Z0]
    const v₀=[v₁(0, q₀), v₂(0, q₀), v₃(0, q₀)]


    function hamiltonian(t, q)
        A1*q[1] + A2*q[2] + A3*q[3] - B1*log(q[1]) - B2*log(q[2]) - B3*log(q[3])
    end

    function casimir(t, q)
        log(q[1]) + log(q[2]) + log(q[3])
    end


    function lotka_volterra_3d_v(t, q, v)
        v[1] = v₁(t, q)
        v[2] = v₂(t, q)
        v[3] = v₃(t, q)
        nothing
    end


    function lotka_volterra_3d_ode(q₀=q₀)
        ODE(lotka_volterra_3d_v, q₀; h=hamiltonian)
    end


    function compute_energy_error(t, q::DataSeries{T}) where {T}
        h = SDataSeries(T, q.nt)
        e = SDataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i])
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_casimir_error(t, q::DataSeries{T}) where {T}
        c = SDataSeries(T, q.nt)
        e = SDataSeries(T, q.nt)

        for i in axes(q,2)
            c[i] = casimir(t[i], q[:,i])
            e[i] = (c[i] - c[0]) / c[0]
        end

        (c, e)
    end

end
