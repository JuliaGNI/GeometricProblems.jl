@doc raw"""
# Massless charged particle in 2D with "singular" vector potential

The Lagrangian is given by
```math
L(x, \dot{x}) = A(x) \cdot \dot{x} - \phi (x) ,
```
with the "singular" (one-component) magnetic vector potential
```math
A(x) = - A_0 \, x_2 \, \big( 1 + 2 x_1^2 + \tfrac{2}{3} x_2^2 \big) \begin{pmatrix}
1 \\
0 \\
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
\hphantom{-} 0 & - 1 \\
+ 1 & \hphantom{+} 0 \\
\end{pmatrix} \nabla \phi (x) .
```

This is the same physical system as [`GeometricProblems.MasslessChargedParticle`](@ref): the
magnetic field ``B``, electric field ``E`` and hence the dynamics are identical. Only the vector
potential differs by a gauge transformation, chosen here so that its second component vanishes
(``A_2 = 0``). This "singular" one-form yields the same Euler-Lagrange equations but a different
variational integrator, and is the form required by the degenerate variational integrator (DVI)
method of GeometricIntegrators.

The problem is available in explicit (`odeproblem`), implicit (`iodeproblem`, `idaeproblem`,
`idaeproblem_spark`) and Lagrangian (`lodeproblem`, `ldaeproblem`) form. The latter two are built
from the Lagrangian ``L(x, \dot{x}) = \vartheta (x) \cdot \dot{x} - H(x)`` with one-form
``\vartheta = A`` and the symplectic two-form
```math
\Omega_{ij} (x) = \frac{\partial \vartheta_i}{\partial x_j} - \frac{\partial \vartheta_j}{\partial x_i}
= \begin{pmatrix}
\hphantom{+} 0 & - B(x) \\
+ B(x) & \hphantom{-} 0 \\
\end{pmatrix} ,
```
so that the Euler-Lagrange equations read ``\Omega (x) \, \dot{x} = - \nabla \phi (x)``. Since the
gauge transformation shifts ``\vartheta`` only by an exact one-form, ``\Omega`` is the same as for
[`GeometricProblems.MasslessChargedParticle`](@ref), while the Lagrangian differs.

"""
module MasslessChargedParticleSingular

    # components of the vector potential (singular gauge, A₂ = 0)
    A₁(q, params) = - params[:A₀] * q[2] * (1 + 2 * q[1]^2 + 2 * q[2]^2 / 3)
    A₂(q, params) = zero(eltype(q))

    # z-componend of the magnetic field (same as the standard vector potential)
    B(q, params) = params[:A₀] * (1 + 2 * q[1]^2 + 2 * q[2]^2)

    # derivatives of the one-form components
    dϑ₁dx₁(t, q, params) = - 4 * params[:A₀] * q[1] * q[2]
    dϑ₁dx₂(t, q, params) = - params[:A₀] * (1 + 2 * q[1]^2 + 2 * q[2]^2)
    dϑ₂dx₁(t, q, params) = zero(eltype(q))
    dϑ₂dx₂(t, q, params) = zero(eltype(q))


    include("massless_charged_particle_common.jl")

end
