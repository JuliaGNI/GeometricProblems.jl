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
\hphantom{-} 0 & - 1 \\
+ 1 & \hphantom{+} 0 \\
\end{pmatrix} \nabla \phi (x) .
```

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
so that the Euler-Lagrange equations read ``\Omega (x) \, \dot{x} = - \nabla \phi (x)``.

See [`GeometricProblems.MasslessChargedParticleSingular`](@ref) for a formulation with the same
magnetic field but a "singular" (one-component) vector potential, suitable for degenerate
variational integrators.

"""
module MasslessChargedParticle

    # components of the vector potential
    A₁(q, params) = - params[:A₀] * q[2] * (1 + q[1]^2 + q[2]^2) / 2
    A₂(q, params) = + params[:A₀] * q[1] * (1 + q[1]^2 + q[2]^2) / 2

    # z-componend of the magnetic field
    B(q, params) = params[:A₀] * (1 + 2 * q[1]^2 + 2 * q[2]^2)

    # derivatives of the one-form components
    dϑ₁dx₁(t, q, params) = - params[:A₀] * q[1] * q[2]
    dϑ₁dx₂(t, q, params) = - params[:A₀] * (1 + q[1]^2 + 3 * q[2]^2) / 2
    dϑ₂dx₁(t, q, params) = + params[:A₀] * (1 + 3 * q[1]^2 + q[2]^2) / 2
    dϑ₂dx₂(t, q, params) = + params[:A₀] * q[1] * q[2]


    include("massless_charged_particle_common.jl")

end
