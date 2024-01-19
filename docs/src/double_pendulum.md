# Double Pendulum

The double pendulum consists of two pendula, one attached to the origin at ``(x,y) = (0,0)``, and the second attached to the first.
Each pendulum consists of a point mass ``m_i`` attached to a massless rod of length ``l_i`` with ``i \in (1,2)``.
All motion is assumed to be frictionless.

![Double Pendulum](images/double-pendulum.png)

The dynamics of the system is most naturally described in terms of the angles ``\theta_i`` between the rods ``l_i`` and the vertical axis ``y``.
In terms of these angles, the cartesian coordinates are given by
```math
\begin{align*}
x_1 &= l_1 \sin\theta_1 , \\
x_2 &= l_1 \sin\theta_1 + l_2 \sin\theta_2 , \\
y_1 &= - l_1 \cos\theta_1 , \\
y_2 &= -l_1 \cos\theta_1 - l_2 \cos\theta_2 .
\end{align*}
```

In terms of the generalized coordinates ``\theta_i``, the Lagrangian reads
```math
\begin{align*}
L (\theta_1, \theta_2, \dot{\theta}_1, \dot{\theta}_2)
 = \frac{1}{2} (m_1 + m_2) l_1^2 \dot{\theta}_1^2 
 &+ \frac{1}{2} m_2 l_2^2 \dot{\theta}_2^2
 + m_2 l_1 l_2 \dot{\theta}_1 \dot{\theta}_2 \cos(\theta_1 - \theta_2) \\
 &+ g (m_1 + m_2) l_1 \cos\theta_1
 + g m_2 l_2 \cos\theta_2 .
\end{align*}
```

The canonical conjugate momenta ``p_i`` are obtained from the Lagrangian as
```math
\begin{align*}
p_1 &= \frac{\partial L}{\partial \dot{\theta}_1} = (m_1 + m_2) l_1^2 \dot{\theta}_1 + m_2 l_1 l_2 \dot{\theta}_2 \cos(\theta_1 - \theta_2), \\
p_2 &= \frac{\partial L}{\partial \dot{\theta}_2} = m_2 l_2^2 \dot{\theta}_2 + m_2 l_1 l_2 \dot{\theta}_1 \cos(\theta_1 - \theta_2) .
\end{align*}
```

After solving these relations for the generalized velocities ``\dot{\theta}_i``,
```math
\begin{align*}
\dot{\theta}_1 &= \frac{l_2 p_{\theta_1} - l_1 p_{\theta_2} \cos(\theta_1 - \theta_2)}{l_1^2 l_2 \left[ m_1 + m_2 \sin^2(\theta_1 - \theta_2) \right] } \\
\dot{\theta}_2 &= \frac{(m_1 + m_2) l_1 p_{\theta_2} - m_2 l_2 p_{\theta_1} \cos(\theta_1 - \theta_2)}{m_2 l_1 l_2^2 \left[ m_1 + m_2 \sin^2 (\theta_1 - \theta_2) \right] } ,
\end{align*}
```
the Hamiltonian can be obtained via the Legendre transform,
```math
H = \sum_{i=1}^2 \dot{\theta}_i p_i - L ,
```
as
```math
\begin{align*}
H &= \frac{m_2 l_2^2 p^2_{\theta_1} + (m_1 + m_2) l_1^2 p^2_{\theta_2} - 2 m_2 l_1 l_2 p_{\theta_1} p_{\theta_2} \cos(\theta_1 - \theta_2)}{2 m_2 l_1^2 l_2^2 \left[ m_1 + m_2 \sin^2(\theta_1 - \theta_2) \right] } \\
  & \qquad\qquad \vphantom{\frac{l}{l}} - g (m_1 + m_2) l_1 \cos\theta_1 - g m_2 l_2 \cos\theta_2 .
\end{align*}
```


## Library functions

```@docs
GeometricProblems.DoublePendulum
```

```@autodocs
Modules = [GeometricProblems.DoublePendulum]
Order   = [:constant, :type, :macro, :function]
```
