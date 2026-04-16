# Harmonic Oscillator

The equations of motion of a harmonic Oscillator with spring constant $k$ and mass $m$ are given by
```math
m \, \ddot{x} (t) = - k \, x (t) , \qquad
x(t_0) = x_0 , \qquad
\dot{x} (t_0) = \dot{x}_0 .
```
Its solution is obtained as
```math
x (t) = A \, \cos( \omega (t - t_0) - \varphi )
```
with the frequency $\omega = \sqrt{k / m}$, amplitude $A = x_0 \, \sqrt{1 + \dot{x}_0^2 / (m k x_0^2)}$, and phase $\varphi = \arctan(\dot{x}_0 / (m \omega x_0))$.
The energy of the system is
```math
E(x, \dot{x}) = \frac{m}{2} \dot{x}^2 + \frac{k}{2} x^2 .
```

## Hamiltonian Formulation

In canonical Hamiltonian coordinates $(q = x, \, p = m \dot{x})$, the Hamiltonian of the harmonic oscillator is given by
```math
H(q, p) = \frac{1}{2m} p^2 + \frac{k}{2} q^2 ,
```
such that Hamilton's equations of motion read
```math
\dot{q} =   \frac{\partial H}{\partial p}(q, p) = \frac{p}{m} , \qquad
\dot{p} = - \frac{\partial H}{\partial q}(q, p) = - k q , \qquad
q(t_0) = q_0 , \qquad
p(t_0) = p_0 .
```
and their solution is given by
```math
q (t) = A \, \cos( \omega (t - t_0) - \varphi ) , \qquad
p (t) = - \omega A \, \sin( \omega (t - t_0) - \varphi ) ,
```
with
```math
A = q_0 \, \sqrt{1 + p_0^2 / (k q_0^2)} , \qquad
\varphi = \arctan(p_0 / (\omega q_0)) .
```

```@eval
using GeometricProblems.HarmonicOscillator
using CairoMakie

H(q,p) = hamiltonian(0, q, p, HarmonicOscillator.default_parameters)

levels = [0.0, 0.05, 0.25, 0.65, 1.2, collect(2:7)...]
ticks = [-2, -1, 0, +1, +2]
fig = Figure(size = (600, 600), figure_padding = 5, fontsize = 24)
ax = Axis(
    fig[1,1],
    aspect = 1,
    xlabel = "q",
    ylabel = "p",
    xticks = ticks,
    yticks = ticks,
    title = "Harmonic Oscillator Hamiltonian",
)

z = range(-2.5, +2.5; length = 100)
q = (z' .* one.(z))[:]
p = (one.(z)' .* z)[:]

co = contourf!(ax, q, p, H.(q,p); levels = levels)
cb = Colorbar(fig[1, 2], co)

xlims!(ax, minimum(z), maximum(z))
ylims!(ax, minimum(z), maximum(z))

save("harmonic-oscillator-hamiltonian.svg", fig)

nothing
```

![](harmonic-oscillator-hamiltonian.svg)

## Lagrangian Formulation


## Dynamics


## Library

```@autodocs
Modules = [GeometricProblems.HarmonicOscillator]
Order   = [:constant, :type, :macro, :function]
```
