# Mathematical Pendulum

The mathematical pendulum consists of a point mass ``m`` attached to the origin ``(x,y) = (0,0)`` by a massless rod of length ``l``.
All motion is assumed to be frictionless and confined to a plane, and the pendulum moves in a homogeneous gravitational field of strength ``g``.

```@example
HTML("""<object type="image/svg+xml" class="display-light-only" data=$(joinpath(Main.buildpath, "images/pendulum.png"))></object>""") # hide
```

```@example
HTML("""<object type="image/svg+xml" class="display-dark-only" data=$(joinpath(Main.buildpath, "images/pendulum_dark.png"))></object>""") # hide
```

The configuration of the system is described by a single angle ``\theta``, measured from the upward vertical, in terms of which the cartesian coordinates of the mass are
```math
x = l \sin\theta , \qquad
y = l \cos\theta .
```

The Lagrangian is the difference of the kinetic energy ``T = \tfrac{1}{2} m l^2 \dot\theta^2`` and the potential energy ``V = m g y = m g l \cos\theta``,
```math
L (\theta, \dot\theta) = \frac{1}{2} m l^2 \dot\theta^2 - m g l \cos\theta .
```

The canonical conjugate momentum is obtained from the Lagrangian as
```math
p = \frac{\partial L}{\partial \dot\theta} = m l^2 \dot\theta ,
```
so that the Hamiltonian, obtained via the Legendre transform ``H = \dot\theta \, p - L``, reads
```math
H (\theta, p) = \frac{p^2}{2 m l^2} + m g l \cos\theta .
```

Hamilton's equations of motion are therefore
```math
\dot\theta =   \frac{\partial H}{\partial p}      = \frac{p}{m l^2} , \qquad
\dot p     = - \frac{\partial H}{\partial \theta} = m g l \sin\theta ,
```
which are equivalent to the familiar second-order form
```math
\ddot\theta = \frac{g}{l} \sin\theta .
```

The default parameters ``l = m = g = 1`` reduce the Hamiltonian to the normalised form ``H(\theta, p) = \tfrac{1}{2} p^2 + \cos\theta``.
Note that, with this sign convention, ``\theta = 0`` (the upright position) is an unstable equilibrium, whereas ``\theta = \pi`` (the mass hanging straight down) is stable.

The phase-space portrait of the Hamiltonian, with its characteristic separatrix between the librating (oscillating) and rotating regimes, is shown below.

```@eval
using GeometricProblems.Pendulum
using CairoMakie

H(q, p) = hamiltonian(0, q, p, Pendulum.default_parameters())

fig = Figure(size = (700, 400), figure_padding = 5, fontsize = 20)
ax = Axis(
    fig[1, 1],
    xlabel = "θ",
    ylabel = "p",
    title = "Pendulum Hamiltonian",
)

zq = range(-2π, 2π; length = 200)
zp = range(-3, 3; length = 200)
q = (zq' .* one.(zp))[:]
p = (one.(zq)' .* zp)[:]

co = contourf!(ax, q, p, H.(q, p); levels = 20)
Colorbar(fig[1, 2], co)

xlims!(ax, minimum(zq), maximum(zq))
ylims!(ax, minimum(zp), maximum(zp))

save("pendulum-hamiltonian.svg", fig)

nothing
```

![](pendulum-hamiltonian.svg)

## Hamiltonian dynamics

The canonical Hamiltonian dynamics is provided by `hodeproblem`, which can be integrated with any of the integrators from `GeometricIntegrators`.

```@example pendulum
using GeometricProblems.Pendulum
using GeometricIntegrators
using CairoMakie
using Logging # hide
Logging.disable_logging(Logging.Warn) # hide

prob = hodeproblem([1.0], [0.0]; timespan = (0.0, 20.0), timestep = 0.05)
sol = integrate(prob, ImplicitMidpoint())

t = [sol.t[n]    for n in axes(sol.q, 1)]
θ = [sol.q[n][1] for n in axes(sol.q, 1)]
p = [sol.p[n][1] for n in axes(sol.q, 1)]

fig = Figure(size = (800, 400), fontsize = 18)
ax = Axis(fig[1, 1], xlabel = "t", ylabel = "θ, p", title = "Hamiltonian pendulum")
lines!(ax, t, θ, label = "θ")
lines!(ax, t, p, label = "p")
axislegend(ax)
fig
```

## Ensembles

The function `hodeensemble` constructs an `HODEEnsemble` holding many instances of the problem at once.
In its default form it samples a grid of initial conditions over given bounds in ``(\theta, p)``, sharing a single set of parameters.
Integrating the ensemble and plotting all trajectories in phase space reproduces the portrait above:

```@example pendulum
ens = hodeensemble([-3.0], [3.0], [-2.5], [2.5], [7], [5]; timespan = (0.0, 10.0), timestep = 0.05)
sols = integrate(ens, ImplicitMidpoint())

fig = Figure(size = (600, 500), fontsize = 18)
ax = Axis(fig[1, 1], xlabel = "θ", ylabel = "p", title = "Ensemble over initial conditions")
for i in 1:length(sols)
    lines!(ax,
        [sols[i].q[n][1] for n in axes(sols[i].q, 1)],
        [sols[i].p[n][1] for n in axes(sols[i].q, 1)])
end
fig
```

The pendulum length ``l`` (as well as the mass ``m`` and the gravitational acceleration ``g``) is a parameter of the problem.
A second method of `hodeensemble` takes a single initial condition together with a vector of parameter sets, which is convenient for studying how the dynamics depends on the parameters.
Here we fix the initial condition and vary the length ``l``, which changes the oscillation period ``\propto \sqrt{l / g}``:

```@example pendulum
lengths = [0.5, 1.0, 2.0, 4.0]
params  = [(l = L, m = 1.0, g = 1.0) for L in lengths]

ens  = hodeensemble([2.0], [0.0], params; timespan = (0.0, 20.0), timestep = 0.05)
sols = integrate(ens, ImplicitMidpoint())

fig = Figure(size = (800, 400), fontsize = 18)
ax = Axis(fig[1, 1], xlabel = "t", ylabel = "θ", title = "Effect of the pendulum length l")
for i in 1:length(sols)
    lines!(ax,
        [sols[i].t[n]    for n in axes(sols[i].q, 1)],
        [sols[i].q[n][1] for n in axes(sols[i].q, 1)];
        label = "l = $(lengths[i])")
end
axislegend(ax)
fig
```

## Library functions

```@docs
GeometricProblems.Pendulum
```

```@autodocs
Modules = [GeometricProblems.Pendulum]
Order   = [:constant, :type, :macro, :function]
```
