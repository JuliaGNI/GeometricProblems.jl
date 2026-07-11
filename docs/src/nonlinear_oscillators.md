# Nonlinear Oscillators

Phase-space trajectories of the four oscillators, obtained by integrating their default Hamiltonian
problems:

```@example
using GeometricIntegrators: integrate, Gauss
using CairoMakie

import GeometricProblems.DuffingOscillator as duffing
import GeometricProblems.LennardJonesOscillator as lj
import GeometricProblems.MorseOscillator as morse
import GeometricProblems.MathewsLakshmananOscillator as ml

models = [("Duffing", duffing), ("Lennard-Jones", lj), ("Morse", morse), ("Mathews–Lakshmanan", ml)]

fig = Figure(size = (800, 700))
for (i, (name, m)) in enumerate(models)
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1
    ax = Axis(fig[row, col]; xlabel = "q", ylabel = "p", title = name)
    sol = integrate(m.hodeproblem(), Gauss(2))
    lines!(ax, sol.q[:, 1], sol.p[:, 1])
end
fig
```

## Duffing Oscillator

```@docs
GeometricProblems.DuffingOscillator
```

```@autodocs
Modules = [GeometricProblems.DuffingOscillator]
Order   = [:constant, :type, :macro, :function]
```

## Lennard-Jones Oscillator

```@docs
GeometricProblems.LennardJonesOscillator
```

```@autodocs
Modules = [GeometricProblems.LennardJonesOscillator]
Order   = [:constant, :type, :macro, :function]
```

## Mathews-Lakshmanan Oscillator

```@docs
GeometricProblems.MathewsLakshmananOscillator
```

```@autodocs
Modules = [GeometricProblems.MathewsLakshmananOscillator]
Order   = [:constant, :type, :macro, :function]
```

## Morse Oscillator

```@docs
GeometricProblems.MorseOscillator
```

```@autodocs
Modules = [GeometricProblems.MorseOscillator]
Order   = [:constant, :type, :macro, :function]
```
