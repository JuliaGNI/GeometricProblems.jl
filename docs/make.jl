using Documenter

makedocs(
    sitename = "GeometricProblems.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = ["Home" => "index.md",
             "Diagnostics"          => "diagnostics.md",
             "Exponential Growth"   => "exponential_growth.md",
             "Harmonic Oscillator"  => "harmonic_oscillator.md",
             "Lorenz Attractor"     => "lorenz_attractor.md",
             "Lotka-Volterra 2d"    => "lotka_volterra_2d.md",
             "Lotka-Volterra 2d (singular)"  => "lotka_volterra_2d_singular.md",
             "Lotka-Volterra 3d"    => "lotka_volterra_3d.md",
             "Pendulum"             => "pendulum.md",
             "Point Vortices"       => "point_vortices.md"
            ]
)

deploydocs(
    repo   = "github.com/DDMGNI/GeometricProblems.jl.git",
)
