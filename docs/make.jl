using Documenter
using GeometricProblems
using DocumenterCitations 

# if the docs are generated with github actions, then this changes the path; see: https://github.com/JuliaDocs/Documenter.jl/issues/921 
const buildpath = haskey(ENV, "CI") ? ".." : ""

const bib = CitationBibliography(joinpath(@__DIR__, "src", "GeometricProblems.bib"))

makedocs(;
    plugins = [bib],
    sitename = "GeometricProblems.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = [
            "assets/extra_styles.css",
        ],
    ),
    pages = ["Home" => "index.md",
             "Diagnostics"                 => "diagnostics.md",
             "Coupled Harmonic Oscillator" => "coupled_harmonic_oscillator.md",
             "Double Pendulum"             => "double_pendulum.md",
             "Harmonic Oscillator"         => "harmonic_oscillator.md",
             "Hénon-Heiles System"         => "henon_heiles.md",
             "Kepler Problem"              => "kepler_problem.md",
             "Lorenz Attractor"            => "lorenz_attractor.md",
             "Lotka-Volterra 2d"           => "lotka_volterra_2d.md",
             "Lotka-Volterra 3d"           => "lotka_volterra_3d.md",
             "Lotka-Volterra 4d"           => "lotka_volterra_4d.md",
             "Massless Charged Particle"   => "massless_charged_particle.md",
             "Mathematical Pendulum"       => "pendulum.md",
             "Nonlinear Oscillators"       => "nonlinear_oscillators.md",
             "Point Vortices"              => "point_vortices.md",
             "Inner Solar System"          => "inner_solar_system.md",
             "Outer Solar System"          => "outer_solar_system.md",
             "Toda Lattice"                => "toda_lattice.md",
            ]
)

deploydocs(
    repo   = "github.com/JuliaGNI/GeometricProblems.jl",
    devurl = "latest",
    devbranch = "main",
)
