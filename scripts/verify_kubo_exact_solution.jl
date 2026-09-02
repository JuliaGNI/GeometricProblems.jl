#
# Verification of the Kubo oscillator's closed-form solutions.
#
#     julia --project=. scripts/verify_kubo_exact_solution.jl
#
# `exact_solution` and `exact_mean_energy` are reference solutions: everything measured against
# them — a convergence order, an energy drift — is only as trustworthy as they are, and neither
# can be checked by reading it. The four claims settled here are
#
#   1. the pathwise solution solves the deterministic damped oscillator in the reparametrised
#      time θ, which with the θ-reparametrisation is what makes it the solution of the SDE;
#   2. it reproduces the initial condition, and t₀ enters only as the elapsed time t - t₀;
#   3. without damping it conserves the energy exactly along every sample path;
#   4. `exact_mean_energy` is the average of that pathwise energy over the Wiener increment.
#
# Exits non-zero if any check fails.
#

using GeometricProblems.KuboOscillator

const failures = String[]

function check(claim, passed, detail)
    println(passed ? "  PASS  " : "  FAIL  ", claim, "   ", detail)
    passed || push!(failures, claim)
    passed
end

header(title) = println("\n", title, "\n", "-"^length(title))

# The pathwise solution at ν = 0 is the deterministic one as a function of θ, since then
# θ = t - t₀. That is the function the ODE residual below is taken of.
deterministic(θ, q₀, p₀, γ) = exact_solution(θ, 0.0, q₀, p₀, 0.0, (ν = 0.0, γ = γ))

const q₀ = 2.0
const p₀ = 0.3
const H₀ = (q₀^2 + p₀^2) / 2

# The module defaults plus two more strongly damped settings, so that a check is not confined to
# the γ = 0.001 of §4.1, where every damped term is nearly invisible.
const PARAMETER_SETS = (KuboOscillator.default_parameters(),
    KuboOscillator.damped_parameters(),
    (ν = 0.5, γ = 0.3),
    (ν = 0.2, γ = 1.0))

# 1. The deterministic ODE, by central differences
#
# q' = p, p' = -q - γp. The step is chosen so that the central-difference truncation error,
# O(h²) times a third derivative of order one, sits near the round-off floor.

header("1. the pathwise solution solves q' = p, p' = -q - γp in θ")

for γ in (0.0, 0.001, 0.5, 1.0, 1.5, 1.9999)
    h = 1e-5
    residual = 0.0

    for θ in (0.0, 0.7, 3.1, 11.0)
        qm, pm = deterministic(θ - h, q₀, p₀, γ)
        qp, pp = deterministic(θ + h, q₀, p₀, γ)
        q, p = deterministic(θ, q₀, p₀, γ)

        residual = max(residual, abs((qp - qm) / 2h - p))
        residual = max(residual, abs((pp - pm) / 2h - (-q - γ * p)))
    end

    check("γ = $γ", residual < 1e-8, "max residual $(round(residual, sigdigits = 2))")
end

# 2. Initial condition and the role of t₀

header("2. the initial condition is reproduced and t₀ is the elapsed time")

for params in PARAMETER_SETS
    for t₀ in (0.0, 1.5, -2.0)
        q, p = exact_solution(t₀, 0.0, q₀, p₀, t₀, params)
        check("(q, p) at t = t₀ = $t₀, ν = $(params.ν)",
            abs(q - q₀) < 1e-14 && abs(p - p₀) < 1e-14,
            "error $(round(max(abs(q - q₀), abs(p - p₀)), sigdigits = 2))")
    end

    # Shifting t and t₀ together must not move the solution.
    shifted = exact_solution(4.5 + 1.5, 0.4, q₀, p₀, 1.5, params)
    base = exact_solution(4.5, 0.4, q₀, p₀, 0.0, params)
    check("t₀ shift invariance, ν = $(params.ν)",
        maximum(abs.(shifted .- base)) < 1e-14,
        "difference $(round(maximum(abs.(shifted .- base)), sigdigits = 2))")
end

# 3. Pathwise energy conservation without damping

header("3. the undamped energy is conserved along every sample path")

let params = KuboOscillator.default_parameters()
    for (t, W) in ((0.0, 0.0), (1.3, -0.8), (7.0, 2.4), (50.0, -11.7))
        q, p = exact_solution(t, W, q₀, p₀, 0.0, params)
        H = (q^2 + p^2) / 2
        check("t = $t, W = $W", abs(H - H₀) < 1e-13,
            "|H - H₀| = $(round(abs(H - H₀), sigdigits = 2))")
    end
end

# 4. E(H) against a quadrature of the pathwise energy
#
# W(t) - W(t₀) is Gaussian with variance t - t₀, so E(H) is the integral of the pathwise energy
# against that density. The integrand is smooth and Gaussian-weighted, so the trapezoidal rule
# converges spectrally; ±10 standard deviations over 20001 nodes puts the quadrature error well
# below the tolerance asserted here.

header("4. exact_mean_energy is the noise average of the pathwise energy")

function quadrature_mean_energy(t, params; nodes = 20001, width = 10)
    σ = sqrt(t)
    z = range(-width * σ, width * σ; length = nodes)
    density = exp.(-(z .^ 2) ./ (2 * σ^2)) ./ (σ * sqrt(2π))
    energy = [(qp = exact_solution(t, W, q₀, p₀, 0.0, params); (qp[1]^2 + qp[2]^2) / 2)
              for W in z]
    step(z) * sum(density .* energy)
end

for params in PARAMETER_SETS
    for t in (0.5, 2.0, 10.0, 50.0)
        closed = exact_mean_energy(t, q₀, p₀, 0.0, params)
        quad = quadrature_mean_energy(t, params)
        error = abs(closed - quad) / abs(quad)
        check("ν = $(params.ν), γ = $(get(params, :γ, 0.0)), t = $t", error < 1e-10,
            "relative error $(round(error, sigdigits = 2))")
    end
end

println()
if isempty(failures)
    println("verify_kubo_exact_solution.jl: all checks passed")
else
    println("verify_kubo_exact_solution.jl: $(length(failures)) check(s) FAILED")
    foreach(f -> println("  ", f), failures)
    exit(1)
end
