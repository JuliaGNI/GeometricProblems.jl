#
# Verification of the shadow Hamiltonian of symplectic Euler on the harmonic oscillator.
#
#     julia --project=examples scripts/verify_symplectic_euler_shadow_hamiltonian.jl
#
# `examples/harmonic_oscillator.jl` plots a shadow Hamiltonian H̃ with the run of one of the two
# symplectic Euler methods. For H = p²/(2m) + k q²/2 it is, to second order in the step h,
#
#     H̃ = H + s h/2 H_p H_q + h²/12 (H_pp H_q² + H_qq H_p²),
#     H_p = p/m,  H_q = k q,  H_pp = 1/m,  H_qq = k,
#
# with s = +1 for the method that evaluates at (q_{n+1}, p_n), `SymplecticEulerA` in
# GeometricIntegrators 0.19, and s = -1 for the one at (q_n, p_{n+1}), `SymplecticEulerB`. A
# correct H̃ is constant along a trajectory up to the third order of the step, so the spread of H̃
# along a run falls by a factor 8 when h halves, and an error in either term shows as a factor 4 or
# less. The claims settled here are
#
#   1. the sign s and the h² term above hold for `SymplecticEulerA` and `SymplecticEulerB`;
#   2. the opposite sign fails for each method, so the check can tell the two letters apart;
#   3. the h² term `(k² q² + p²) / 12`, which lacks the factors 1/m and k/m², is wrong at order h².
#
# Exits non-zero if any check fails.
#

using GeometricIntegrators
using GeometricProblems.HarmonicOscillator

const failures = String[]

function check(claim, passed, detail)
    println(passed ? "  PASS  " : "  FAIL  ", claim, "   ", detail)
    passed || push!(failures, claim)
end

H(q, p, params) = hamiltonian(0, q, p, params)
H_p(q, p, params) = p / params.m
H_q(q, p, params) = params.k * q
h2_correct(q, p, params) = (params.k^2 * q^2 / params.m + params.k * p^2 / params.m^2) / 12
h2_unscaled(q, p, params) = (params.k^2 * q^2 + p^2) / 12

function shadow(s, h2, params)
    (q, p, h) -> H(q, p, params) + s * h / 2 * H_p(q, p, params) * H_q(q, p, params) +
                 h^2 * h2(q, p, params)
end

# Spread of f along a run of `method` with step h.
function spread(method, h, f, params)
    prob = podeproblem(timespan = (0.0, 40.0), timestep = h, parameters = params)
    sol = integrate(prob, method)
    vals = [f(sol.q[n][1], sol.p[n][1], h) for n in eachindex(sol.q)]
    maximum(vals) - minimum(vals)
end

# Factor by which the spread falls when h halves: 8 for a third-order error, 4 for a second-order one.
function ratio(method, f, params)
    spread(method, 0.04, f, params) / spread(method, 0.02, f, params)
end

# The default m = 1 cannot tell a factor 1/m from its absence, so a second set has m ≠ 1.
for params in (default_parameters(), (m = 2.5, k = 0.7, ω = √(0.7 / 2.5)))
    println("Shadow Hamiltonian of symplectic Euler, H = p²/(2m) + k q²/2, m = $(params.m), k = $(params.k)")
    for (name, method, s) in (("SymplecticEulerA", SymplecticEulerA(), +1),
        ("SymplecticEulerB", SymplecticEulerB(), -1))
        claim(text) = "$name, m = $(params.m): $text"
        r_correct = ratio(method, shadow(s, h2_correct, params), params)
        r_unscaled = ratio(method, shadow(s, h2_unscaled, params), params)
        spread_wrong = spread(method, 0.02, shadow(-s, h2_correct, params), params)
        spread_right = spread(method, 0.02, shadow(s, h2_correct, params), params)
        check(claim("s = $s and the h² term hold"), r_correct > 7,
            "spread ratio $(round(r_correct; digits = 2))")
        check(claim("the h² term (k² q² + p²) / 12 is wrong at order h²"), r_unscaled < 5,
            "spread ratio $(round(r_unscaled; digits = 2))")
        check(claim("the opposite sign fails"), spread_wrong > 100 * spread_right,
            "spread $(spread_wrong) against $(spread_right)")
    end
end

if isempty(failures)
    println("COMPLETED")
else
    println("\n", length(failures), " check(s) failed")
    exit(1)
end
