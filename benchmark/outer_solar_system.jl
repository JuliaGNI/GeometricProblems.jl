# Benchmark for the OuterSolarSystem equations of motion.
#
# Run with
#
#     julia --project=. benchmark/outer_solar_system.jl
#
# It answers three questions:
#
#   1. What does symbolic construction cost, and how do `simplify` and `cse` affect it?
#   2. How large is the generated code, and how long does it take to compile and to evaluate,
#      compared with the hand-written vector fields?
#   3. Do the hand-written and generated functions agree?
#
# `simplify` and `cse` are keyword arguments of `EulerLagrange.LagrangianSystem` /
# `HamiltonianSystem`. `simplify` is cheap on the Lagrangian itself but rewrites the sum of
# N(N-1)/2 inverse-distance terms into a common-denominator form whose derivatives are far more
# expensive; `cse` turns on common subexpression elimination in `Symbolics.build_function`, which
# is what keeps each pairwise distance from being recomputed at every occurrence.

using EulerLagrange
using GeometricProblems.OuterSolarSystem
using LinearAlgebra
using Printf

const OSS = OuterSolarSystem
const Symbolics = EulerLagrange.Symbolics

const N_DIM = OSS.N_DIM
const ALL_MASSES = OSS.default_parameters().m
const G = OSS.default_parameters().G

"Parameters for the first `n` bodies."
sub_parameters(n) = (G=G, m=ALL_MASSES[1:n])

# `OuterSolarSystem.lagrangian` is fixed at `N_BODIES` bodies, so the scaling sweep below needs a
# version parameterised by `n`. It is otherwise identical, including the explicit scalar indexing
# (which `Symbolics.scalarize` handles far better than slicing a reshaped array).
function sub_lagrangian(t, q, w, params, n)
    kinetic = zero(eltype(w))
    for i in 1:n
        kinetic += params.m[i] * (w[N_DIM*i-2]^2 + w[N_DIM*i-1]^2 + w[N_DIM*i]^2)
    end
    potential = zero(eltype(q))
    for i in 2:n, j in 1:i-1
        Δ₁ = q[N_DIM*i-2] - q[N_DIM*j-2]
        Δ₂ = q[N_DIM*i-1] - q[N_DIM*j-1]
        Δ₃ = q[N_DIM*i] - q[N_DIM*j]
        potential += params.G * params.m[i] * params.m[j] / sqrt(Δ₁^2 + Δ₂^2 + Δ₃^2)
    end
    return kinetic / 2 + potential
end

"Initial positions and velocities for the first `n` bodies."
function sub_state(n)
    q = OSS.q₀[1:N_DIM*n]
    p = OSS.p₀[1:N_DIM*n]
    v = [p[N_DIM*i-2+k] / ALL_MASSES[i] for i in 1:n for k in 0:2]
    return q, p, v
end

function elapsed(f)
    GC.gc()
    start = time()
    result = f()
    return time() - start, result
end

"""
    best_of(f, samples, repetitions)

Seconds per call to `f`, taken as the best of `samples` batches of `repetitions` calls each.

A single force evaluation here is well under a microsecond, which is the resolution of `time()`, so
each batch has to loop inside the timed region. `f` must already have been called once, so that
compilation is not measured.
"""
function best_of(f, samples=20, repetitions=10_000)
    best = Inf
    for _ in 1:samples
        GC.gc()
        start = time()
        for _ in 1:repetitions
            f()
        end
        best = min(best, (time() - start) / repetitions)
    end
    return best
end

# Warm the compiler on the smallest case of the sweep itself, on all four branches. A trivial
# one-degree-of-freedom system is not enough: `simplify=true` on a sum of inverse distances goes
# through `PolyForm`, whose compilation would otherwise land entirely in the first table cell.
let n = 2
    t, x, v = lagrangian_variables(N_DIM * n)
    sparams = symbolize(sub_parameters(n))
    lag = sub_lagrangian(t, x, v, sparams, n)
    for simplify in (true, false), cse in (true, false)
        LagrangianSystem(lag, t, x, v, sparams; simplify=simplify, cse=cse)
    end
end

# ------------------------------------------------------------------------------------------------
println("\n", "="^96)
println("1. Symbolic construction time  (LagrangianSystem, seconds)")
println("="^96)
@printf("%8s %8s   %12s %12s %12s %12s\n",
    "bodies", "dof", "simp/cse", "simp/nocse", "nosimp/cse", "nosimp/nocse")

for n in 2:OSS.N_BODIES
    dof = N_DIM * n
    parameters = sub_parameters(n)
    times = Float64[]
    for simplify in (true, false), cse in (true, false)
        # `simplify=true, cse=false` is the pre-fix default and grows very steeply; skip the
        # combinations that would dominate the whole benchmark run.
        if simplify && n > 5
            push!(times, NaN)
            continue
        end
        seconds, _ = elapsed() do
            t, x, v = lagrangian_variables(dof)
            sparams = symbolize(parameters)
            LagrangianSystem(sub_lagrangian(t, x, v, sparams, n), t, x, v, sparams;
                simplify=simplify, cse=cse)
        end
        push!(times, seconds)
    end
    fmt(x) = isnan(x) ? "     skipped" : @sprintf("%12.2f", x)
    @printf("%8d %8d   %s %s %s %s\n", n, dof, fmt(times[1]), fmt(times[2]), fmt(times[3]), fmt(times[4]))
end
println("\n(`simp/nocse` was the default before the code-generation fixes; combinations with",
    "\n `simplify=true` are skipped above 5 bodies because they take several minutes.)")

# ------------------------------------------------------------------------------------------------
println("\n", "="^96)
println("2. Generated code size at $(N_DIM * OSS.N_BODIES) degrees of freedom  (characters)")
println("="^96)

let n = OSS.N_BODIES, dof = N_DIM * n
    parameters = sub_parameters(n)
    t, x, v = lagrangian_variables(dof)
    sparams = symbolize(parameters)
    lag = OSS.lagrangian(t, x, v, sparams)
    Ls = Symbolics.Num(Symbolics.scalarize(lag))

    Dt, Dx, Dv = EulerLagrange.lagrangian_derivatives(t, x, v)
    Dz = vcat(Dx, Dv)
    ϑ = [Symbolics.expand_derivatives(dv(Ls)) for dv in Dv]
    f = [Symbolics.expand_derivatives(dx(Ls)) for dx in Dx]

    # `ω = dθ_L` for the Lagrange one-form θ_L = ϑᵢ dqⁱ, which in the (x, v) coordinates spanned by
    # `Dz` has the components of ϑ along dx and none along dv. This mirrors what
    # `EulerLagrange.LagrangianSystem` builds; note that antisymmetrising the *full* gradient ∂L/∂z
    # instead — as EulerLagrange did before 0.5 — yields d(dL) ≡ 0.
    ϑz = vcat(ϑ, zero(ϑ))
    ω = Symbolics.expand_derivatives.([Dz[i](ϑz[j]) - Dz[j](ϑz[i]) for i in eachindex(Dz, ϑz), j in eachindex(Dz, ϑz)])

    equations = EulerLagrange.substitute_lagrangian_variables(
        (f=f, ϑ=ϑ, ω=ω), x, collect(Dt.(x)), v)

    Symbolics.@variables X[1:dof] V[1:dof]

    @printf("%6s %10s   %14s %14s %8s\n", "", "shape", "cse=false", "cse=true", "ratio")
    for key in (:f, :ϑ, :ω)
        expression = getproperty(equations, key)
        without = Symbolics.build_function(expression, t, X, V, sparams...; nanmath=false, cse=false)[2]
        with = Symbolics.build_function(expression, t, X, V, sparams...; nanmath=false, cse=true)[2]
        n_without, n_with = length(string(without)), length(string(with))
        @printf("%6s %10s   %14d %14d %8.2f\n",
            key, string(size(expression)), n_without, n_with, n_without / n_with)
    end
end

# ------------------------------------------------------------------------------------------------
println("\n", "="^96)
println("3. Evaluation cost at $(N_DIM * OSS.N_BODIES) degrees of freedom")
println("="^96)

let n = OSS.N_BODIES, dof = N_DIM * n
    parameters = sub_parameters(n)
    q, p, v = sub_state(n)

    build_seconds, lag_sys = elapsed(() -> OSS.lagrangian_system(parameters))
    @printf("  lagrangian_system(...)              %8.2f s   (simplify=false, cse=true)\n", build_seconds)
    build_seconds, ham_sys = elapsed(() -> OSS.hamiltonian_system(parameters))
    @printf("  hamiltonian_system(...)             %8.2f s\n\n", build_seconds)

    lagrangian_functions = functions(lag_sys)
    hamiltonian_functions = functions(ham_sys)

    out = zeros(dof)
    @printf("%26s %13s %13s %8s\n", "", "generated", "hand-written", "speedup")

    pairs = (
        ("f  (=-∂V/∂q, Lagrangian)",
            () -> lagrangian_functions.f(out, 0.0, q, v, parameters),
            () -> OSS.outer_solar_system_f(out, 0.0, q, v, parameters)),
        ("f  (=-∂H/∂q, Hamiltonian)",
            () -> hamiltonian_functions.f(out, 0.0, q, p, parameters),
            () -> OSS.outer_solar_system_f(out, 0.0, q, p, parameters)),
        ("v  (=∂H/∂p)",
            () -> hamiltonian_functions.v(out, 0.0, q, p, parameters),
            () -> OSS.outer_solar_system_v(out, 0.0, q, p, parameters)),
        ("ϑ  (=∂L/∂q̇)",
            () -> lagrangian_functions.ϑ(out, 0.0, q, v, parameters),
            () -> OSS.outer_solar_system_ϑ(out, 0.0, q, v, parameters)),
    )

    for (label, generated, handwritten) in pairs
        # the first call includes compilation of the runtime-generated function
        first_generated, _ = elapsed(generated)
        first_handwritten, _ = elapsed(handwritten)
        steady_generated = best_of(generated)
        steady_handwritten = best_of(handwritten)
        @printf("%26s %10.3f µs %10.3f µs %8.2f×   (first call: %6.3f s vs %6.3f s)\n",
            label, 1e6 * steady_generated, 1e6 * steady_handwritten,
            steady_generated / steady_handwritten, first_generated, first_handwritten)
    end

    println("""
  `f` is where the hand-written version wins: it computes each of the N(N-1)/2 distances once,
  whereas the generated code re-derives shared subexpressions across components even with CSE on.
  For `v` and `ϑ` the generated code is the faster one — they are linear in p resp. q̇, so there is
  nothing to share, and the generated form is fully unrolled straight-line code while the
  hand-written one loops over bodies and unpacks `params`.""")

    println()
    @printf("  L  generated = %+.17e\n", lagrangian_functions.L(0.0, q, v, parameters))
    @printf("  L  hand-written = %+.17e\n", OSS.lagrangian(0.0, q, v, parameters))
    @printf("  H  generated = %+.17e\n", hamiltonian_functions.H(0.0, q, p, parameters))
    @printf("  H  hand-written = %+.17e\n", OSS.hamiltonian(0.0, q, p, parameters))

    println("\n  Agreement (max absolute difference against the scale of the quantity):")
    for (label, gen_call, hand_call, sizes, args) in (
        ("f (Lagrangian)", lagrangian_functions.f, OSS.outer_solar_system_f, dof, (q, v)),
        ("v", hamiltonian_functions.v, OSS.outer_solar_system_v, dof, (q, p)),
        ("ϑ", lagrangian_functions.ϑ, OSS.outer_solar_system_ϑ, dof, (q, v)),
    )
        a = zeros(sizes)
        b = zeros(sizes)
        gen_call(a, 0.0, args..., parameters)
        hand_call(b, 0.0, args..., parameters)
        @printf("%20s   max|Δ| = %.3e   scale = %.3e\n", label, maximum(abs, a - b), maximum(abs, a))
    end
end

println()
