# Timing helpers shared by the scripts in this directory.
#
# These are `include`d rather than packaged: the benchmark scripts are run directly with
# `julia --project=…` and there is no `benchmark/Project.toml`, so they can only use what the main or
# the test project provides. BenchmarkTools is not a dependency anywhere in this repository, hence
# the hand-rolled `percall` below — the point of which is that the functions being measured here take
# a few nanoseconds, which is close enough to the resolution of `time()` that a fixed repetition
# count produces meaningless numbers.

# A sink for the results of the timed calls. Without consuming the output, a cheap generated
# function inlined into the timing loop can be optimised away entirely, which is what produced
# sub-nanosecond "0.0000 µs" readings and meaningless NaN/Inf ratios in an earlier version of
# `simplify_evaluation.jl`. `@noinline` keeps the accumulation itself from being elided.
const SINK = Ref(0.0)

@noinline sink!(x::Number) = (SINK[] += x; nothing)
@noinline sink!(x::AbstractArray) = (SINK[] += @inbounds x[begin]; nothing)
@noinline sink!(::Nothing) = nothing

"The wall-clock resolution these measurements are built on, estimated empirically."
function timer_resolution(samples=1000)
    smallest = Inf
    for _ in 1:samples
        a = time()
        b = time()
        d = b - a
        d > 0 && (smallest = min(smallest, d))
    end
    return smallest
end

const TIMER_RES = timer_resolution()

"""
    percall(f; target_seconds, samples) -> (seconds_per_call, reliable)

Seconds per call to `f`, the best of `samples` batches.

The number of repetitions per batch is calibrated so that each batch runs for at least
`target_seconds`, which keeps the measured interval far above the timer resolution — a single
evaluation of these functions is a few nanoseconds, so a fixed repetition count is not enough.
`reliable` is false when even the calibrated batch could not be made long enough to trust.

`f` must already have been called once, so that compilation is not measured.
"""
function percall(f; target_seconds=0.05, samples=7)
    # calibrate the repetition count
    reps = 64
    while reps < 1 << 30
        start = time()
        for _ in 1:reps
            sink!(f())
        end
        elapsed = time() - start
        elapsed >= target_seconds && break
        # grow by the shortfall, with headroom, rather than doubling blindly
        factor = elapsed <= 0 ? 16 : clamp(ceil(Int, 1.5 * target_seconds / elapsed), 2, 64)
        reps *= factor
    end

    best = Inf
    total = 0.0
    for _ in 1:samples
        GC.gc()
        start = time()
        for _ in 1:reps
            sink!(f())
        end
        elapsed = time() - start
        total = max(total, elapsed)
        best = min(best, elapsed / reps)
    end
    # a batch must span many timer ticks for the per-call figure to mean anything
    return best, total > 1000 * TIMER_RES
end

"""
    elapsed(f) -> (seconds, result)

Wall-clock time of a *single* call to `f`, together with what it returned.

For one-off costs — symbolic construction, or the first call to a generated function, which is
dominated by compiling it. Both are far above the timer resolution, so no repetition is needed; the
`GC.gc()` beforehand keeps a collection triggered by earlier work out of the measurement.
"""
function elapsed(f)
    GC.gc()
    start = time()
    result = f()
    return time() - start, result
end
