using GeometricIntegrators: integrate
using Logging: Warn, current_logger, with_logger, handle_message, min_enabled_level, shouldlog
using Test

# Integrate, and assert that the nonlinear solver stayed quiet.
#
# An implicit step whose Newton iteration cannot reach a root still returns a solution rather than
# throwing — `ThreeBody.sympnets_initial_condition` is such a case, since it ends in a collision —
# so on its own it fails no assertion while filling the log with line-search warnings. Asserting
# silence turns that into a test failure; `verbosity = 0` would only hide it.
#
# Warnings are attributed by originating module rather than caught with `@test_nowarn`, both to
# ignore the unrelated `GeometricEquations` warning about single-sample ensembles and because
# SimpleSolvers is only an indirect dependency (via GeometricIntegrators), so it must not be
# imported. The `:SimpleSolvers` symbol is therefore load-bearing: should those messages move to
# another module, this filter would match nothing and all 33 assertions built on it would pass
# vacuously rather than failing. The positive control in `three_body_tests.jl` is what keeps that
# from going unnoticed — it integrates a case known to complain and asserts the filter still
# catches it.

const SIMPLESOLVERS = :SimpleSolvers

# Integrate while capturing the log, and return the distinct SimpleSolvers messages alongside the
# solution. Records from every other module are forwarded to the enclosing logger rather than
# swallowed, so that unrelated warnings stay visible and stay catchable by an enclosing
# `@test_nowarn`.
function integrate_capturing_messages(problem, method)
    logger = TestLogger(min_level = Warn)
    solution = with_logger(() -> integrate(problem, method), logger)

    outer = current_logger()
    for l in logger.logs
        nameof(l._module) === SIMPLESOLVERS && continue
        min_enabled_level(outer) <= l.level || continue
        shouldlog(outer, l.level, l._module, l.group, l.id) || continue
        handle_message(outer, l.level, l.message, l._module, l.group, l.id, l.file, l.line;
            l.kwargs...)
    end

    unique(l.message for l in logger.logs if nameof(l._module) === SIMPLESOLVERS), solution
end

"Distinct SimpleSolvers messages emitted while integrating `problem` with `method`."
simplesolvers_messages(problem, method) = first(integrate_capturing_messages(problem, method))

function integrate_quietly(problem, method)
    messages, solution = integrate_capturing_messages(problem, method)
    # Assert on the distinct messages rather than their count, so that a failure reports what the
    # solver actually complained about.
    @test isempty(messages)
    solution
end
