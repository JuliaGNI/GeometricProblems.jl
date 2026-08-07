using GeometricIntegrators: integrate
using Logging: Warn, with_logger
using Test

# Integrate and assert that the nonlinear solver stayed quiet. An implicit step whose Newton
# iteration cannot reach a root — as every `ThreeBody` step used to do when the default initial
# condition was still a member of the collisional `initial_conditions` grid, see
# `ThreeBody.sympnets_initial_condition` — still returns a solution, so it fails no assertion here
# while flooding the log with line-search warnings.
#
# Warnings are attributed by originating module rather than caught with `@test_nowarn`, both to
# ignore the unrelated `GeometricEquations` warning about single-sample ensembles and because
# SimpleSolvers is only an indirect dependency (via GeometricIntegrators), so it must not be
# imported. The `:SimpleSolvers` symbol is therefore load-bearing: should those messages ever move
# to another module, this filter matches nothing and the assertion passes vacuously rather than
# failing. `test/runtests.jl`'s suite-wide count is the backstop for that.
#
# This asserts silence rather than suppressing it, which is only possible from SimpleSolvers 0.10
# on: before that a line search built its own `Options` and never saw the solver's `verbosity`, so
# a converging solve reported its round-off floor as a warning and which orbits tripped it depended
# on platform floating-point details.
function integrate_quietly(problem, method)
    logger = TestLogger(min_level = Warn)
    solution = with_logger(() -> integrate(problem, method), logger)
    # Assert on the distinct messages rather than their count, so that a failure reports what the
    # solver actually complained about.
    @test isempty(unique(l.message for l in logger.logs if nameof(l._module) === :SimpleSolvers))
    solution
end
