using GeometricIntegrators: integrate
using Logging: Warn, with_logger
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
# another module, this filter matches nothing and the assertion passes vacuously rather than
# failing.
function integrate_quietly(problem, method)
    logger = TestLogger(min_level = Warn)
    solution = with_logger(() -> integrate(problem, method), logger)
    # Assert on the distinct messages rather than their count, so that a failure reports what the
    # solver actually complained about.
    @test isempty(unique(l.message for l in logger.logs if nameof(l._module) === :SimpleSolvers))
    solution
end
