using Aqua
using GeometricProblems
using Test

# Package-level quality assurance: method ambiguities, unbound type parameters, undefined exports,
# the agreement between `Project.toml` and `test/Project.toml`, stale dependencies, `[compat]`
# bounds, type piracy and persistent tasks. The rest of the suite exercises one problem module at a
# time and cannot see these, as each is a property of the package as a whole.
#
# On Julia 1.11 the persistent-tasks wrapper's `Pkg.precompile` also builds the Makie and
# PoincareInvariants extensions after the package loads. That took 218 s and 92 s in two runs,
# more than Aqua's default `tmax` of 30 s. A real persistent task blocks forever, so a larger
# `tmax` hides none.
Aqua.test_all(
    GeometricProblems;
    stale_deps = false,                    # issue #116: run as @test_broken below
    persistent_tasks = (; tmax = 300)
)

# `test_stale_deps` takes no `broken` keyword, so its check is run here directly.
# issue #116: Documenter is in [deps] of Project.toml and nothing under src/ or ext/ loads it
@testset "Stale dependencies" begin
    @test_broken isempty(Aqua.find_stale_deps(Base.PkgId(GeometricProblems)))  # issue #116
end
