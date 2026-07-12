using Test
using GeometricProblems
# GeometricEquations re-exports GeometricBase.timespan / .timestep; use it here because
# GeometricBase is not a direct test dependency.
using GeometricEquations: timespan, timestep

# Regression tests for issues #83 and #82.
#
# PR #81 mistakenly renamed the module-level default constants `tspan`/`tstep` to
# `timespan`/`timestep`, which shadowed the `GeometricBase.timespan` /
# `GeometricBase.timestep` functions inside the problem modules (issue #83). The
# fix renames those constants to `DEFAULT_TIMESPAN` / `DEFAULT_TIMESTEP`. These
# tests guard against a regression: in every problem module the names
# `timespan`/`timestep` must still resolve to the GeometricBase functions (never a
# shadowing constant), and the default values must live in `DEFAULT_TIMESPAN` /
# `DEFAULT_TIMESTEP`.

# All problem modules that define default-timespan/timestep constants.
const PROBLEM_MODULES = (
    :ABCFlow,
    :CoupledHarmonicOscillator,
    :DoublePendulum,
    :DuffingOscillator,
    :KuboOscillator,
    :LennardJonesOscillator,
    :LinearWave,
    :LorenzAttractor,
    :LotkaVolterra2d,
    :LotkaVolterra3d,
    :LotkaVolterra4d,
    :LotkaVolterra4dLagrangian,
    :MasslessChargedParticle,
    :MathewsLakshmananOscillator,
    :MorseOscillator,
    :PointVortices,
    :PointVorticesLinear,
    :RigidBody,
    :ThreeBody,
    :TodaLattice,
)

@testset "Default timespan/timestep constants ($name)" for name in PROBLEM_MODULES
    @test isdefined(GeometricProblems, name)
    M = getproperty(GeometricProblems, name)

    # The default timespan lives in DEFAULT_TIMESPAN and is a valid time interval.
    @test isdefined(M, :DEFAULT_TIMESPAN)
    @test getproperty(M, :DEFAULT_TIMESPAN) isa Tuple
    @test length(getproperty(M, :DEFAULT_TIMESPAN)) == 2

    # If a default timestep constant exists it must be named DEFAULT_TIMESTEP and be a number.
    if isdefined(M, :DEFAULT_TIMESTEP)
        @test getproperty(M, :DEFAULT_TIMESTEP) isa Number
    end

    # The names `timespan`/`timestep` must NOT be shadowed by a module constant:
    # where they resolve at all, they must be the GeometricBase functions.
    if isdefined(M, :timespan)
        @test getproperty(M, :timespan) === timespan
    end
    if isdefined(M, :timestep)
        @test getproperty(M, :timestep) === timestep
    end
end
