module MasslessChargedParticlePoincareInvariants

using PoincareInvariants

import GeometricProblems.MasslessChargedParticle
import GeometricProblems.MasslessChargedParticleSingular

# `massless_charged_particle_common.jl` is included into both gauges of the massless charged
# particle, so both declare these names and both need methods here. The bodies are identical up to
# the module they resolve their equations in, hence the loop.
#
# The invariants are built over the module's own one- and two-form: both are already in-place with
# the `form(out, t, z, p)` argument order PoincareInvariants 0.5 calls, so they are passed through
# unwrapped. `D = 2`: the Lagrangian is degenerate, so the loop and the surface live in the
# two-dimensional configuration space alone, with the momentum determined by ϑ(q).
#
# The point count `N` is the number of phase space samples. The `plan` keywords default to
# PoincareInvariants' own defaults rather than naming a plan, so that upstream stays free to change
# them: currently `FirstFourierPlan`, which takes any number of points and wants a periodic
# parameterisation (which `f_loop` is), and `SecondChebyshevPlan`, which samples at Padua points and
# rounds `N` up to the next Padua number.
#
# `f_loop` and `f_surface` stay in `src/`, since they need nothing from this package.
for M in (:MasslessChargedParticle, :MasslessChargedParticleSingular)
    @eval begin
        function $M.poincare_invariant_1st(N; DT = Float64, plan = PoincareInvariants.DEFAULT_FIRST_PLAN)
            FirstPI{DT, 2}($M.massless_charged_particle_ϑ, N, plan)
        end

        function $M.poincare_invariant_2nd(N; DT = Float64, plan = PoincareInvariants.DEFAULT_SECOND_PLAN)
            SecondPI{DT, 2}($M.massless_charged_particle_ω, N, plan)
        end
    end
end

end
