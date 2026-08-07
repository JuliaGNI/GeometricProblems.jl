module LotkaVolterra2dPoincareInvariants

using PoincareInvariants

import GeometricProblems.LotkaVolterra2d
import GeometricProblems.LotkaVolterra2dGauge
import GeometricProblems.LotkaVolterra2dSingular
import GeometricProblems.LotkaVolterra2dSymmetric

# `lotka_volterra_2d_equations.jl` is included into all four Lotka-Volterra 2d modules, so all four
# export `ode_poincare_invariant_1st`/`iode_poincare_invariant_1st` and all four need a method here.
# The bodies are identical up to the module they resolve their equations in, hence the loop.
#
# These methods are **dead**, and were dead under the `Requires.@require` block this replaces:
# `PoincareInvariant1st` was removed in PoincareInvariants 0.5.0, which moved the package from
# SciML to JuliaGNI and reworked the interface (`FirstPoincareInvariant`/`FirstPI`, in-place forms
# `form(out, t, z, p)`, and `compute!(pinv, integrate(PIEnsembleProblem(prob, pinv, init), method))`),
# and `lotka_volterra_2d_ode`/`lotka_volterra_2d_iode` are pre-0.7 constructor names — they are
# `odeproblem`/`iodeproblem` now. Calling either method therefore throws `UndefVarError`.
#
# They are carried over unchanged on purpose. The Poincaré-invariant support needs a complete
# makeover and the upstream interface is still being tuned, so 0.5 is not the target to write
# against and a port done now would only be rewritten. What the move buys is visibility: an
# extension is precompiled, so this file is at least loaded and parsed by CI, which is exactly what
# `Requires` deferred to load time and nothing ever triggered.
for M in (:LotkaVolterra2d, :LotkaVolterra2dGauge, :LotkaVolterra2dSingular, :LotkaVolterra2dSymmetric)
    @eval begin
        function $M.ode_poincare_invariant_1st(timestep, nloop, ntime, nsave, DT = Float64)
            PoincareInvariant1st($M.lotka_volterra_2d_ode, $M.f_loop, $M.ϑ,
                timestep, 2, nloop, ntime, nsave, DT)
        end

        function $M.iode_poincare_invariant_1st(timestep, nloop, ntime, nsave, DT = Float64)
            PoincareInvariant1st($M.lotka_volterra_2d_iode, $M.f_loop, $M.ϑ,
                timestep, 2, nloop, ntime, nsave, DT)
        end
    end
end

end
