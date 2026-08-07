module LotkaVolterra2dPoincareInvariants

using PoincareInvariants

import GeometricProblems.LotkaVolterra2d
import GeometricProblems.LotkaVolterra2dGauge
import GeometricProblems.LotkaVolterra2dSingular
import GeometricProblems.LotkaVolterra2dSymmetric

# `lotka_volterra_2d_equations.jl` is included into all four Lotka-Volterra 2d modules, so all four
# declare these four names and all four need methods here. The bodies are identical up to the
# module they resolve their equations in, hence the loop.
#
# `ode_loop`/`iode_loop` wrap the parent module's `initial_conditions_loop` into a problem; they are
# unexported scaffolding for the invariants and are used by nothing else. `f_loop` and
# `initial_conditions_loop` themselves stay in `src/`, since they need nothing from this package.
#
# All four are **dead** and throw `UndefVarError` when called:
# `PoincareInvariant1st` is not defined by PoincareInvariants 0.5 — the first invariant is
# `FirstPoincareInvariant`/`FirstPI` there, forms are in-place (`form(out, t, z, p)`), and the entry
# point is `compute!(pinv, integrate(PIEnsembleProblem(prob, pinv, init), method))` — and
# `lotka_volterra_2d_ode`/`lotka_volterra_2d_iode` are now `odeproblem`/`iodeproblem`.
#
# The bodies are left unrepaired deliberately: the Poincaré-invariant support needs a makeover and
# the PoincareInvariants interface is still being tuned, so 0.5 is not the interface to target.
# What *is* exercised is that this file loads: `test/Project.toml` lists PoincareInvariants, so the
# extension is precompiled during the test run and `test/poincare_invariants_tests.jl` asserts that
# it resolved and attached its methods to all four modules.
for M in (:LotkaVolterra2d, :LotkaVolterra2dGauge, :LotkaVolterra2dSingular, :LotkaVolterra2dSymmetric)
    @eval begin
        function $M.ode_loop(n)
            $M.lotka_volterra_2d_ode($M.initial_conditions_loop(n))
        end

        function $M.iode_loop(n)
            $M.lotka_volterra_2d_iode($M.initial_conditions_loop(n))
        end

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
