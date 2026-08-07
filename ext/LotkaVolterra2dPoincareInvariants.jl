module LotkaVolterra2dPoincareInvariants

using PoincareInvariants

import GeometricProblems.LotkaVolterra2d
import GeometricProblems.LotkaVolterra2dGauge
import GeometricProblems.LotkaVolterra2dSingular
import GeometricProblems.LotkaVolterra2dSymmetric

# `lotka_volterra_2d_equations.jl` is included into all four Lotka-Volterra 2d modules, so all four
# declare these five names and all four need methods here. The bodies are identical up to the
# module they resolve their equations in, hence the loop.
#
# `initial_conditions_loop` samples the loop the parent module's `f_loop` parameterises, and
# `ode_loop`/`iode_loop` wrap it into a problem. All three are unexported scaffolding for the
# invariants and are used by nothing else.
#
# Everything except `initial_conditions_loop` is **dead** and throws `UndefVarError` when called:
# `PoincareInvariant1st` is not defined by PoincareInvariants 0.5 — the first invariant is
# `FirstPoincareInvariant`/`FirstPI` there, forms are in-place (`form(out, t, z, p)`), and the entry
# point is `compute!(pinv, integrate(PIEnsembleProblem(prob, pinv, init), method))` — and
# `lotka_volterra_2d_ode`/`lotka_volterra_2d_iode` are now `odeproblem`/`iodeproblem`.
#
# The bodies are left unrepaired deliberately: the Poincaré-invariant support needs a makeover and
# the PoincareInvariants interface is still being tuned, so 0.5 is not the interface to target.
# Nothing exercises this file yet — an extension is precompiled only in an environment that has its
# trigger package, and `test/Project.toml` does not list PoincareInvariants. Adding it there is what
# would put this code in front of CI.
for M in (:LotkaVolterra2d, :LotkaVolterra2dGauge, :LotkaVolterra2dSingular, :LotkaVolterra2dSymmetric)
    @eval begin
        function $M.initial_conditions_loop(n)
            q₀ = zeros(2, n)

            for i in axes(q₀, 2)
                q₀[:, i] .= $M.f_loop(i, n)
            end

            return q₀
        end

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
