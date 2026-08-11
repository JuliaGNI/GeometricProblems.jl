module LotkaVolterra2dPoincareInvariants

using PoincareInvariants

import GeometricProblems.LotkaVolterra2d
import GeometricProblems.LotkaVolterra2dGauge
import GeometricProblems.LotkaVolterra2dSingular
import GeometricProblems.LotkaVolterra2dSymmetric

# `lotka_volterra_2d_equations.jl` is included into all four Lotka-Volterra 2d modules, so all four
# declare these names and all four need methods here. The bodies are identical up to the module
# they resolve their equations in, hence the loop.
#
# `poincare_invariant_1st`/`poincare_invariant_2nd` are the live interface. They build the
# invariants of PoincareInvariants 0.5 over the module's own one- and two-form: both are already
# in-place with the `form(out, t, z, p)` argument order the package calls, so they are passed
# through unwrapped. `D = 2`: these are *degenerate* Lagrangian systems, whose loop and surface
# live in the two-dimensional configuration space alone, with the momentum determined by ϑ(q).
#
# The point count `N` is the number of phase space samples. The first invariant's Fourier plan
# takes any number of them and wants a periodic parameterisation, which `f_loop` is; the second
# invariant's Chebyshev plan samples at Padua points and rounds `N` up to the next Padua number.
#
# `f_loop`, `f_surface` and `initial_conditions_loop` stay in `src/`, since they need nothing from
# this package.
#
# `ode_loop`/`iode_loop` and the two `*_poincare_invariant_1st` below are the pre-0.4 interface,
# superseded by the two constructors above. All four are **dead** and throw `UndefVarError` when
# called: `PoincareInvariant1st` is not defined by PoincareInvariants 0.5, and
# `lotka_volterra_2d_ode`/`lotka_volterra_2d_iode` are now `odeproblem`/`iodeproblem`. They are
# left unrepaired: removing them would be breaking, and nothing calls them.
# `test/poincare_invariants_tests.jl` pins both halves — that the new constructors work, and that
# the old names still throw.
for M in (:LotkaVolterra2d, :LotkaVolterra2dGauge, :LotkaVolterra2dSingular, :LotkaVolterra2dSymmetric)
    @eval begin
        function $M.poincare_invariant_1st(N; DT = Float64, plan = FirstFourierPlan)
            FirstPI{DT, 2}($M.lotka_volterra_2d_ϑ, N, plan)
        end

        function $M.poincare_invariant_2nd(N; DT = Float64, plan = SecondChebyshevPlan)
            SecondPI{DT, 2}($M.lotka_volterra_2d_ω, N, plan)
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
