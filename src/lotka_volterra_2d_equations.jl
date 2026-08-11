
using GeometricEquations
using GeometricSolutions

export odeproblem,  daeproblem,
       podeproblem, pdaeproblem,
       iodeproblem, idaeproblem,
       lodeproblem, ldaeproblem,
       hodeproblem, hdaeproblem,
       iodeproblem_dg, ldaeproblem_slrk, idaeproblem_spark

export poincare_invariant_1st,
       poincare_invariant_2nd

export ode_poincare_invariant_1st,
       iode_poincare_invariant_1st

export compute_energy_error

const Δt = 0.01
const nt = 1000
const DEFAULT_TIMESPAN = (0.0, Δt*nt)
const DEFAULT_TIMESTEP = Δt

default_parameters(::Type{T}=Float64) where {T} = (a₁=T(-1.0), a₂=T(-1.0), b₁=T(1.0), b₂=T(2.0))
const reference_solution = [2.576489958858641, 1.5388112243762107]

const t₀ = DEFAULT_TIMESPAN[begin]
const q₀ = [2.0, 1.0]
const v₀ = [v₁(0, q₀, default_parameters()), v₂(0, q₀, default_parameters())]


# Centre and semi-axes of the ellipse in phase space that `f_loop` traces, for the first Poincaré
# invariant. `f_surface` reads them too, so that the two parameterisations cannot drift apart, and
# `test/poincare_invariants_tests.jl` reads them to build the disc that the loop bounds.
const loop_centre = (1.0, 1.0)
const loop_radii = (0.2, 0.3)

function f_loop(s)
   x0, y0 = loop_centre
   rx, ry = loop_radii

   xs = x0 + rx*cos(2π*s)
   ys = y0 + ry*sin(2π*s)

   qs = [xs, ys]

   return qs
end

function f_loop(i, n)
   f_loop(i/n)
end

# A rectangular patch of phase space for the second Poincaré invariant, concentric with `f_loop` and
# spanning half its semi-axes in each direction: x ∈ 1 ± 0.1, y ∈ 1 ± 0.15. It lies inside the loop
# but is *not* the region the loop bounds, so its invariant is a different number from the loop's;
# `test/poincare_invariants_tests.jl` integrates over both, and over the disc as well.
# `PoincareInvariants.getpoints` samples it over the unit square, at the points its plan prescribes
# (Padua points for the Chebyshev plan, a regular grid for the finite-difference one).
function f_surface(s, t)
   x0, y0 = loop_centre
   rx, ry = loop_radii

   xs = x0 + rx*(s - 0.5)
   ys = y0 + ry*(t - 0.5)

   qs = [xs, ys]

   return qs
end

# Samples the loop `f_loop` parameterises at `n` equidistant points. Together with `f_loop` this is
# the live half of the Poincaré-invariant scaffolding, and it depends on nothing optional, so it
# stays here rather than moving into the extension with the invariants themselves.
function initial_conditions_loop(n)
   q₀ = zeros(2, n)

   for i in axes(q₀,2)
       q₀[:,i] .= f_loop(i, n)
   end

   return q₀
end


compute_energy_error(t, q, params) = compute_invariant_error(t, q, params, hamiltonian)


"Creates an ODE object for the Lotka-Volterra 2D model."
function odeproblem(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    ODEProblem(lotka_volterra_2d_v, timespan, timestep, q₀; parameters=parameters, invariants=(h=hamiltonian,))
end

"Creates a Hamiltonian ODE object for the Lotka-Volterra 2D model."
function hodeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    HODEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, hamiltonian, timespan, timestep, q₀, p₀;
                parameters=parameters)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model."
function iodeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    IODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_g, timespan, timestep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian,), v̄=lotka_volterra_2d_v)
end

"Creates a partitioned ODE object for the Lotka-Volterra 2D model."
function podeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    PODEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, timespan, timestep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian,))
end

"Creates a variational ODE object for the Lotka-Volterra 2D model."
function lodeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    LODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_g, lotka_volterra_2d_ω, lagrangian, timespan, timestep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian,), v̄=lotka_volterra_2d_v)
end

"Creates a DAE object for the Lotka-Volterra 2D model."
function daeproblem(q₀=vcat(q₀,v₀), λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    DAEProblem(lotka_volterra_2d_v_dae, lotka_volterra_2d_u_dae, lotka_volterra_2d_ϕ_dae, timespan, timestep, q₀, λ₀;
                parameters=parameters, invariants=(h=hamiltonian,))
end

"Creates a Hamiltonian DAE object for the Lotka-Volterra 2D model."
function hdaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    HDAEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ,
                hamiltonian, timespan, timestep, q₀, p₀, λ₀; parameters=parameters)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function idaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    IDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian,),
                v̄=lotka_volterra_2d_v)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function idaeproblem_spark(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    IDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a partitioned DAE object for the Lotka-Volterra 2D model."
function pdaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    PDAEProblem(lotka_volterra_2d_v_ham, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model."
function ldaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    LDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ_lode,
                lotka_volterra_2d_ω, lagrangian,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model for use with SLRK integrators."
function ldaeproblem_slrk(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    LDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ,
                lotka_volterra_2d_ω, lagrangian,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model for use with DG integrators."
function iodeproblem_dg(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters=default_parameters())
    IODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f, lotka_volterra_2d_g,
                timespan, timestep, q₀, p₀; parameters=parameters, invariants=(h=hamiltonian,), v̄=lotka_volterra_2d_v)
end


# The Poincaré invariants are implemented in the `LotkaVolterra2dPoincareInvariants` extension
# (loaded with PoincareInvariants), in the same shape as the `*Plots` extensions. `f_loop`,
# `f_surface` and `initial_conditions_loop` above stay here instead: they parameterise and sample
# the loop and the surface in phase space, need nothing optional to do it, and are read off the
# module by the extension.
#
@doc raw"""
    poincare_invariant_1st(N; DT = Float64, plan = PoincareInvariants.DEFAULT_FIRST_PLAN)

Set up the first Poincaré invariant
```math
I_1 (t) = \oint_{\gamma_t} \vartheta ,
```
the integral of this module's one-form ``\vartheta`` over a loop ``\gamma_t``, sampled at `N` points.
The Lagrangian is degenerate, so the momentum is not an independent coordinate but determined by
``p = \vartheta(q)``, and the loop lives in the two-dimensional configuration space alone.

This is implemented in the `LotkaVolterra2dPoincareInvariants` extension and becomes available once
PoincareInvariants is loaded. Sample a loop with `f_loop`, advect the sample points with
`PIEnsembleProblem` and evaluate with `compute!`:

```julia
pinv = poincare_invariant_1st(200)
prob = iodeproblem(; timespan = (0.0, 1E2), timestep = 1E-1)
sol  = integrate(PIEnsembleProblem(prob, pinv, f_loop), VPRKGauss(2))
I₁   = compute!(pinv, sol, parameters(prob))
```

Note that `I₁` is preserved by a variational integrator only to the order of the discretisation, not
exactly: the numerical solution satisfies ``p = \vartheta(q)`` only up to the truncation error.

The four gauges of the Lotka-Volterra 2d model differ by gauge transformations, which change
``\vartheta`` by an exact form. Since the integral of an exact form over a closed loop vanishes, all
four give the same ``I_1`` for the same loop.

See also `poincare_invariant_2nd`.
"""
function poincare_invariant_1st end

@doc raw"""
    poincare_invariant_2nd(N; DT = Float64, plan = PoincareInvariants.DEFAULT_SECOND_PLAN)

Set up the second Poincaré invariant
```math
I_2 (t) = \int_{\sigma_t} \omega ,
```
the integral of this module's two-form ``\omega`` over a surface ``\sigma_t``, sampled at `N` points.
As for the first invariant the surface lives in the two-dimensional configuration space alone. The
default plan samples at Padua points and rounds `N` up to the next Padua number, so the invariant
may use slightly more points than requested; `getpointnum` reports how many.

This is implemented in the `LotkaVolterra2dPoincareInvariants` extension and becomes available once
PoincareInvariants is loaded. It is used exactly like the first invariant, with `f_surface` in place
of `f_loop`. If the surface is the region a loop bounds, then ``I_2`` over the surface and ``I_1``
over the loop agree by Stokes' theorem — `f_surface` is *not* that region for `f_loop`, but lies
inside it.

Unlike ``\vartheta``, ``\omega = -d\vartheta`` is gauge invariant, so all four gauges of the
Lotka-Volterra 2d model share the same two-form.

See also `poincare_invariant_1st`.
"""
function poincare_invariant_2nd end


# The four names below are the pre-0.4 interface. They are **dead** and throw `UndefVarError` when
# called: the invariant needs `PoincareInvariant1st`, which PoincareInvariants 0.5 does not define,
# and the two `*_loop` wrappers call `lotka_volterra_2d_ode`/`lotka_volterra_2d_iode`, which are now
# `odeproblem`/`iodeproblem`. They are superseded by the two constructors above and kept only so
# that the export list stays backwards compatible within 0.8; the extension has the details.
function ode_loop end
function iode_loop end
function ode_poincare_invariant_1st end
function iode_poincare_invariant_1st end
