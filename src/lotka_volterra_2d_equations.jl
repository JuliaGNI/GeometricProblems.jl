
using GeometricEquations
using GeometricSolutions

export odeproblem,  daeproblem,
       podeproblem, pdaeproblem,
       iodeproblem, idaeproblem,
       lodeproblem, ldaeproblem,
       hodeproblem, hdaeproblem,
       iodeproblem_dg, ldaeproblem_slrk, idaeproblem_spark

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


function f_loop(s)
   rx = 0.2
   ry = 0.3
   x0 = 1.0
   y0 = 1.0

   xs = x0 + rx*cos(2π*s)
   ys = y0 + ry*sin(2π*s)

   qs = [xs, ys]

   return qs
end

function f_loop(i, n)
   f_loop(i/n)
end

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


function ode_loop(n)
   lotka_volterra_2d_ode(initial_conditions_loop(n))
end

function iode_loop(n)
   lotka_volterra_2d_iode(initial_conditions_loop(n))
end


# The first Poincaré invariant is implemented in the `LotkaVolterra2dPoincareInvariants` extension
# (loaded with PoincareInvariants), in the same shape as the `*Plots` extensions.
#
# Both methods are **dead** and were dead before the move: their bodies call `PoincareInvariant1st`,
# removed in PoincareInvariants 0.5.0, on `lotka_volterra_2d_ode`/`lotka_volterra_2d_iode`, which
# are pre-0.7 constructor names (`odeproblem`/`iodeproblem` today). They are carried over unchanged
# rather than repaired: the whole Poincaré-invariant integration is slated for a makeover, and the
# PoincareInvariants interface is still being tuned, so 0.5 is not the target to write against.
# Moving them out of `Requires.@require` and into an extension at least puts them where
# precompilation — and therefore CI — can see them.
function ode_poincare_invariant_1st end
function iode_poincare_invariant_1st end
