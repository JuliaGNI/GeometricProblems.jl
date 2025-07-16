
using GeometricEquations
using Requires

import ..Diagnostics: compute_invariant_error

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
const timespan = (0.0, Δt*nt)

const default_parameters = (a₁=-1.0, a₂=-1.0, b₁=1.0, b₂=2.0)
const reference_solution = [2.576489958858641, 1.5388112243762107]

const t₀ = timespan[begin]
const q₀ = [2.0, 1.0]
const v₀ = [v₁(0, q₀, default_parameters), v₂(0, q₀, default_parameters)]


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


compute_energy_error(t,q,params) = compute_invariant_error(t,q, (t,q) -> hamiltonian(t,q,params))


"Creates an ODE object for the Lotka-Volterra 2D model."
function odeproblem(q₀=q₀; timespan=timespan, timestep=Δt, parameters=default_parameters)
    ODEProblem(lotka_volterra_2d_v, timespan, timestep, q₀; parameters=parameters, invariants=(h=hamiltonian,))
end

"Creates a Hamiltonian ODE object for the Lotka-Volterra 2D model."
function hodeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    HODEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, hamiltonian_pode, timespan, timestep, q₀, p₀;
                parameters=parameters)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model."
function iodeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    IODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_g, timespan, timestep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_2d_v)
end

"Creates a partitioned ODE object for the Lotka-Volterra 2D model."
function podeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    PODEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, timespan, timestep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian_pode,))
end

"Creates a variational ODE object for the Lotka-Volterra 2D model."
function lodeproblem(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    LODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_g, lagrangian, lotka_volterra_2d_ω, timespan, timestep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_2d_v)
end

"Creates a DAE object for the Lotka-Volterra 2D model."
function daeproblem(q₀=vcat(q₀,v₀), λ₀=zero(q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    DAEProblem(lotka_volterra_2d_v_dae, lotka_volterra_2d_u_dae, lotka_volterra_2d_ϕ_dae, timespan, timestep, q₀, λ₀;
                parameters=parameters, invariants=(h=hamiltonian,))
end

"Creates a Hamiltonian DAE object for the Lotka-Volterra 2D model."
function hdaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    HDAEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, 
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ,
                hamiltonian_pode, timespan, timestep, q₀, p₀, λ₀; parameters=parameters)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function idaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    IDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function idaeproblem_spark(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    IDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a partitioned DAE object for the Lotka-Volterra 2D model."
function pdaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    PDAEProblem(lotka_volterra_2d_v_ham, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_pode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model."
function ldaeproblem(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    LDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ_lode,
                lagrangian, lotka_volterra_2d_ω,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model for use with SLRK integrators."
function ldaeproblem_slrk(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    LDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ,
                lagrangian, lotka_volterra_2d_ω,
                timespan, timestep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model for use with DG integrators."
function iodeproblem_dg(q₀=q₀, p₀=ϑ(t₀, q₀); timespan=timespan, timestep=Δt, parameters=default_parameters)
    IODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f, lotka_volterra_2d_g,
                timespan, timestep, q₀, p₀; parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_2d_v)
end


function ode_loop(n)
   lotka_volterra_2d_ode(initial_conditions_loop(n))
end

function iode_loop(n)
   lotka_volterra_2d_iode(initial_conditions_loop(n))
end


function __init__()
    @require PoincareInvariants = "26663084-47d3-540f-bd97-40ca743aafa4" begin

        function ode_poincare_invariant_1st(timestep, nloop, ntime, nsave, DT=Float64)
            PoincareInvariant1st(lotka_volterra_2d_ode, f_loop, ϑ, timestep, 2, nloop, ntime, nsave, DT)
        end

        function iode_poincare_invariant_1st(timestep, nloop, ntime, nsave, DT=Float64)
            PoincareInvariant1st(lotka_volterra_2d_iode, f_loop, ϑ, timestep, 2, nloop, ntime, nsave, DT)
        end
        
    end
end
