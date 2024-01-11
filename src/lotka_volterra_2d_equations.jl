
using GeometricEquations
using Requires

import ..Diagnostics: compute_invariant_error

export lotka_volterra_2d_ode,  lotka_volterra_2d_dae,
       lotka_volterra_2d_pode, lotka_volterra_2d_pdae,
       lotka_volterra_2d_iode, lotka_volterra_2d_idae,
       lotka_volterra_2d_lode, lotka_volterra_2d_ldae,
       lotka_volterra_2d_hode, lotka_volterra_2d_hdae,
       lotka_volterra_2d_dg,   lotka_volterra_2d_slrk,
       lotka_volterra_2d_idae_spark

export lotka_volterra_2d_ode_poincare_invariant_1st,
       lotka_volterra_2d_iode_poincare_invariant_1st

export compute_energy_error

const Δt = 0.01
const nt = 1000
const tspan = (0.0, Δt*nt)

const default_parameters = (a₁=-1.0, a₂=-1.0, b₁=1.0, b₂=2.0)
const reference_solution = [2.576489958858641, 1.5388112243762107]

const t₀ = tspan[begin]
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
function lotka_volterra_2d_ode(q₀=q₀; tspan=tspan, tstep=Δt, parameters=default_parameters)
    ODEProblem(lotka_volterra_2d_v, tspan, tstep, q₀; parameters=parameters, invariants=(h=hamiltonian,))
end

"Creates a Hamiltonian ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_hode(q₀=q₀, p₀=ϑ(t₀, q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    HODEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, hamiltonian_pode, tspan, tstep, q₀, p₀;
                parameters=parameters)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_iode(q₀=q₀, p₀=ϑ(t₀, q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    IODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_g, tspan, tstep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_2d_v)
end

"Creates a partitioned ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_pode(q₀=q₀, p₀=ϑ(t₀, q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    PODEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, tspan, tstep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian_pode,))
end

"Creates a variational ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_lode(q₀=q₀, p₀=ϑ(t₀, q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    LODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_g, lagrangian, lotka_volterra_2d_ω, tspan, tstep, q₀, p₀;
                parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_2d_v)
end

"Creates a DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_dae(q₀=vcat(q₀,v₀), λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    DAEProblem(lotka_volterra_2d_v_dae, lotka_volterra_2d_u_dae, lotka_volterra_2d_ϕ_dae, tspan, tstep, q₀, λ₀;
                parameters=parameters, invariants=(h=hamiltonian,))
end

"Creates a Hamiltonian DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_hdae(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    HDAEProblem(lotka_volterra_2d_v, lotka_volterra_2d_f, 
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ,
                hamiltonian_pode, tspan, tstep, q₀, p₀, λ₀; parameters=parameters)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_idae(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    IDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                tspan, tstep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_idae_spark(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    IDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                tspan, tstep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a partitioned DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_pdae(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    PDAEProblem(lotka_volterra_2d_v_ham, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                tspan, tstep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_pode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_ldae(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    LDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ_lode,
                lagrangian, lotka_volterra_2d_ω,
                tspan, tstep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model for use with SLRK integrators."
function lotka_volterra_2d_slrk(q₀=q₀, p₀=ϑ(t₀, q₀), λ₀=zero(q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    LDAEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
                lotka_volterra_2d_u, lotka_volterra_2d_g, lotka_volterra_2d_ϕ,
                lotka_volterra_2d_ū, lotka_volterra_2d_ḡ, lotka_volterra_2d_ψ,
                lagrangian, lotka_volterra_2d_ω,
                tspan, tstep, q₀, p₀, λ₀; parameters=parameters, invariants=(h=hamiltonian_iode,),
                v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model for use with DG integrators."
function lotka_volterra_2d_dg(q₀=q₀, p₀=ϑ(t₀, q₀); tspan=tspan, tstep=Δt, parameters=default_parameters)
    IODEProblem(lotka_volterra_2d_ϑ, lotka_volterra_2d_f, lotka_volterra_2d_g,
                tspan, tstep, q₀, p₀; parameters=parameters, invariants=(h=hamiltonian_iode,), v̄=lotka_volterra_2d_v)
end


function lotka_volterra_2d_loop_ode(n)
   lotka_volterra_2d_ode(initial_conditions_loop(n))
end

function lotka_volterra_2d_loop_iode(n)
   lotka_volterra_2d_iode(initial_conditions_loop(n))
end


function __init__()
    @require PoincareInvariants = "26663084-47d3-540f-bd97-40ca743aafa4" begin

        function lotka_volterra_2d_ode_poincare_invariant_1st(tstep, nloop, ntime, nsave, DT=Float64)
            PoincareInvariant1st(lotka_volterra_2d_ode, f_loop, ϑ, tstep, 2, nloop, ntime, nsave, DT)
        end

        function lotka_volterra_2d_iode_poincare_invariant_1st(tstep, nloop, ntime, nsave, DT=Float64)
            PoincareInvariant1st(lotka_volterra_2d_iode, f_loop, ϑ, tstep, 2, nloop, ntime, nsave, DT)
        end
        
    end
end
