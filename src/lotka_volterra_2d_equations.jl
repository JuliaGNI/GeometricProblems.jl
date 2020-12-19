
using GeometricIntegrators.Equations
using Requires

import ..Diagnostics: compute_invariant_error

export lotka_volterra_2d_ode,  lotka_volterra_2d_dae,
       lotka_volterra_2d_pode, lotka_volterra_2d_pdae,
       lotka_volterra_2d_iode, lotka_volterra_2d_idae,
       lotka_volterra_2d_vode, lotka_volterra_2d_vdae,
       lotka_volterra_2d_hode, lotka_volterra_2d_hdae,
       lotka_volterra_2d_dg,   lotka_volterra_2d_slrk,
       lotka_volterra_2d_idae_spark

export lotka_volterra_2d_ode_poincare_invariant_1st,
       lotka_volterra_2d_iode_poincare_invariant_1st

export compute_energy_error

const Δt = 0.1
const nt = 1000

const parameters = (a₁=-1.0, a₂=-1.0, b₁=1.0, b₂=2.0)
const reference_solution = [0.34520671417508214, 1.1140647339845946]

const q₀ = [2.0, 1.0]
const v₀=[v₁(0, q₀, parameters), v₂(0, q₀, parameters)]


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
function lotka_volterra_2d_ode(q₀=q₀; params=parameters)
    ODE(lotka_volterra_2d_v, q₀; parameters=params, h=hamiltonian)
end

"Creates a Hamiltonian ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_hode(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀); params=parameters)
    HODE(lotka_volterra_2d_v, lotka_volterra_2d_f, hamiltonian, q₀, p₀;
         parameters=params)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_iode(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀); params=parameters)
    IODE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
         lotka_volterra_2d_g, q₀, p₀;
         parameters=params, h=hamiltonian, v̄=lotka_volterra_2d_v)
end

"Creates a partitioned ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_pode(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀); params=parameters)
    PODE(lotka_volterra_2d_v, lotka_volterra_2d_f,
         q₀, p₀; parameters=params, h=hamiltonian)
end

"Creates a variational ODE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_vode(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀); params=parameters)
    VODE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
         lotka_volterra_2d_g, q₀, p₀;
         parameters=params, h=hamiltonian, v̄=lotka_volterra_2d_v,
         Ω=lotka_volterra_2d_ω, ∇H=lotka_volterra_2d_dH)
end

"Creates a DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_dae(q₀=vcat(q₀,v₀), λ₀=zero(q₀); params=parameters)
    DAE(lotka_volterra_2d_v_dae, lotka_volterra_2d_u_dae, lotka_volterra_2d_ϕ_dae, q₀, λ₀;
        parameters=params, h=hamiltonian)
end

"Creates a Hamiltonian DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_hdae(q₀=q₀, p₀=ϑ(0, q₀), λ₀=zero(q₀); params=parameters)
    HDAE(lotka_volterra_2d_v, lotka_volterra_2d_f,
         lotka_volterra_2d_u, lotka_volterra_2d_g,
         lotka_volterra_2d_u̅, lotka_volterra_2d_g̅,
         lotka_volterra_2d_ϕ, lotka_volterra_2d_ψ,
         hamiltonian, q₀, p₀, λ₀;
         parameters=params)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_idae(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀), λ₀=zero(q₀); params=parameters)
    IDAE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
         lotka_volterra_2d_u, lotka_volterra_2d_g,
         lotka_volterra_2d_ϕ, q₀, p₀, λ₀;
         parameters=params, h=hamiltonian, v̄=lotka_volterra_2d_v)
end

"Creates an implicit DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_idae_spark(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀), λ₀=zero(q₀); params=parameters)
    IDAE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
         lotka_volterra_2d_u, lotka_volterra_2d_g,
         lotka_volterra_2d_ϕ, q₀, p₀, λ₀;
         parameters=params, h=hamiltonian,
         v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a partitioned DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_pdae(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀), λ₀=zero(q₀); params=parameters)
    PDAE(lotka_volterra_2d_v_ham, lotka_volterra_2d_f_ham,
         lotka_volterra_2d_u, lotka_volterra_2d_g,
         lotka_volterra_2d_ϕ, q₀, p₀, λ₀;
         parameters=params, h=hamiltonian,
         v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model."
function lotka_volterra_2d_vdae(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀), λ₀=zero(q₀); params=parameters)
    VDAE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f_ham,
         lotka_volterra_2d_g, lotka_volterra_2d_g̅,
         lotka_volterra_2d_ϕ, lotka_volterra_2d_ψ,
         q₀, p₀, λ₀; parameters=params, h=hamiltonian,
         v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates a variational DAE object for the Lotka-Volterra 2D model for use with SLRK integrators."
function lotka_volterra_2d_slrk(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀), λ₀=zero(q₀); params=parameters)
    VDAE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
            lotka_volterra_2d_g, lotka_volterra_2d_g̅,
            lotka_volterra_2d_ϕ, lotka_volterra_2d_ψ,
            q₀, p₀, λ₀; parameters=params, h=hamiltonian,
            v̄=lotka_volterra_2d_v, f̄=lotka_volterra_2d_f)
end

"Creates an implicit ODE object for the Lotka-Volterra 2D model for use with DG integrators."
function lotka_volterra_2d_dg(q₀=q₀, p₀=lotka_volterra_2d_pᵢ(q₀); params=parameters)
    IODE(lotka_volterra_2d_ϑ, lotka_volterra_2d_f,
            lotka_volterra_2d_g, q₀, p₀;
            parameters=params, h=hamiltonian, v̄=lotka_volterra_2d_v)
end


function lotka_volterra_2d_loop_ode(n)
   lotka_volterra_2d_ode(initial_conditions_loop(n))
end

function lotka_volterra_2d_loop_iode(n)
   lotka_volterra_2d_iode(initial_conditions_loop(n))
end


function __init__()
    @require PoincareInvariants = "26663084-47d3-540f-bd97-40ca743aafa4" begin

        function lotka_volterra_2d_ode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
            PoincareInvariant1st(lotka_volterra_2d_ode, f_loop, ϑ, Δt, 2, nloop, ntime, nsave, DT)
        end

        function lotka_volterra_2d_iode_poincare_invariant_1st(Δt, nloop, ntime, nsave, DT=Float64)
            PoincareInvariant1st(lotka_volterra_2d_iode, f_loop, ϑ, Δt, 2, nloop, ntime, nsave, DT)
        end
        
    end
end
