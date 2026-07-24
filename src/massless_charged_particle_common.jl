
using GeometricEquations
using GeometricSolutions

export ϑ, A, B, ϕ, E, hamiltonian
export odeproblem, iodeproblem, idaeproblem, idaeproblem_spark
export compute_energy_error, compute_momentum_error


# default simulation parameters
const Δt = 0.2
const nt = 5000
const DEFAULT_TIMESPAN = (0.0, Δt*nt)
const DEFAULT_TIMESTEP = Δt

# default initial conditions and parameters
const q₀ = [1.0, 1.0]

default_parameters(::Type{T}=Float64) where {T} = (A₀ = T(1.0), E₀ = T(1.0))

# vector potential (components A₁, A₂ are defined by the including module)
A(q, params) = [A₁(q, params), A₂(q, params)]

# electrostatic potential
ϕ(q, params) = params[:E₀] * (cos(q[1]) + sin(q[2]))

# components of the electric field
E₁(q, params) = + params[:E₀] * sin(q[1])
E₂(q, params) = - params[:E₀] * cos(q[2])

E(q, params) = [E₁(q, params), E₂(q, params)]

# components of the velocity (B is defined by the including module)
v₁(t, q, params) = + E₂(q, params) / B(q, params)
v₂(t, q, params) = - E₁(q, params) / B(q, params)

# components of the one-form (symplectic potential)
ϑ₁(t, q, params) = A₁(q, params)
ϑ₂(t, q, params) = A₂(q, params)

ϑ(t, q, params) = [ϑ₁(t, q, params), ϑ₂(t, q, params)]

function ϑ(t, q::AbstractVector, params::NamedTuple, k::Int)
    if k == 1
        ϑ₁(t, q, params)
    elseif k == 2
        ϑ₂(t, q, params)
    else
        throw(BoundsError(ϑ,k))
    end
end

# components of the force (derivatives dϑᵢdxⱼ are defined by the including module)
f₁(v, t, q, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₂dx₁(t, q, params) * v[2]
f₂(v, t, q, params) = dϑ₁dx₂(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

g₁(v, t, q, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₁dx₂(t, q, params) * v[2]
g₂(v, t, q, params) = dϑ₂dx₁(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

# Hamiltonian (total energy)
hamiltonian(t, q, params) = ϕ(q, params)

# components of the gradient of the Hamiltonian
dHd₁(t, q, params) = - E₁(q, params)
dHd₂(t, q, params) = - E₂(q, params)


function massless_charged_particle_dH(dH, t, q, params)
    dH[1] = dHd₁(t, q, params)
    dH[2] = dHd₂(t, q, params)
    nothing
end

function massless_charged_particle_v(v, t, q, params)
    v[1] = v₁(t, q, params)
    v[2] = v₂(t, q, params)
    nothing
end

function massless_charged_particle_v(v, t, q, p, params)
    massless_charged_particle_v(v, t, q, params)
end

function massless_charged_particle_ϑ(Θ, t, q, params)
    Θ[1] = ϑ₁(t, q, params)
    Θ[2] = ϑ₂(t, q, params)
    nothing
end

massless_charged_particle_ϑ(Θ, t, q, v, params) = massless_charged_particle_ϑ(Θ, t, q, params)

function massless_charged_particle_f(f, t, q, v, params)
    f[1] = f₁(v, t, q, params) - dHd₁(t, q, params)
    f[2] = f₂(v, t, q, params) - dHd₂(t, q, params)
    nothing
end

function massless_charged_particle_f̄(f, t, q, v, params)
    f[1] = - dHd₁(t, q, params)
    f[2] = - dHd₂(t, q, params)
    nothing
end

function massless_charged_particle_g(g, t, q, v, params)
    g[1] = f₁(v, t, q, params)
    g[2] = f₂(v, t, q, params)
    nothing
end

massless_charged_particle_g(g, t, q, p, v, params) = massless_charged_particle_g(g, t, q, v, params)
# IDAE evaluates the constraint force g with the multiplier λ (mirrors LotkaVolterra2d).
massless_charged_particle_g(g, t, q, v, p, λ, params) = massless_charged_particle_g(g, t, q, λ, params)

function massless_charged_particle_u(u, t, q, v, params)
    u .= v
    nothing
end

massless_charged_particle_u(u, t, q, p, v, params) = massless_charged_particle_u(u, t, q, v, params)
# IDAE projects along the multiplier λ (u .= λ), mirroring LotkaVolterra2d.
massless_charged_particle_u(u, t, q, v, p, λ, params) = massless_charged_particle_u(u, t, q, λ, params)

function massless_charged_particle_ϕ(ϕ, t, q, p, params)
    ϕ[1] = p[1] - ϑ₁(t,q,params)
    ϕ[2] = p[2] - ϑ₂(t,q,params)
    nothing
end
# IDAE calls ϕ with an extra velocity slot; the constraint p = ϑ(q) does not depend on it.
massless_charged_particle_ϕ(ϕ, t, q, v, p, params) = massless_charged_particle_ϕ(ϕ, t, q, p, params)

function massless_charged_particle_ψ(ψ, t, q, p, v, f, params)
    ψ[1] = f[1] - g₁(v,t,q,params)
    ψ[2] = f[2] - g₂(v,t,q,params)
    nothing
end



"Creates an ODE object for the massless charged particle in 2D."
function odeproblem(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters = default_parameters())
    ODEProblem(massless_charged_particle_v, timespan, timestep, q₀; invariants=(h=hamiltonian,), parameters=parameters)
end

"Creates an implicit ODE object for the massless charged particle in 2D."
function iodeproblem(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters = default_parameters())
    IODEProblem(massless_charged_particle_ϑ, massless_charged_particle_f,
            massless_charged_particle_g,
            timespan, timestep, q₀, ϑ(0., q₀, parameters);
            v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
            invariants=(h=hamiltonian,), parameters=parameters)
end

"Creates an implicit DAE object for the massless charged particle in 2D."
function idaeproblem(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters = default_parameters())
    IDAEProblem(massless_charged_particle_ϑ, massless_charged_particle_f,
            massless_charged_particle_u, massless_charged_particle_g,
            massless_charged_particle_ϕ,
            timespan, timestep, q₀, ϑ(0., q₀, parameters), zero(q₀);
            v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
            invariants=(h=hamiltonian,), parameters=parameters)
end

"Creates an implicit DAE object for the massless charged particle in 2D."
function idaeproblem_spark(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters = default_parameters())
    IDAEProblem(massless_charged_particle_ϑ, massless_charged_particle_f̄,
            massless_charged_particle_u, massless_charged_particle_g,
            massless_charged_particle_ϕ,
            timespan, timestep, q₀, ϑ(0., q₀, parameters), zero(q₀);
            v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
            invariants=(h=hamiltonian,), parameters=parameters)
end


compute_energy_error(t,q,params) = compute_invariant_error(t,q, (t,q) -> hamiltonian(t,q,params))
compute_momentum_error(t,q,p,params::NamedTuple) = compute_momentum_error(t, q, p, (t,q,k) -> ϑ(t,q,params,k))


export plot_solution, plot_phase_portrait, plot_traces

# Plot functions are implemented in the `MasslessChargedParticlePlots` extension (loaded with Makie).
function plot_solution end
function plot_phase_portrait end
function plot_traces end
