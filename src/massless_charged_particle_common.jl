
using GeometricEquations
using GeometricSolutions

export ϑ, ω, A, B, ϕ, E, hamiltonian, lagrangian
export odeproblem, iodeproblem, lodeproblem, idaeproblem, idaeproblem_spark, ldaeproblem
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

# the symplectic two-form (derivatives dϑᵢdxⱼ are defined by the including module).
# Convention (shared with LotkaVolterra2d/4d and PointVortices): Ωᵢⱼ = ∂ϑᵢ/∂qⱼ - ∂ϑⱼ/∂qᵢ, i.e.
# Ω = -dϑ. With it, the Euler-Lagrange equations of L = ϑ⋅v - H read Ω v = -∇H, and Ω₁₂ = -B(q)
# for both gauges.
ω₁₂(t, q, params) = dϑ₁dx₂(t, q, params) - dϑ₂dx₁(t, q, params)

function ω(t, q, params)
    Ω₁₂ = ω₁₂(t, q, params)
    Z = zero(Ω₁₂)
    [Z (+Ω₁₂); (-Ω₁₂) Z]
end

# components of the force (derivatives dϑᵢdxⱼ are defined by the including module)
f₁(v, t, q, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₂dx₁(t, q, params) * v[2]
f₂(v, t, q, params) = dϑ₁dx₂(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

g₁(v, t, q, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₁dx₂(t, q, params) * v[2]
g₂(v, t, q, params) = dϑ₂dx₁(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

# Hamiltonian (total energy)
hamiltonian(t, q, params) = ϕ(q, params)

# The Hamiltonian depends only on q; this method covers the iode/lode/idae/ldae contexts, where
# the invariant is evaluated with an additional velocity (or momentum) slot.
hamiltonian(t, q, p, params) = hamiltonian(t, q, params)

# Lagrangian L = ϑ(q)⋅v - H(q). It is linear in v and hence degenerate (∂²L/∂v² = 0).
function lagrangian(t, q, v, params)
    ϑ₁(t, q, params) * v[1] + ϑ₂(t, q, params) * v[2] - hamiltonian(t, q, params)
end

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

function massless_charged_particle_ω(Ω, t, q, params)
    Ω₁₂ = ω₁₂(t, q, params)
    Ω[1,1] = 0
    Ω[1,2] = + Ω₁₂
    Ω[2,1] = - Ω₁₂
    Ω[2,2] = 0
    nothing
end

# LODE/LDAE evaluate the symplectic matrix with an extra velocity slot; ω depends only on q.
massless_charged_particle_ω(Ω, t, q, v, params) = massless_charged_particle_ω(Ω, t, q, params)

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

# The secondary projection force of the LDAE. In contrast to g, which contracts the columns of ∇ϑ
# (f₁/f₂), ḡ contracts its rows (g₁/g₂); mirrors LotkaVolterra2d.
function massless_charged_particle_ḡ(g, t, q, λ, params)
    g[1] = g₁(λ, t, q, params)
    g[2] = g₂(λ, t, q, params)
    nothing
end

massless_charged_particle_ḡ(g, t, q, p, λ, params) = massless_charged_particle_ḡ(g, t, q, λ, params)
massless_charged_particle_ḡ(g, t, q, v, p, λ, params) = massless_charged_particle_ḡ(g, t, q, λ, params)

function massless_charged_particle_u(u, t, q, v, params)
    u .= v
    nothing
end

massless_charged_particle_u(u, t, q, p, v, params) = massless_charged_particle_u(u, t, q, v, params)
# IDAE projects along the multiplier λ (u .= λ), mirroring LotkaVolterra2d.
massless_charged_particle_u(u, t, q, v, p, λ, params) = massless_charged_particle_u(u, t, q, λ, params)

# The secondary projection of the LDAE, again along the multiplier λ (mirrors LotkaVolterra2d).
function massless_charged_particle_ū(u, t, q, λ, params)
    u .= λ
    nothing
end

massless_charged_particle_ū(u, t, q, p, λ, params) = massless_charged_particle_ū(u, t, q, λ, params)
massless_charged_particle_ū(u, t, q, v, p, λ, params) = massless_charged_particle_ū(u, t, q, λ, params)

function massless_charged_particle_ϕ(ϕ, t, q, p, params)
    ϕ[1] = p[1] - ϑ₁(t,q,params)
    ϕ[2] = p[2] - ϑ₂(t,q,params)
    nothing
end
# IDAE calls ϕ with an extra velocity slot; the constraint p = ϑ(q) does not depend on it.
massless_charged_particle_ϕ(ϕ, t, q, v, p, params) = massless_charged_particle_ϕ(ϕ, t, q, p, params)

# The secondary constraint ψ = ṗ - q̇⋅∇ϑ, obtained by differentiating ϕ = p - ϑ(q) in time.
function massless_charged_particle_ψ(ψ, t, q, p, q̇, ṗ, params)
    ψ[1] = ṗ[1] - g₁(q̇,t,q,params)
    ψ[2] = ṗ[2] - g₂(q̇,t,q,params)
    nothing
end

# LDAE evaluates ψ with an additional velocity slot (mirrors LotkaVolterra2d).
massless_charged_particle_ψ(ψ, t, q, v, p, q̇, ṗ, params) = massless_charged_particle_ψ(ψ, t, q, p, q̇, ṗ, params)



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

"Creates a variational (Lagrangian) ODE object for the massless charged particle in 2D."
function lodeproblem(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters = default_parameters())
    LODEProblem(massless_charged_particle_ϑ, massless_charged_particle_f,
            massless_charged_particle_g, massless_charged_particle_ω, lagrangian,
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

"Creates a variational DAE object for the massless charged particle in 2D."
function ldaeproblem(q₀=q₀; timespan=DEFAULT_TIMESPAN, timestep=DEFAULT_TIMESTEP, parameters = default_parameters())
    LDAEProblem(massless_charged_particle_ϑ, massless_charged_particle_f,
            massless_charged_particle_u, massless_charged_particle_g,
            massless_charged_particle_ϕ,
            massless_charged_particle_ū, massless_charged_particle_ḡ,
            massless_charged_particle_ψ,
            massless_charged_particle_ω, lagrangian,
            timespan, timestep, q₀, ϑ(0., q₀, parameters), zero(q₀);
            v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
            invariants=(h=hamiltonian,), parameters=parameters)
end


compute_energy_error(t, q, params) = compute_invariant_error(t, q, params, hamiltonian)

function compute_momentum_error(t, q::DataSeries{T}, p::DataSeries{T}, params::NamedTuple) where {T}
    err = DataSeries(zero(p[begin]), ntime(p))
    for i in axes(p, 1)
        err[i] = p[i] .- ϑ(t[i], q[i], params)
    end
    return err
end


export plot_solution, plot_phase_portrait, plot_traces

# Plot functions are implemented in the `MasslessChargedParticlePlots` extension (loaded with Makie).
function plot_solution end
function plot_phase_portrait end
function plot_traces end
