using CairoMakie
using GeometricIntegrators
using GeometricProblems.HarmonicOscillator

# Generate a single PDF showing the last solution step or save one frame per time step (animation mode)
const animation = false

# Simulation parameters
const Δt = 0.5
const tspan_long = (0.0, 1000.0)
const tspan_short_ee = (0.0, 50.0)
const tspan_short_ie = (0.0, 100.0)

const params = default_parameters()

# Shadow Hamiltonian for Symplectic Euler B
_H̃1(q, p) = -params.k * p * q / 2
_H̃2(q, p) = (params.k^2 * q^2 + p^2) / 12
H̃(q, p) = hamiltonian(0, q, p, params) - Δt * _H̃1(q, p) + Δt^2 * _H̃2(q, p)


# Extract (t, q, p, h) arrays from an ODE solution where state = [position, velocity]
function ode_arrays(sol)
    t = sol.t
    q = [x[1] for x in sol.q]
    p = [x[2] for x in sol.q]
    h = hamiltonian.(0, q, p, Ref(params))
    return t, q, p, h
end

# Extract (t, q, p, h) arrays from a PODE/HODE solution
function pode_arrays(sol)
    t = sol.t
    q = [x[1] for x in sol.q]
    p = [x[1] for x in sol.p]
    h = hamiltonian.(0, q, p, Ref(params))
    return t, q, p, h
end

# Save one PDF frame per time step (animation mode)
function save_solution_frames(t, q, p, h, dir, labels, hmod = nothing)
    mkpath(dir)
    println("Saving $(length(t)) frames to $(dir)/")
    for i in eachindex(t)
        fig = plot_solution(i, t, q, p, h, labels, hmod)
        save(joinpath(dir, "harmonic-oscillator-$(lpad(string(i), 4, '0')).pdf"), fig)
    end
end

# Save a single PDF showing the last solution step
function save_solution_plot(t, q, p, h, filename, labels, hmod = nothing)
    println("Saving plot to $(filename)")
    fig = plot_solution(lastindex(t), t, q, p, h, labels, hmod)
    save(filename, fig)
end

function save_solution(t, q, p, h, name, labels, hmod = nothing)
    if animation
        save_solution_frames(t, q, p, h, name, labels, hmod)
    else
        save_solution_plot(t, q, p, h, "$(name).pdf", labels, hmod)
    end
end


## Analytic solution — 10 full periods sampled at 100 points per period
T_period = 2π / params.ω
ana_prob = podeproblem(timespan = (0.0, T_period * 10), timestep = T_period / 100)
ana_sol = exact_solution(ana_prob)
t_ana, q_ana, p_ana, h_ana = pode_arrays(ana_sol)
save_solution(t_ana, q_ana, p_ana, h_ana, "harmonic-oscillator-analytic", labels_ode)


## Explicit Euler (ODE)
sol_ee = integrate(odeproblem(timespan = tspan_short_ee, timestep = 0.1), ExplicitEuler())
t_ee, q_ee, p_ee, h_ee = ode_arrays(sol_ee)
save_solution(t_ee, q_ee, p_ee, h_ee, "harmonic-oscillator-explicit-euler", labels_ode)


## Implicit Euler (ODE)
sol_ie = integrate(odeproblem(timespan = tspan_short_ie, timestep = 0.1), ImplicitEuler())
t_ie, q_ie, p_ie, h_ie = ode_arrays(sol_ie)
save_solution(t_ie, q_ie, p_ie, h_ie, "harmonic-oscillator-implicit-euler", labels_ode)


## Symplectic Euler A (PODE)
sol_sea = integrate(podeproblem(timespan = tspan_long, timestep = Δt), SymplecticEulerA())
t_sea, q_sea, p_sea, h_sea = pode_arrays(sol_sea)
save_solution(t_sea, q_sea, p_sea, h_sea, "harmonic-oscillator-symplectic-euler-a", labels_hamiltonian)


## Symplectic Euler B (PODE) — also plots the shadow Hamiltonian H̃
sol_seb = integrate(podeproblem(timespan = tspan_long, timestep = Δt), SymplecticEulerB())
t_seb, q_seb, p_seb, h_seb = pode_arrays(sol_seb)
save_solution(t_seb, q_seb, p_seb, h_seb, "harmonic-oscillator-symplectic-euler-b", labels_hamiltonian, H̃)
