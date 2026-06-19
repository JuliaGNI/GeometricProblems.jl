using CairoMakie
using GeometricIntegrators
using GeometricProblems.Pendulum

# Generate a single PDF showing the last solution step or save one frame per time step (animation mode)
const animation = false

# Simulation parameters
const Δt = 0.1
const tspan_short = (0.0, 50.0)
const tspan_long = (0.0, 1000.0)


# Extract (t, q, p, h) arrays from an ODE solution where state = [angle, angular velocity]
function ode_arrays(sol)
    t = sol.t
    q = [x[1] for x in sol.q]
    p = [x[2] for x in sol.q]
    h = hamiltonian.(t, q, p)
    return t, q, p, h
end

# Extract (t, q, p, h) arrays from a PODE solution
function pode_arrays(sol)
    t = sol.t
    q = [x[1] for x in sol.q]
    p = [x[1] for x in sol.p]
    h = hamiltonian.(t, q, p)
    return t, q, p, h
end

# Save one PDF frame per time step (animation mode)
function save_solution_frames(t, q, p, h, dir, labels, hmod = nothing)
    mkpath(dir)
    println("Saving $(length(t)) frames to $(dir)/")
    for i in eachindex(t)
        fig = plot_solution(i, t, q, p, h, labels, hmod)
        save(joinpath(dir, "pendulum-$(lpad(string(i), 4, '0')).pdf"), fig)
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


## Explicit Euler (ODE)
sol_ee = integrate(odeproblem(timespan = tspan_short, timestep = Δt), ExplicitEuler())
t_ee, q_ee, p_ee, h_ee = ode_arrays(sol_ee)
save_solution(t_ee, q_ee, p_ee, h_ee, "pendulum-explicit-euler", labels_ode)


## Implicit Euler (ODE)
sol_ie = integrate(odeproblem(timespan = tspan_short, timestep = Δt), ImplicitEuler())
t_ie, q_ie, p_ie, h_ie = ode_arrays(sol_ie)
save_solution(t_ie, q_ie, p_ie, h_ie, "pendulum-implicit-euler", labels_ode)


## Symplectic Euler A (PODE)
sol_sea = integrate(podeproblem(timespan = tspan_long, timestep = Δt), SymplecticEulerA())
t_sea, q_sea, p_sea, h_sea = pode_arrays(sol_sea)
save_solution(t_sea, q_sea, p_sea, h_sea, "pendulum-symplectic-euler-a", labels_hamiltonian)


## Symplectic Euler B (PODE)
sol_seb = integrate(podeproblem(timespan = tspan_long, timestep = Δt), SymplecticEulerB())
t_seb, q_seb, p_seb, h_seb = pode_arrays(sol_seb)
save_solution(t_seb, q_seb, p_seb, h_seb, "pendulum-symplectic-euler-b", labels_hamiltonian)
