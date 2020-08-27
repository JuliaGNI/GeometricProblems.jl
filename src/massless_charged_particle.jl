@doc raw"""
# Massless charged particle in 2D

The Lagrangian is given by
```math
L(x, \dot{x}) = A(x) \cdot \dot{x} - \phi (x) ,
```
with magnetic vector potential
```math
A(x) = \frac{A_0}{2} \big( 1 + x_1^2 + x_2^2 \big) \begin{pmatrix}
- x_2 \\
+ x_1 \\
\end{pmatrix} ,
```
electrostatic potential
```math
\phi(x) =  E_0 \, \big( \cos (x_1) + \sin(x_2) \big) ,
```
and magnetic and electric fields
```math
\begin{aligned}
B(x) &= \nabla \times A(x) = A_0 \, (1 + 2 x_1^2 + 2 x_2^2) , \\
E(x) &= - \nabla \phi(x) = E_0 \, \big( \sin x_1, \, - \cos x_2 \big)^T .
\end{aligned}
```

The Hamiltonian form of the equations of motion reads
```math
\dot{x} = \frac{1}{B(x)} \begin{pmatrix}
\hphantom{-} 0 & + 1 \\
- 1 & \hphantom{+} 0 \\
\end{pmatrix} \nabla \phi (x) .
```

"""
module MasslessChargedParticle

    using ..Diagnostics
    using GeometricIntegrators.Equations
    using GeometricIntegrators.Solutions
    
    using LaTeXStrings
    using Plots
    using Plots.PlotMeasures
    # using RecipesBase

    export ϑ, A, B, ϕ, E, hamiltonian
    export massless_charged_particle_ode, massless_charged_particle_iode
    export compute_energy_error

    # default simulation parameters
    Δt = 0.2
    nt = 5000

    # default initial conditions and parameters
    q₀ = [1.0, 1.0]
    parameters = (A₀ = 1.0, E₀ = 1.0)

    # components of the vector potential
    A₁(q, params) = - params[:A₀] * q[2] * (1 + q[1]^2 + q[2]^2) / 2
    A₂(q, params) = + params[:A₀] * q[1] * (1 + q[1]^2 + q[2]^2) / 2

    A(q, params) = [A₁(q, params), A₂(q, params)]

    # z-componend of the magnetic field
    B(q, params) = params[:A₀] * (1 + 2 * q[1]^2 + 2 * q[2]^2)

    # electrostatic potential
    ϕ(q, params) = params[:E₀] * (cos(q[1]) + sin(q[2]))
    # ϕ(q, params) = E₀ * (q[1]^2 + q[2]^2)

    # components of the electric field
    E₁(q, params) = + params[:E₀] * sin(q[1])
    E₂(q, params) = - params[:E₀] * cos(q[2])
    # E₁(q, params) = + params[:E₀] * q[1]
    # E₂(q, params) = + params[:E₀] * q[2]

    E(q, params) = [E₁(q, params), E₂(q, params)]

    # components of the velocity
    v₁(t, q, params) = - E₂(q, params) / B(q, params)
    v₂(t, q, params) = + E₁(q, params) / B(q, params)

    # components of the one-form (symplectic potential)
    ϑ₁(t, q, params) = A₁(q, params)
    ϑ₂(t, q, params) = A₂(q, params)

    ϑ(t, q, params) = [ϑ₁(t, q, params), ϑ₂(t, q, params)]

    # derivatives of the one-form components
    dϑ₁dx₁(t, q, params) = - params[:A₀] * q[1] * q[2]
    dϑ₁dx₂(t, q, params) = - params[:A₀] * (1 + q[1]^2 + 3 * q[2]^2) / 2
    dϑ₂dx₁(t, q, params) = + params[:A₀] * (1 + 3 * q[1]^2 + q[2]^2) / 2
    dϑ₂dx₂(t, q, params) = + params[:A₀] * q[1] * q[2]

    # components of the force
    f₁(t, q, v, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₂dx₁(t, q, params) * v[2]
    f₂(t, q, v, params) = dϑ₁dx₂(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

    # Hamiltonian (total energy)
    hamiltonian(t, q, params) = ϕ(q, params)

    # components of the gradient of the Hamiltonian
    dHd₁(t, q, params) = - E₁(q, params)
    dHd₂(t, q, params) = - E₂(q, params)


    function massless_charged_particle_dH(t, q, dH, params)
        dH[1] = dHd₁(t, q, params)
        dH[2] = dHd₂(t, q, params)
        nothing
    end 

    function massless_charged_particle_v(t, q, v, params)
        v[1] = v₁(t, q, params)
        v[2] = v₂(t, q, params)
        nothing
    end

    function massless_charged_particle_ϑ(t, q, Θ, params)
        Θ[1] = ϑ₁(t, q, params)
        Θ[2] = ϑ₂(t, q, params)
        nothing 
    end

    massless_charged_particle_ϑ(t, q, v, Θ, params) = massless_charged_particle_ϑ(t, q, Θ, params)

    function massless_charged_particle_f(t, q, v, f, params)
        f[1] = f₁(t, q, v, params) - dHd₁(t, q, params)
        f[2] = f₂(t, q, v, params) - dHd₂(t, q, params)
        nothing
    end

    function massless_charged_particle_g(t, q, v, g, params)
        g[1] = f₁(t, q, v, params)
        g[2] = f₂(t, q, v, params)
        nothing
    end


    "Creates an ODE object for the massless charged particle in 2D."
    function massless_charged_particle_ode(q₀=q₀; params=parameters)
        ODE(massless_charged_particle_v, q₀; h=hamiltonian, parameters=params)
    end

    "Creates an implicit ODE object for the massless charged particle in 2D."
    function massless_charged_particle_iode(q₀=q₀; params=parameters)
        IODE(massless_charged_particle_ϑ, massless_charged_particle_f,
                massless_charged_particle_g, q₀, ϑ(0., q₀, params);
                v=massless_charged_particle_v, h=hamiltonian, parameters=params)
    end


    compute_energy_error(t,q,params) = compute_invariant_error(t,q, (t,q) -> hamiltonian(t,q,params))


    """
    Plots the solution of a massless charged particle together with the energy error.

    Arguments:
    * `sol <: Solution`
    * `params <: NamedTuple`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `xlims=:auto`: xlims for solution plot
    * `ylims=:auto`: ylims for solution plot
    * `latex=true`: use LaTeX guides
    """
    @userplot Plot_Massless_Charged_Particle
    @recipe function f(p::Plot_Massless_Charged_Particle; nplot=1, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        params = p.args[2]

        H, ΔH = compute_energy_error(sol.t, sol.q, params);

        size   := (800,300)
        layout := @layout [solPlot{0.4w,1.0h} EPlot]
        legend := :none

        right_margin := 10mm
        bottom_margin := 10mm

        guidefont := font(14)
        tickfont  := font(10)

        @series begin
            subplot := 1

            if sol.nt ≤ 200
                markersize := 5
            else
                markersize  := 1
                markercolor := 1
                linecolor   := 1
                markerstrokewidth := 1
                markerstrokecolor := 1
            end

            if latex
                xguide := L"x_1"
                yguide := L"x_2"
            else
                xguide := "x₁"
                yguide := "x₂"
            end
            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end

        @series begin
            subplot := 2
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            ylims  := :auto
            yformatter := :scientific
            sol.t[0:nplot:end], ΔH[0:nplot:end]
        end
    end


    """
    Plots the solution of a massless charged particle.

    Arguments:
    * `sol <: Solution`
    * `params <: NamedTuple`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `xlims=:auto`: xlims for solution plot
    * `ylims=:auto`: ylims for solution plot
    * `latex=true`: use LaTeX guides
    """
    @userplot Plot_Massless_Charged_Particle_Solution
    @recipe function f(p::Plot_Massless_Charged_Particle_Solution; nplot=1, xlims=:auto, ylims=:auto, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        params = p.args[2]

        if sol.nt ≤ 200
            markersize := 5
        else
            markersize  := 1
            markercolor := 1
            linecolor   := 1
            markerstrokewidth := 1
            markerstrokecolor := 1
        end

        legend := :none
        size := (400,400)

        # solution
        @series begin
            # seriestype := :scatter
            if latex
                xguide := L"x_1"
                yguide := L"x_2"
            else
                xguide := "x₁"
                yguide := "x₂"
            end
            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            guidefont := font(18)
            tickfont := font(12)
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end
    end

    """
    Plots time traces of the solution of a massless charged particle and its energy error.

    Arguments:
    * `sol <: Solution`
    * `params <: NamedTuple`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `latex=true`: use LaTeX guides
    """
    @userplot Plot_Massless_Charged_Particle_Traces
    @recipe function f(p::Plot_Massless_Charged_Particle_Traces; nplot=1, latex=true)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        params = p.args[2]

        H, ΔH = compute_energy_error(sol.t, sol.q, params);

        size   := (800,600)
        legend := :none
        guidefont := font(18)
        tickfont  := font(12)
        right_margin := 10mm

        # traces
        layout := @layout [x₁Plot
                            x₂Plot
                            EPlot]

        if latex
            ylabels = (L"x_1", L"x_2")
        else
            ylabels = ("x₁", "x₂")
        end

        for i in 1:2
            @series begin
                subplot := i
                yguide := ylabels[i]
                xlims  := (sol.t[0], Inf)
                xaxis := false
                sol.t[0:nplot:end], sol.q[i,0:nplot:end]
            end
        end

        @series begin
            subplot := 3
            if latex
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            else
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            sol.t[0:nplot:end], ΔH[0:nplot:end]
        end
    end

end
