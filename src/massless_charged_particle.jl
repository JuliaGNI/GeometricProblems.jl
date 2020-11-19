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

    import ..Diagnostics: compute_invariant_error, compute_momentum_error

    export ϑ, A, B, ϕ, E, hamiltonian
    export massless_charged_particle_ode, massless_charged_particle_iode,
           massless_charged_particle_idae, massless_charged_particle_idae_spark
    export compute_energy_error, compute_momentum_error

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

    function ϑ(t::Number, q::AbstractVector, params::NamedTuple, k::Int)
        if k == 1
            ϑ₁(t, q, params)
        elseif k == 2
            ϑ₂(t, q, params)
        else
            nothing
        end
    end

    # derivatives of the one-form components
    dϑ₁dx₁(t, q, params) = - params[:A₀] * q[1] * q[2]
    dϑ₁dx₂(t, q, params) = - params[:A₀] * (1 + q[1]^2 + 3 * q[2]^2) / 2
    dϑ₂dx₁(t, q, params) = + params[:A₀] * (1 + 3 * q[1]^2 + q[2]^2) / 2
    dϑ₂dx₂(t, q, params) = + params[:A₀] * q[1] * q[2]

    # components of the force
    f₁(t, q, v, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₂dx₁(t, q, params) * v[2]
    f₂(t, q, v, params) = dϑ₁dx₂(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

    g₁(t, q, v, params) = dϑ₁dx₁(t, q, params) * v[1] + dϑ₁dx₂(t, q, params) * v[2]
    g₂(t, q, v, params) = dϑ₂dx₁(t, q, params) * v[1] + dϑ₂dx₂(t, q, params) * v[2]

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

    function massless_charged_particle_f̄(t, q, v, f, params)
        f[1] = - dHd₁(t, q, params)
        f[2] = - dHd₂(t, q, params)
        nothing
    end

    function massless_charged_particle_g(t, q, v, g, params)
        g[1] = f₁(t, q, v, params)
        g[2] = f₂(t, q, v, params)
        nothing
    end

    massless_charged_particle_g(t, q, p, v, g, params) = massless_charged_particle_g(t, q, v, g, params)

    function massless_charged_particle_u(t, q, v, u, params)
        u .= v
        nothing
    end

    massless_charged_particle_u(t, q, p, v, u, params) = massless_charged_particle_u(t, q, v, u, params)

    function massless_charged_particle_ϕ(t, q, p, ϕ, params)
        ϕ[1] = p[1] - ϑ₁(t,q,params)
        ϕ[2] = p[2] - ϑ₂(t,q,params)
        nothing
    end

    function massless_charged_particle_ψ(t, q, p, v, f, ψ, params)
        ψ[1] = f[1] - g₁(t,q,v,params)
        ψ[2] = f[2] - g₂(t,q,v,params)
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
                v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
                parameters=params, h=hamiltonian)
    end

    "Creates an implicit DAE object for the massless charged particle in 2D."
    function massless_charged_particle_idae(q₀=q₀; params=parameters)
        IDAE(massless_charged_particle_ϑ, massless_charged_particle_f,
                massless_charged_particle_u, massless_charged_particle_g,
                massless_charged_particle_ϕ, q₀, ϑ(0., q₀, params), zero(q₀);
                v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
                h=hamiltonian, parameters=params)
    end

    "Creates an implicit DAE object for the massless charged particle in 2D."
    function massless_charged_particle_idae_spark(q₀=q₀; params=parameters)
        IDAE(massless_charged_particle_ϑ, massless_charged_particle_f̄,
                massless_charged_particle_u, massless_charged_particle_g,
                massless_charged_particle_ϕ, q₀, ϑ(0., q₀, params), zero(q₀);
                v̄=massless_charged_particle_v, f̄=massless_charged_particle_f,
                h=hamiltonian, parameters=params)
    end



    compute_energy_error(t,q,params) = compute_invariant_error(t,q, (t,q) -> hamiltonian(t,q,params))
    compute_momentum_error(t,q,p,params::NamedTuple) = compute_momentum_error(t, q, p, (t,q,k) -> ϑ(t,q,params,k))


    """
    Plots the solution of a massless charged particle together with the energy error.

    Arguments:
    * `sol <: Solution`
    * `params <: NamedTuple`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `xlims=:auto`: xlims for solution plot
    * `ylims=:auto`: ylims for solution plot
    """
    @userplot Plot_Massless_Charged_Particle
    @recipe function f(p::Plot_Massless_Charged_Particle; xlims=:auto, ylims=:auto, elims=:auto, nplot=1)
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
                markersize  := .5
                markercolor := 1
                linecolor   := 1
                markerstrokewidth := .5
                markerstrokecolor := 1
            end
            if nplot > 1
                seriestype := :scatter
            end

            if backend() == Plots.PlotlyJSBackend()
                xguide := "x₁"
                yguide := "x₂"
            else
                xguide := L"x_1"
                yguide := L"x_2"
            end
            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end

        @series begin
            subplot := 2
            if backend() == Plots.PlotlyJSBackend()
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            else
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            end
            xlims  := (sol.t[0], Inf)
            ylims  := elims
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
    """
    @userplot Plot_Massless_Charged_Particle_Solution
    @recipe function f(p::Plot_Massless_Charged_Particle_Solution; nplot=1, xlims=:auto, ylims=:auto)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        params = p.args[2]

        if sol.nt ≤ 200
            markersize := 5
        else
            markersize  := .5
            markercolor := 1
            linecolor   := 1
            markerstrokewidth := .5
            markerstrokecolor := 1
        end
        if nplot > 1
            seriestype := :scatter
        end

        legend := :none
        size := (400,400)

        if backend() == Plots.PlotlyJSBackend() || backend() == Plots.GRBackend()
            guidefont := font(18)
            tickfont  := font(12)
        else
            guidefont := font(14)
            tickfont  := font(10)
        end

        # solution
        @series begin
            if backend() == Plots.PlotlyJSBackend()
                xguide := "x₁"
                yguide := "x₂"
            else
                xguide := L"x_1"
                yguide := L"x_2"
            end

            xlims  := xlims
            ylims  := ylims
            aspect_ratio := 1
            sol.q[1,0:nplot:end], sol.q[2,0:nplot:end]
        end
    end


    """
    Plots time traces of the energy error of a massless charged particle.

    Arguments:
    * `sol <: Solution`
    * `params <: NamedTuple`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    """
    @userplot Plot_Massless_Charged_Particle_Energy_Error
    @recipe function f(p::Plot_Massless_Charged_Particle_Energy_Error; nplot=1)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        params = p.args[2]

        H, ΔH = compute_energy_error(sol.t, sol.q, params);

        size   := (800,200)
        legend := :none
        right_margin := 10mm

        if backend() == Plots.GRBackend()
            left_margin := 10mm
            bottom_margin := 10mm
        end

        guidefont := font(12)
        tickfont  := font(10)

        @series begin
            if backend() == Plots.PlotlyJSBackend()
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            else
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            end

            xlims  := (sol.t[0], Inf)
            yformatter := :scientific
            sol.t[0:nplot:end], ΔH[0:nplot:end]
        end
    end


    """
    Plots time traces of the momentum error of a massless charged particle.

    Arguments:
    * `sol <: Solution`
    * `params <: NamedTuple`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    * `k=0`: index of momentum component (0 plots all components)
    """
    @userplot Plot_Massless_Charged_Particle_Momentum_Error
    @recipe function f(p::Plot_Massless_Charged_Particle_Momentum_Error; k=0, nplot=1)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        params = p.args[2]

        Δp = compute_momentum_error(sol.t, sol.q, sol.p, params)

        ntrace = (k == 0 ?    Δp.nd  :    1 )
        trange = (k == 0 ? (1:Δp.nd) : (k:k))

        size   := (800,200*ntrace)
        legend := :none
        layout := (ntrace,1)

        right_margin := 10mm

        if ntrace == 1 && backend() == Plots.GRBackend()
            left_margin := 10mm
        end

        guidefont := font(12)
        tickfont  := font(10)

        for i in trange
            @series begin
                if k == 0
                    subplot := i
                end
                if i == Δp.nd || k ≠ 0
                    if backend() == Plots.PlotlyJSBackend()
                        xguide := "t"
                    else
                        xguide := L"t"
                    end
                else
                    xaxis := false
                    bottom_margin := 10mm
                end
                # if i < Δp.nd && k == 0
                #     if backend() == Plots.PGFPlotsXBackend()
                #         bottom_margin := -6mm
                #     elseif backend() == Plots.GRBackend()
                #         bottom_margin := -3mm
                #     end
                # end
                # if i > 1 && k == 0
                #     if backend() == Plots.PGFPlotsXBackend()
                #         top_margin := -6mm
                #     elseif backend() == Plots.GRBackend()
                #         top_margin := -3mm
                #     end
                # end
                if backend() == Plots.PlotlyJSBackend()
                    yguide := "p" * subscript(i) * "(t) - ϑ" * subscript(i) * "(t)"
                else
                    yguide := latexstring("p_$i (t) - \\vartheta_$i (t)")
                end
                xlims  := (sol.t[0], Inf)
                yformatter := :scientific
                sol.t[0:nplot:end], Δp[i,0:nplot:end]
            end
        end
    end


    """
    Plots time traces of the solution of a massless charged particle trajectory and its energy error.

    Arguments:
    * `sol <: Solution`
    * `params <: NamedTuple`

    Keyword aguments:
    * `nplot=1`: plot every `nplot`th time step
    """
    @userplot Plot_Massless_Charged_Particle_Traces
    @recipe function f(p::Plot_Massless_Charged_Particle_Traces; xlims=:auto, ylims=:auto, elims=:auto, nplot=1)
        if length(p.args) != 2 || !(typeof(p.args[1]) <: Solution) || !(typeof(p.args[2]) <: NamedTuple)
            error("Massless charged particle plots should be given two arguments: a solution and a parameter tuple. Got: $(typeof(p.args))")
        end
        sol = p.args[1]
        params = p.args[2]

        H, ΔH = compute_energy_error(sol.t, sol.q, params);

        size   := (800,600)
        legend := :none
        right_margin := 10mm

        # traces
        layout := @layout [x₁Plot
                            x₂Plot
                            EPlot]

        if backend() == Plots.PlotlyJSBackend() || backend() == Plots.GRBackend()
            guidefont := font(14)
            tickfont  := font(10)
        else
            guidefont := font(12)
            tickfont  := font(10)
        end

        if backend() == Plots.PlotlyJSBackend()
            ylabels = ("x₁", "x₂")
        else
            ylabels = (L"x_1", L"x_2")
        end

        lims = (xlims, ylims, elims)

        for i in 1:2
            @series begin
                subplot := i
                yguide := ylabels[i]
                xlims  := (sol.t[0], Inf)
                ylims  := lims[i]
                xaxis  := false

                # if i == 1 && backend() == Plots.PGFPlotsXBackend()
                #     bottom_margin := -4mm
                # elseif i == 1 && backend() == Plots.GRBackend()
                #     bottom_margin := -2mm
                # end
                # if i == 2 && backend() == Plots.PGFPlotsXBackend()
                #     top_margin := -2mm
                #     bottom_margin := -2mm
                # elseif i == 2 && backend() == Plots.GRBackend()
                #     top_margin := -1mm
                #     bottom_margin := -1mm
                # end

                sol.t[0:nplot:end], sol.q[i,0:nplot:end]
            end
        end

        @series begin
            subplot := 3
            # if backend() == Plots.PGFPlotsXBackend()
            #     top_margin := -4mm
            # elseif backend() == Plots.GRBackend()
            #     top_margin := -2mm
            # end
            if backend() == Plots.PlotlyJSBackend()
                xguide := "t"
                yguide := "[H(t) - H(0)] / H(0)"
            else
                xguide := L"t"
                yguide := L"[H(t) - H(0)] / H(0)"
            end
            xlims := (sol.t[0], Inf)
            ylims := elims
            yformatter := :scientific
            sol.t[0:nplot:end], ΔH[0:nplot:end]
        end
    end

end
