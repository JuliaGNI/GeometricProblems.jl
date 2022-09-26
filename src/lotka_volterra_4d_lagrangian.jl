@doc raw"""

Implements equations for Lagrangian Lotka-Volterra models in 4d of the form
```math
L(q, \dot{q}) = \frac{1}{2} (\log q)^T A \frac{\dot{q}}{q} + q^T B \dot{q} - H (q) ,
```
with hamiltonian
```math
H(q) = a_1 q^1 + a_2 q^2 + a_3 q^3 + a_4 q^4 + b_1 \log q^1 + b_2 \log q^2 + b_3 \log q^3 + b_4 \log q^4 + b_5 \log q^5 .
```

"""
module LotkaVolterra4dLagrangian

    using GeometricEquations
    using GeometricSolutions
    using ModelingToolkit
    using Parameters
    using RuntimeGeneratedFunctions
    using Symbolics

    RuntimeGeneratedFunctions.init(@__MODULE__)

    # export hamiltonian, ϑ, ϑ₁, ϑ₂, ω

    export lotka_volterra_4d_ode, 
           lotka_volterra_4d_iode, lotka_volterra_4d_idae,
           lotka_volterra_4d_lode
        #    lotka_volterra_4d_ldae,
        #    lotka_volterra_4d_dg


    # include("lotka_volterra_4d_plots.jl")
       

    const Δt = 0.01
    const nt = 1000
    const tspan = (0.0, Δt*nt)

    const q₀ = [ 2.0,  1.0,  1.0,  1.0]
    const a₀ = [ 1.0,  1.0,  1.0,  1.0]
    const b₀ = [-1.0, -2.0, -1.0, -1.0]

    const default_parameters = (a₁=1.0, a₂=1.0, a₃=1.0, a₄=1.0, b₁=-1.0, b₂=-2.0, b₃=-1.0, b₄=-1.0)
    const reference_solution = [0.5988695239096916, 2.068567531039674, 0.2804351458645534, 1.258449091830993]

    const A_antisym = 1//2 .* [
                0 -1 +1 -1
               +1  0 -1 +1
               -1 +1  0 -1
               +1 -1 +1  0]

    const A_positive = [
                0  0  1  0
                1  0  0  1
                0  1  0  0
                1  0  1  0]

    const A_upper = [
                0 -1 +1 -1
                0  0 -1 +1
                0  0  0 -1
                0  0  0  0]

    const A_lower = [
                0  0  0  0
               +1  0  0  0
               -1 +1  0  0
               +1 -1 +1  0]
    
    const A_quasicanonical_antisym = 1//2 .* [
                0 -1 +1 -1
               +1  0 -1  0
               -1 +1  0 -1
               +1  0 +1  0]

    const A_quasicanonical_reduced = 1//2 .* [
                0  0 +1  0
               +2  0 -2  0
               -1  0  0  0
               +2  0 +2  0]

    const B = [ 0  1  1  1
                1  0  1  1
                1  1  0  1
                1  1  1  0]


    H(x, a, b) = sum(a .* x) + sum(b .* log.(x))
    K(x, v, A, B) = sum(log.(x) .* (A * (v ./ x))) + sum(x .* (B * v))
    L(x, v, A, B, a, b) = K(x, v, A, B) - H(x, a, b)


    function substitute_ẋ_with_v!(eqs, ẋ, v)
        for i in eachindex(eqs)
            eqs[i] = substitute(eqs[i], [ẋᵢ=>vᵢ for (ẋᵢ,vᵢ) in zip(ẋ,v)])
        end
    end

    function substitute_variables!(eqs, x, v, X, V)
        for i in eachindex(eqs)
            eqs[i] = substitute(eqs[i], [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
        end
    end


    function get_functions(A,B,a,b)
        @parameters t

        @variables x₁(t), x₂(t), x₃(t), x₄(t)
        @variables v₁(t), v₂(t), v₃(t), v₄(t)
        @variables X[1:4]
        @variables V[1:4]
        @variables P[1:4]
        @variables F[1:4]

        x = [x₁, x₂, x₃, x₄]
        v = [v₁, v₂, v₃, v₄]

        Dt = Differential(t)
        Dx = Differential.(x)
        Dv = Differential.(v)
        DX = Differential.(X)
        DV = Differential.(V)

        let L = simplify(L(x, v, A, B, a, b)), K = simplify(K(x, v, A, B)), H = simplify(H(x, a, b))
            EL = [expand_derivatives(Dx[i](L) - Dt(Dv[i](L))) for i in eachindex(Dx,Dv)]
            ∇H = [expand_derivatives(dx(H)) for dx in Dx]
            f  = [expand_derivatives(dx(L)) for dx in Dx]
            f̄  = [expand_derivatives(dx(K)) for dx in Dx]
            g  = [expand_derivatives(Dt(dv(L))) for dv in Dv]
            ϑ  = [expand_derivatives(dv(L)) for dv in Dv]
            ω  = [expand_derivatives(simplify(Dx[i](ϑ[j]) - Dx[j](ϑ[i]))) for i in eachindex(Dx,ϑ), j in eachindex(Dx,ϑ)]
            Σ  = simplify.(inv(ω))
            ẋ  = simplify.(Σ * ∇H)

            for eq in (EL, ∇H, f, f̄, g, ϑ, ω, Σ, ẋ)
                substitute_ẋ_with_v!(eq, Dt.(x), v)
                substitute_variables!(eq, x, v, X, V)
            end

            H = substitute(H, [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
            L = substitute(L, [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
            ϕ = [P[i] - ϑ[i] for i in eachindex(P,ϑ)]
            ψ = [F[i] - g[i] for i in eachindex(F,g)]

            # component wise code generation
            # code_EL = [build_function(eq, t, X, V)[2] for eq in EL]
            # code_f  = [build_function(eq, t, X, V)[2] for eq in f ]
            # code_f̄  = [build_function(eq, t, X, V)[2] for eq in f̄ ]
            # code_g  = [build_function(eq, t, X, V)[2] for eq in g ]
            # code_∇H = [build_function(eq, t, X)[2]    for eq in ∇H]
            # code_p  = [build_function(eq, t, X)[1]    for eq in ϑ ]
            # code_ϑ  = [build_function(eq, t, X)[2]    for eq in ϑ ]
            # code_ω  = [build_function(eq, t, X)[2]    for eq in ω ]
            # code_P  = [build_function(eq, t, X)[2]    for eq in Σ ]
            # code_ẋ  = [build_function(eq, t, X)[2]    for eq in ẋ ]
            # code_ϕ  = [build_function(eq, t, X, P)[2]        for eq in ϕ ]
            # code_ψ  = [build_function(eq, t, X, V, P, F)[2]  for eq in ψ ]
            # code_H  =  build_function(H,  t, X)

            code_EL = build_function(EL, t, X, V)[2]
            code_f  = build_function(f,  t, X, V)[2]
            code_f̄  = build_function(f̄,  t, X, V)[2]
            code_g  = build_function(g,  t, X, V)[2]
            code_∇H = build_function(∇H, t, X)[2]
            code_p  = build_function(ϑ,  t, X)[1]
            code_ϑ  = build_function(ϑ,  t, X)[2]
            code_ω  = build_function(ω,  t, X)[2]
            code_P  = build_function(Σ,  t, X)[2]
            code_ẋ  = build_function(ẋ,  t, X)[2]
            code_ϕ  = build_function(ϕ,  t, X, P)[2]
            code_ψ  = build_function(ψ,  t, X, V, P, F)[2]
            code_H  = build_function(H,  t, X)
            code_L  = build_function(L,  t, X, V)
            
            return (
                EL = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_EL)),
                ∇H = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_∇H)),
                f  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_f)),
                f̄  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_f̄)),
                g  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_g)),
                p  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_p)),
                ϑ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ϑ)),
                ω  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ω)),
                P  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_P)),
                ẋ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ẋ)),
                ϕ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ϕ)),
                ψ  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_ψ)),
                H  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_H)),
                L  = @RuntimeGeneratedFunction(Symbolics.inject_registered_module_functions(code_L)),
            )
        end
    end


    function lotka_volterra_4d_ode(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀; tspan=tspan, tstep=Δt)
        funcs = get_functions(A,B,a,b)
        GeometricEquations.ODEProblem((t,x,ẋ,params) -> funcs[:ẋ](ẋ,t,x), tspan, tstep, q₀;
                    invariants = (h = funcs[:H],))
    end

    function lotka_volterra_4d_iode(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀; tspan=tspan, tstep=Δt)
        funcs = get_functions(A,B,a,b)
        IODEProblem((t,x,v,ϑ,params) -> funcs[:ϑ](ϑ,t,x),
                    (t,x,v,f,params) -> funcs[:f](f,t,x,v),
                    (t,x,v,f̄,params) -> funcs[:f̄](f̄,t,x,v),
                    tspan, tstep, q₀, funcs[:p](0, q₀);
                    v̄ = (t,x,v,params) -> funcs[:ẋ](v,t,x),
                    invariants = (h = (t,x,v,params) -> funcs[:H](t,x),))
    end

    function lotka_volterra_4d_lode(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀; tspan=tspan, tstep=Δt)
        funcs = get_functions(A,B,a,b)
        LODEProblem((t,x,v,ϑ,params) -> funcs[:ϑ](ϑ,t,x),
                    (t,x,v,f,params) -> funcs[:f](f,t,x,v),
                    (t,x,v,f̄,params) -> funcs[:f̄](f̄,t,x,v),
                    (t,x,v,params)   -> funcs[:L](t,x,v),
                    (t,x,v,ω,params) -> funcs[:ω](ω,t,x),
                    tspan, tstep, q₀, funcs[:p](0, q₀);
                    v̄ = (t,x,v,params) -> funcs[:ẋ](v,t,x),
                    invariants = (h = (t,x,v,params) -> funcs[:H](t,x),))
    end

    function lotka_volterra_4d_idae(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀; tspan=tspan, tstep=Δt)
        funcs = get_functions(A,B,a,b)
        IDAEProblem((t,x,v,ϑ,params) -> funcs[:ϑ](ϑ,t,x,v),
                    (t,x,v,f,params) -> funcs[:f](f,t,x,v),
                    (t,x,p,v,u,params) -> u .= v,
                    (t,x,p,v,f̄,params) -> funcs[:f̄](f̄,t,x,v),
                    (t,x,p,ϕ,params) -> funcs[:ϕ](ϕ,t,x,p),
                    tspan, tstep, q₀, funcs[:p](0, q₀), zero(q₀);
                    v̄ = (t,x,v,params) -> funcs[:ẋ](v,t,x),
                    invariants = (h = (t,x,v,params) -> funcs[:H](t,x),))
    end

    # function lotka_volterra_4d_ldae(q₀=q₀, p₀=ϑ(0, q₀), λ₀=zero(q₀), params=p)
    #     LDAE(lotka_volterra_4d_ϑ, lotka_volterra_4d_f_ham,
    #          lotka_volterra_4d_g, lotka_volterra_4d_g̅,
    #          lotka_volterra_4d_ϕ, lotka_volterra_4d_ψ,
    #          q₀, p₀, λ₀; parameters=params, h=hamiltonian, v=lotka_volterra_4d_v)
    # end

    # function lotka_volterra_4d_dg(q₀=q₀, p₀=ϑ(0, q₀), params=p)
    #     IODE(lotka_volterra_4d_ϑ, lotka_volterra_4d_f,
    #          lotka_volterra_4d_g, q₀, p₀;
    #          parameters=params, h=hamiltonian, v=lotka_volterra_4d_v)
    # end

end
