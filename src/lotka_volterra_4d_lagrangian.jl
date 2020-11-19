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

    using ModelingToolkit
    using Parameters
    using RuntimeGeneratedFunctions
    using SymEngine

    using GeometricIntegrators.Equations
    using GeometricIntegrators.Solutions

    RuntimeGeneratedFunctions.init(@__MODULE__)

    # export hamiltonian, ϑ, ϑ₁, ϑ₂, ω

    export lotka_volterra_4d_ode, 
           lotka_volterra_4d_iode, lotka_volterra_4d_idae,
           lotka_volterra_4d_vode
        #    lotka_volterra_4d_vdae,
        #    lotka_volterra_4d_dg


    # include("lotka_volterra_4d_plots.jl")
       

    Δt = 0.1
    nt = 1000

    const q₀ = [ 2.0,  1.0,  1.0,  1.0]
    const a₀ = [ 1.0,  1.0,  1.0,  1.0]
    const b₀ = [-1.0, -2.0, -1.0, -1.0]

    const reference_solution = [1.6390462434739954, 1.3764800055785835, 0.37903204434372284, 1.4399236281802124]

    const A_antisym = 1//2 * [
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
    
    const A_quasicanonical_antisym = 1//2 * [
                0 -1 +1 -1
               +1  0 -1  0
               -1 +1  0 -1
               +1  0 +1  0]

    const A_quasicanonical_reduced = 1//2 * [
                0  0 +1  0
               +2  0 -2  0
               -1  0  0  0
               +2  0 +2  0]

    const B = [ 0  1  1  1
                1  0  1  1
                1  1  0  1
                1  1  1  0]


    # TODO: Remove H₁ and H₂ as soon as ModelingToolkit is fixed.
    H(x, a, b) = sum(a .* x) + sum(b .* log.(x))
    K(x, v, A, B) = sum(log.(x) .* (A * (v ./ x))) + sum(x .* (B * v))
    L(x, v, A, B, a, b) = K(x, v, A, B) - H(x, a, b)

    # function poisson(x, A, B)
    #     x₁, x₂, x₃, x₄ = symbols("x₁, x₂, x₃, x₄")
    #     Ω = ( [(A[1,1] - A[1,1])/x₁/x₁   (A[1,2] - A[2,1])/x₁/x₂   (A[1,3] - A[3,1])/x₁/x₃   (A[1,4] - A[4,1])/x₁/x₄
    #            (A[2,1] - A[1,2])/x₂/x₁   (A[2,2] - A[2,2])/x₂/x₂   (A[2,3] - A[3,2])/x₂/x₃   (A[2,4] - A[4,2])/x₂/x₄
    #            (A[3,1] - A[1,3])/x₃/x₁   (A[3,2] - A[2,3])/x₃/x₂   (A[3,3] - A[3,3])/x₃/x₃   (A[3,4] - A[4,3])/x₃/x₄
    #            (A[4,1] - A[1,4])/x₄/x₁   (A[4,2] - A[2,4])/x₄/x₂   (A[4,3] - A[3,4])/x₄/x₃   (A[4,4] - A[4,4])/x₄/x₄]
    #        .+ [(B[1,1] - B[1,1])         (B[1,2] - B[2,1])         (B[1,3] - B[3,1])         (B[1,4] - B[4,1])
    #            (B[2,1] - B[1,2])         (B[2,2] - B[2,2])         (B[2,3] - B[3,2])         (B[2,4] - B[4,2])
    #            (B[3,1] - B[1,3])         (B[3,2] - B[2,3])         (B[3,3] - B[3,3])         (B[3,4] - B[4,3])
    #            (B[4,1] - B[1,4])         (B[4,2] - B[2,4])         (B[4,3] - B[3,4])         (B[4,4] - B[4,4])]
    #         )

    #     P = expand.(inv(Ω))

    #     body = convert(Expr, P)
    #     code = quote
    #         let (x₁, x₂, x₃, x₄) = $x
    #             $body
    #         end
    #     end
        
    #     eval(code)
    # end

    function poisson(Ω, X)
        # TODO: Replace SymEngine with ModelingToolkit as soon as it is fixed.
        x₁, x₂, x₃, x₄ = symbols("x₁, x₂, x₃, x₄")
        code = ModelingToolkit.build_function(Ω, X...)[1]
        Ωfnc = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code))
        Ωsym = Ωfnc(x₁, x₂, x₃, x₄)
        Psym = expand.(inv(Ωsym))
        body = convert(Expr, Psym)
        code = quote
            let (x₁, x₂, x₃, x₄) = $X
                $body
            end
        end
        eval(code)
    end


    function substitute_ẋ_with_v!(eqs, ẋ, v)
        for i in eachindex(eqs)
            eqs[i] = substitute(eqs[i], [ẋᵢ=>vᵢ for (ẋᵢ,vᵢ) in zip(ẋ,v)])
        end
    end

    function substitute_inverses!(eqs, x, v, Dv)
        for i in eachindex(eqs)
            for j in eachindex(x,v,Dv)
                eqs[i] = substitute(eqs[i], expand_derivatives(Dv[j](v[j] / x[j]  )) => 1/x[j]  )
                eqs[i] = substitute(eqs[i], expand_derivatives(Dv[j](v[j] / x[j]^2)) => 1/x[j]^2)
                eqs[i] = substitute(eqs[i], inv(x[j]) => 1/x[j])
            end
        end
    end

    function substitute_variables!(eqs, x, v, X, V)
        for i in eachindex(eqs)
            eqs[i] = substitute(eqs[i], [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
            eqs[i] = substitute(eqs[i],  0//1 => Num( 0))
            eqs[i] = substitute(eqs[i],  1//1 => Num( 1))
            eqs[i] = substitute(eqs[i], -1//1 => Num(-1))
        end
    end


    function get_functions(A,B,a,b)
        @parameters t
        @derivatives Dt'~t

        @variables x₁(t), x₂(t), x₃(t), x₄(t)
        @variables v₁(t), v₂(t), v₃(t), v₄(t)
        @variables X₁, X₂, X₃, X₄
        @variables V₁, V₂, V₃, V₄
        @variables P₁, P₂, P₃, P₄
        @variables F₁, F₂, F₃, F₄

        x = [x₁, x₂, x₃, x₄]
        v = [v₁, v₂, v₃, v₄]
        X = [X₁, X₂, X₃, X₄]
        V = [V₁, V₂, V₃, V₄]
        P = [P₁, P₂, P₃, P₄]
        F = [F₁, F₂, F₃, F₄]

        Dx = Differential.(x)
        Dv = Differential.(v)
        DX = Differential.(X)
        DV = Differential.(V)

        let L = L(x, v, A, B, a, b), K = K(x, v, A, B), H = H(x, a, b)
            EL = [expand_derivatives(Dx[i](L) - Dt(Dv[i](L))) for i in eachindex(Dx,Dv)]
            ∇H = [expand_derivatives(Dx[i](H)) for i in eachindex(Dx)]
            f  = [expand_derivatives(Dx[i](L)) for i in eachindex(Dx)]
            f̄  = [expand_derivatives(Dx[i](K)) for i in eachindex(Dx)]
            g  = [expand_derivatives(Dt(Dv[i](L))) for i in eachindex(Dv)]
            ϑ  = [expand_derivatives(Dv[i](L)) for i in eachindex(Dv)]
            ω  = [expand_derivatives(simplify(Dx[i](ϑ[j]) - Dx[j](ϑ[i]))) for i in eachindex(Dx,ϑ), j in eachindex(Dx,ϑ)]

            for eq in (EL, ∇H, f, f̄, g, ϑ, ω, P)
                substitute_ẋ_with_v!(eq, Dt.(x), v)
                substitute_inverses!(eq, x, v, Dv) # TODO: Remove as soon as ModelingToolkit is fixed.
                substitute_variables!(eq, x, v, X, V)
            end

            Σ = poisson(ω, X)
            substitute_variables!(Σ, x, v, X, V)
            ẋ = simplify.(Σ * ∇H)

            H = substitute(H, [z=>Z for (z,Z) in zip([x..., v...], [X..., V...])])
            ϕ = P .- ϑ
            ψ = F .- g

            substitute_variables!(ẋ, x, v, X, V)
            substitute_inverses!(ẋ, X, V, DV)
            
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
            
            return (
                EL = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_EL)),
                ∇H = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_∇H)),
                f  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_f)),
                f̄  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_f̄)),
                g  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_g)),
                p  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_p)),
                ϑ  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_ϑ)),
                ω  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_ω)),
                P  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_P)),
                ẋ  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_ẋ)),
                ϕ  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_ϕ)),
                ψ  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_ψ)),
                H  = @RuntimeGeneratedFunction(ModelingToolkit.inject_registered_module_functions(code_H)),
            )
        end
    end


    function lotka_volterra_4d_ode(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀)
        funcs = get_functions(A,B,a,b)
        ODE((t,x,ẋ) -> funcs[:ẋ](ẋ,t,x), q₀; h=funcs[:H])
    end

    function lotka_volterra_4d_iode(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀)
        funcs = get_functions(A,B,a,b)
        IODE((t,x,v,ϑ) -> funcs[:ϑ](ϑ,t,x),
             (t,x,v,f) -> funcs[:f](f,t,x,v),
             (t,x,v,f̄) -> funcs[:f̄](f̄,t,x,v),
             q₀, funcs[:p](0, q₀);
             v̄=(t,x,v) -> funcs[:ẋ](v,t,x), h=funcs[:H])
    end

    function lotka_volterra_4d_vode(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀)
        funcs = get_functions(A,B,a,b)
        VODE((t,x,v,ϑ) -> funcs[:ϑ](ϑ,t,x),
             (t,x,v,f) -> funcs[:f](f,t,x,v),
             (t,x,v,f̄) -> funcs[:f̄](f̄,t,x,v),
             q₀, funcs[:p](0, q₀);
             v̄  = (t,x,v) -> funcs[:ẋ](v,t,x),
             h  = funcs[:H],
             Ω  = (t,x,ω) -> funcs[:ω](ω,t,x),
             ∇H = (t,x,∇H) -> funcs[:∇H](∇H,t,x))
    end

    function lotka_volterra_4d_idae(q₀=q₀, A=A_antisym, B=zeros(Int, 4, 4), a=a₀, b=b₀)
        funcs = get_functions(A,B,a,b)
        IDAE((t,x,v,ϑ) -> funcs[:ϑ](ϑ,t,x,v),
             (t,x,v,f) -> funcs[:f](f,t,x,v),
             (t,x,p,v,u) -> u .= v,
             (t,x,p,v,f̄) -> funcs[:f̄](f̄,t,x,v),
             (t,x,p,ϕ) -> funcs[:ϕ](ϕ,t,x,p),
             q₀, funcs[:p](0, q₀), zero(q₀);
             v̄=(t,x,v) -> funcs[:ẋ](v,t,x), h=funcs[:H])
    end

    # function lotka_volterra_4d_vdae(q₀=q₀, p₀=ϑ(0, q₀), λ₀=zero(q₀), params=p)
    #     VDAE(lotka_volterra_4d_ϑ, lotka_volterra_4d_f_ham,
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
