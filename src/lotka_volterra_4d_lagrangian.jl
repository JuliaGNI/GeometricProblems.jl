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

using EulerLagrange
using GeometricSolutions
using LinearAlgebra
using Parameters

# export hamiltonian, ϑ, ϑ₁, ϑ₂, ω

export lotka_volterra_4d_ode,
    lotka_volterra_4d_lode,
    lotka_volterra_4d_ldae


# include("lotka_volterra_4d_plots.jl")


const Δt = 0.01
const nt = 1000
const timespan = (0.0, Δt * nt)

const q₀ = [2.0, 1.0, 1.0, 1.0]

const default_parameters = (a₁=1.0, a₂=1.0, a₃=1.0, a₄=1.0, b₁=-1.0, b₂=-2.0, b₃=-1.0, b₄=-1.0)
const reference_solution = [0.5988695239096916, 2.068567531039674, 0.2804351458645534, 1.258449091830993]

const A_antisym = 1 // 2 .* [
    0 -1 +1 -1
    +1 0 -1 +1
    -1 +1 0 -1
    +1 -1 +1 0]

const A_positive = [
    0 0 1 0
    1 0 0 1
    0 1 0 0
    1 0 1 0]

const A_upper = [
    0 -1 +1 -1
    0 0 -1 +1
    0 0 0 -1
    0 0 0 0]

const A_lower = [
    0 0 0 0
    +1 0 0 0
    -1 +1 0 0
    +1 -1 +1 0]

const A_quasicanonical_antisym = 1 // 2 .* [
    0 -1 +1 -1
    +1 0 -1 0
    -1 +1 0 -1
    +1 0 +1 0]

const A_quasicanonical_reduced = 1 // 2 .* [
    0 0 +1 0
    +2 0 -2 0
    -1 0 0 0
    +2 0 +2 0]

const B = [0 1 1 1
    1 0 1 1
    1 1 0 1
    1 1 1 0]

const A_default = A_antisym

const B_default = zero(B)


function get_parameters(p)
    (a=[p.a₁, p.a₂, p.a₃, p.a₄],
        b=[p.b₁, p.b₂, p.b₃, p.b₄])
end

H(x, a, b) = a ⋅ x + b ⋅ log.(x)
K(x, v, A, B) = log.(x) ⋅ (A * (v ./ x)) + x ⋅ (B * v)
L(x, v, A, B, a, b) = K(x, v, A, B) - H(x, a, b)


function lagrangian_system(A, B, parameters)
    t, x, v = lagrangian_variables(4)
    sparams = symbolize(parameters)

    Ks = K(x, v, A, B)
    Hs = H(x, get_parameters(sparams)...)

    DegenerateLagrangianSystem(Ks, Hs, t, x, v, sparams)
end

function initial_momentum(lag_sys, t₀, q₀, params)
    functions(lag_sys).p(t₀, q₀, zero(q₀), params)
end


function odeproblem(q₀=q₀, A=A_default, B=B_default; timespan=timespan, timestep=Δt, parameters=default_parameters)
    lag_sys = lagrangian_system(A, B, parameters)
    ODEProblem(lag_sys, timespan, timestep, q₀; parameters=parameters)
end

function lodeproblem(q₀=q₀, A=A_default, B=B_default; timespan=timespan, timestep=Δt, parameters=default_parameters)
    lag_sys = lagrangian_system(A, B, parameters)
    p₀ = initial_momentum(lag_sys, timespan[begin], q₀, parameters)
    LODEProblem(lag_sys, timespan, timestep, q₀, p₀; parameters=parameters)
end

function ldaeproblem(q₀=q₀, A=A_default, B=B_default; timespan=timespan, timestep=Δt, parameters=default_parameters)
    lag_sys = lagrangian_system(A, B, parameters)
    p₀ = initial_momentum(lag_sys, timespan[begin], q₀, parameters)
    LDAEProblem(lag_sys, timespan, timestep, q₀, p₀, zero(q₀); parameters=parameters)
end

end
