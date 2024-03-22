@doc raw"""
    OuterSolarSystem

The OuterSolarSystem module provides function "hodeproblem" and "lodeproblem"(not working at the moment).
    "
"""

using EulerLagrange
using LinearAlgebra
using Parameters
using GeometricProblems

const tspan = (0.0, 3)
const tstep = 0.5

const G = 2.95912208286e-4

#Sun
const m₁ =  1.0
const q₁ = [0., 0., 0.]
const q̇₁ = [0., 0., 0.]
const p₁ = [0., 0., 0.]

#Jupiter
const m₂ = 0.000954786104043 
const q₂ = [-3.5023653, -3.8169847, -1.5507963]
const q̇₂ = [0.00565429, -0.0041249, -0.00190589]
const p₂ = m₂*q̇₂
#Saturn 
const m₃ = 0.000285583733151
const q₃ = [9.0755314,-3.0458353,-1.6483708]
const q̇₃ = [0.00168318,0.00483525,0.00192462]
const p₃  = m₃ * q̇₃ 

#Uranus
const m₄ = 0.0000437273164546
const q₄ = [8.3101420,-16.2901086,-7.2521278]
const q̇₄ = [0.00354178,0.00137102,0.00055029]
const p₄ = m₄*q̇₄

#Neptune
const m₅ = 0.0000517759138449
const q₅ = [11.4707666,-25.7294829,-10.8169456]
const q̇₅ = [0.00288930,0.00114527,0.00039677]
const p₅ = m₅*q̇₅

#Pluto
const m₆ = 1/(1.3e8)
const q₆ = [-15.5387357,-25.2225594,-3.1902382]
const q̇₆ = [0.00276725,-0.00170702,-0.00136504]
const p₆ = m₆*q̇₆


const m = [m₁,m₂,m₃,m₄,m₅,m₆]
const q₀ = [q₁;q₂;q₃;q₄;q₅;q₆]
const q̇₀ = [q̇₁;q̇₂;q̇₃;q̇₄;q̇₅;q̇₆]
const p₀ = [p₁;p₂;p₃;p₄;p₅;p₆]

const default_parameters = (
    G = 2.95912208286e-4,
    m = [m₁,m₂,m₃,m₄,m₅,m₆]
)

function hamiltonian(t, q, p, params;d=3,n=6)
    @unpack G, m = params

    q = reshape(q,d,n)
    p = reshape(p,d,n)

    T = 0
    for i in 1:n
        for j in 1:i
            T = T + G*(m[i]*m[j])/norm(q[:,i]-q[:,j])
        end
    end

    1/2 * sum([p[:,i]'*p[:,i]/m[i] for i in 1:n]) - T 
end


function lagrangian(t, q, q̇, params;d=3,n=6)
    @unpack G, m = params
    
    q = reshape(q,d,n)
    q̇ = reshape(q̇,d,n)

    T = 0
    for i in 1:n
        for j in 1:i
            T = T + G*(m[i]*m[j])/norm(q[:,i]-q[:,j])
        end
    end

    1/2 * sum([m[i]*q̇[:,i]'*q̇[:,i] for i in 1:n]) + T 
end

function hodeproblem(q₀ = q₀, p₀ = p₀; tspan = tspan, tstep = tstep, params = default_parameters)
    t, q, p = hamiltonian_variables(18)
    sparams = symbolize(params)
    ham_sys = EulerLagrange.HamiltonianSystem(hamiltonian(t, q, p, sparams), t, q, p, sparams)
    HODEProblem(ham_sys, tspan, tstep, q₀, p₀; parameters = params)
end

function lodeproblem(q₀ = q₀, p₀ = p₀; tspan = tspan, tstep = tstep, params = default_parameters)
    t, x, v = lagrangian_variables(18)
    sparams = symbolize(params)
    lag_sys = EulerLagrange.LagrangianSystem(lagrangian(t, x, v, sparams), t, x, v, sparams)
    LODEProblem(lag_sys, tspan, tstep, q₀, p₀; v̄ = v̄, parameters = params)
end

hode = hodeproblem()
# lode = lodeproblem()
emsol=integrate(hode, ExplicitMidpoint())
Rksol=integrate(hode, RK4())
