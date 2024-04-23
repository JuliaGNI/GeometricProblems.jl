using GeometricProblems.LinearWave: h, s, u₀, ∂h, ∂s, ∂u₀
using ForwardDiff
using Test

const test_point₁ = 0.5
const test_point₂ = 1.5
const test_point₃ = 2.5

function test_h_derivative(test_point)
    h_autodiff = ForwardDiff.derivative(h, test_point)
    @test h_autodiff ≈ h(test_point)
end

function test_s_derivative(test_point)
    s_closure(ξ) = s(ξ, lw.μ)
    s_autodiff = ForwardDiff.derivative(s_closure, test_point)
    @test s_autodiff ≈ ∂s(test_point, lw.μ)
end

function test_u₀_derivative(test_point)
    u₀_closure(ξ) = u₀(ξ, lw.μ)
    u₀_autodiff = ForwardDiff.derivative(u₀_closure, test_point)
    @test u₀_autodiff ≈ ∂u₀(test_point, lw.μ)
end

function test_all_derivatives(test_point)
    test_h_derivative(test_point)
    test_s_derivative(test_point)
    test_u₀_derivative(test_point)
end

test_all_derivatives(test_point₁)
test_all_derivatives(test_point₂)
test_all_derivatives(test_point₃)