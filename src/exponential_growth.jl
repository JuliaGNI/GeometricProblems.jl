@doc raw"""
# Exponential Growth

"""
module ExponentialGrowth

    using GeometricEquations

    export exponential_growth_ode

    function exponential_growth_ode_f(f, t, x)
        f[1] = x[1]
        nothing
    end

    function exponential_growth_ode(x₀=[1.0])
        ODE(exponential_growth_ode_f, x₀)
    end

end
