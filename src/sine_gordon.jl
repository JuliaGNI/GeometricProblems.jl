@doc raw"""
Applying multi-sypmplectic integrators to the nonlinear wave equation, sine-Gordon equation.
u_tt - u_xx = sin(u) = f'(u)
"""

module NonlinearWave

    const l=4*pi

    function f(u)
        1 .- cos.(u)
    end

    function df(u)
        sin.(u)
    end

    function sigma(x,u,ux,ut)
        ux .^ 2 
    end

    function lagrangian(x,u,ux,ut)
        1/2*ut.^2 - 1/2*ux.^2 + 1 - cos.(u)
    end

    function hamiltonian(x,u,ux,ut)
        1/2*ut.^2 + 1/2*ux.^2 - 1 + cos.(u)
    end

    function initial_position(x,u₀,v₀)
        2*exp(-(x-l/2).^2)
    end

    function initial_velocity(x,u₀,v₀)
        2*exp(-(x-l/2).^2) .* 2 * (-(x-l/2))
    end
end

