@doc raw"""
Applying multi-sypmplectic integrators to the nonlinear wave equation, sine-Gordon equation.
u_tt = u_xx - sin(u) 
"""

function f(u,ux,ut)
    1 .- cos.(u)
end

function sigma(u,ux,ut)
    ux .^ 2 
end




