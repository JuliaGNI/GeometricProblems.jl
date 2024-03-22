# Outer Solar System

<!-- 
t_end = 100000 
t_step = 0.5 
hode = hodeproblem(tspan = (0,t_end), tstep = t_step)

sol = integrate(hode, RK4())
ImMi_sol = integrate(hode, ImplicitMidpoint())

using Plots
ImMi_sol_mat = hcat(ImMi_sol.q[0:200000]...)

plot(ImMi_sol_mat[1,:],ImMi_sol_mat[2,:],ImMi_sol_mat[3,:],label = "Sun")
plot!(ImMi_sol_mat[4,:],ImMi_sol_mat[5,:],ImMi_sol_mat[6,:],label ="Jupiter")
plot!(ImMi_sol_mat[7,:],ImMi_sol_mat[8,:],ImMi_sol_mat[9,:],label ="Saturn")
plot!(ImMi_sol_mat[10,:],ImMi_sol_mat[11,:],ImMi_sol_mat[12,:],label ="Uranus")
plot!(ImMi_sol_mat[13,:],ImMi_sol_mat[14,:],ImMi_sol_mat[15,:],label ="Neptune")
plot!(ImMi_sol_mat[16,:],ImMi_sol_mat[17,:],ImMi_sol_mat[18,:],label ="Pluto")
title!("Implicit Midpoint") 
-->

This Outer Solar System problem is simulating the trajactories of 6 planets:Sun,Jupiter,Saturn,Uranus,Neptune and Pluto. 

 ![](images/Solar_system.svg) 

Assuming the Sun is of mass $m_0=1$ while the rest $m_i$ are masses relative to the Sun. 

The Hamiltonian of the system reads 

```math
\begin{equation*}
H(p,q) = \frac{1}{2} \sum_{i=0}^{5} \frac{1}{m_i} p_i^T p_i - G  \sum_{i=1}^{5} \sum_{j=0}^{i-1} \frac{m_i m_j}{\lVert q_i - q_j \rVert},
\end{equation*}
```
where $q_i$ are distances in astronomical units (1[A.U.] =
149597870[km]), and the gravitational constant $G =  2.95912208286e^{−4}$. 

The generalized velocities are derived by 
```math
\begin{equation*}
\dot q_i = \frac{\partial H}{\partial p_i},
\end{equation*}
```
thus the Lagrangian:
```math
\begin{equation*}
L(q,\dot q) = \frac{1}{2} \sum_{i=0}^{5} m_i \dot q_i^T \dot q_i - G  \sum_{i=1}^{5} \sum_{j=0}^{i-1} \frac{m_i m_j}{\lVert q_i - q_j \rVert}\end{equation*}
```

## Library functions
```@docs
GeometricProblems.outer_solar_system
```

```@autodocs
Modules = [GeometricProblems.OuterSolarSystem]
```
