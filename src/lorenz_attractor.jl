
module LorenzAttractor

    using GeometricIntegrators.Equations
    using Plots
    using IJulia

    export lorenz_attractor_ode, plot_lorenz_attractor

    x, y, z = 1., 1., 1.
    σ, ρ, β = 10., 28., 8/3

    const q₀ = [x, y, z]

    function lorenz_attractor_v(t, x, v)
        v[1] = σ * (x[2] - x[1])
        v[2] = x[1] * (ρ - x[3]) - x[2]
        v[3] = x[1] * x[2] - β * x[3]
        nothing
    end


    function lorenz_attractor_ode(q₀=q₀)
        ODE(lorenz_attractor_v, q₀)
    end


    function plot_lorenz_attractor(x)
        # initialize a 3D plot with 1 empty series
        default(size=(900, 750))

        p3d = path3d(1, xlim=(-25,25), ylim=(-25,25), zlim=(0,50),
                        xlab="x", ylab="y", zlab="z",
                        title="Lorenz Attractor", marker=1)

        push!(p3d, x[1,0], x[2,0], x[3,0])

        # initialize plots for planes
        pxz = plot([x[1,0]], [x[3,0]], title="x-z projection")
        pxy = plot([x[1,0]], [x[2,0]], title="x-y projection")
        pyz = plot([x[2,0]], [x[3,0]], title="y-z projection")

        # initialize master plot
        plt = plot(p3d, pxy, pxz, pyz, layout=(2,2), legend=false, fmt=:png)

        # build an animation
        for i in axes(x,2)
            push!(p3d, x[1,i], x[2,i], x[3,i])
            push!(pxz, x[1,i], x[3,i])
            push!(pxy, x[1,i], x[2,i])
            push!(pyz, x[2,i], x[3,i])

            if i % 500 == 0
                IJulia.clear_output(true)
                plt |> display
            end
        end
    end

end
