module LotkaVolterra2d

    using Reexport
    using GeometricIntegrators.Solutions

    @reexport using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem

    using GeometricIntegrators.TestProblems.LotkaVolterra2dProblem: hamiltonian, ϑ, ϑ₁, ϑ₂, ω

    export hamiltonian, ϑ, ϑ₁, ϑ₂, ω
    export compute_energy_error, compute_momentum_error


    function compute_energy_error(t, q::DataSeries{T}) where {T}
        h = SDataSeries(T, q.nt)
        e = SDataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i])
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

    function compute_momentum_error(t, q, p)
        p1_error = zeros(q.nt+1)
        p2_error = zeros(q.nt+1)

        for i in 1:(q.nt+1)
            p1_error[i] = p.d[1,i] - ϑ₁(t.t[i], q.d[:,i])
            p2_error[i] = p.d[2,i] - ϑ₂(t.t[i], q.d[:,i])
        end

        (p1_error, p2_error)
    end

end
