@doc raw"""
# Harmonic Oscillator

"""
module HarmonicOscillator

    using Reexport
    using GeometricIntegrators.Solutions

   @reexport using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem

    using GeometricIntegrators.TestProblems.HarmonicOscillatorProblem: hamiltonian

    export hamiltonian
    export compute_energy_error


    function compute_energy_error(t, q::DataSeries{T}) where {T}
        h = SDataSeries(T, q.nt)
        e = SDataSeries(T, q.nt)

        for i in axes(q,2)
            h[i] = hamiltonian(t[i], q[:,i])
            e[i] = (h[i] - h[0]) / h[0]
        end

        (h, e)
    end

end
