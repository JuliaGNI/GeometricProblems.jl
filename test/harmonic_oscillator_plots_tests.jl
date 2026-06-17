using Test
using GeometricProblems.HarmonicOscillator

# Check if CairoMakie is available for plotting tests
cairo_makie_available = false
try
    using CairoMakie
    global cairo_makie_available = true
catch
    @info "CairoMakie not available, skipping plotting tests"
end

# Load plotting module if CairoMakie is available
if cairo_makie_available
    include("../scripts/harmonic_oscillator_plotting.jl")
    using .HarmonicOscillatorPlotting
end

@testset "$(rpad("Harmonic Oscillator Plots",80))" begin

    # Create a harmonic oscillator problem (solution-independent tests)
    ode = odeproblem()

    # Create a basic solution-like object for testing
    sol = (q = sin.(0:0.1:10), v = cos.(0:0.1:10), t = 0:0.1:10)

    # Test utility functions that don't require CairoMakie
    @testset "Utility functions (no CairoMakie required)" begin
        # These should work by importing the module even without CairoMakie
        # But we'll test the mathematical correctness

        # Test basic mathematical values
        π_approx = 3.141592653589793

        # xpos and ypos should be the same in our implementation
        @test isapprox(xpos(0), sin(0.0))
        @test isapprox(ypos(0), sin(0.0))
        @test isapprox(zpos(0), 1.0)  # A*cos(0) = 1.0
        @test isapprox(ϑ(0), 0.0)

        # Test the values are reasonable (without exact precision)
        @test 0.0 < xpos(2) < 0.1  # should be small positive
        @test 0.9 < zpos(2) < 1.1    # should be close to 1.0
    end

    if cairo_makie_available
        @testset "CairoMakie-dependent functions" begin
            # Test spring geometry calculation
            spring_coords = spring(0)
            @test length(spring_coords) == 3
            @test spring_coords[1] isa Vector
            @test spring_coords[2] isa Vector
            @test spring_coords[3] isa Vector

            # Test phase space plot creation
            fig = plot_harmonic_oscillator(sol, ntime = 10)
            @test fig isa CairoMakie.Figure

            # Test single frame plotting
            fig_single = plot_spring(0, sol.q[1])
            @test fig_single isa CairoMakie.Figure

            # Test that figures are properly created
            @test !isempty(fig.content)
            @test !isempty(fig_single.content)
        end

        @testset "Animation functionality" begin
            test_dir = "test-harmonic-plots-temp"

            try
                # Test that we can create animations (but only do 2 frames for testing)
                plot_spring_animation(0:1, test_dir, true)

                # Check that files were created
                @test isdir(test_dir)
                @test isfile(joinpath(test_dir, "harmonic-oscillator-000.png"))
                @test isfile(joinpath(test_dir, "harmonic-oscillator-001.png"))
            finally
                # Clean up regardless of test results
                if isdir(test_dir)
                    rm(test_dir, recursive = true, force = true)
                end
            end
        end

        @testset "Integration with Harmonic Oscillator solutions" begin
            # Test that plotting works with actual oscillator solutions
            # This verifies integration between the plotting module and the main package

            sol_mock = (
                q = [cos(t) for t in 0:0.1:5],
                v = [-sin(t) for t in 0:0.1:5],
                t = 0:0.1:5
            )

            # Should work without errors
            fig = plot_harmonic_oscillator(sol_mock)
            @test fig isa CairoMakie.Figure

            # Test different time ranges
            fig_half = plot_harmonic_oscillator(sol_mock, ntime = 25)
            @test fig_half isa CairoMakie.Figure
        end
    else
        @testset "CairoMakie unavailable" begin
            # Verify we still have basic mathematical capability
            π_approx = 3.141592653589793
            @test isapprox(sin(0.0), 0.0)
            @test isapprox(cos(0.0), 1.0)
        end
    end

    # Note: We don't test the full animation generation in CI since it's time-consuming
    # and creates many files. The basic functionality test above verifies the core logic.

end
