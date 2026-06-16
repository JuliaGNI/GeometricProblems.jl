using Test
using GeometricProblems.HarmonicOscillator
using GeometricProblems.HarmonicOscillatorPlots

@testset "$(rpad("Harmonic Oscillator Plots",80))" begin

    # Create a harmonic oscillator problem (solution-independent tests)
    ode = odeproblem()

    # Create a basic solution-like object for testing
    sol = (q = sin.(0:0.1:10), v = cos.(0:0.1:10), t = 0:0.1:10)

    # Test that the basic plotting functions work without errors
    @testset "Basic functionality" begin
        # Test utility functions (available without CairoMakie)
        @test xpos(0) ≈ 0
        @test ypos(0) ≈ 0
        @test zpos(0) ≈ 1.0  # A*cos(0) = 1.0
        @test ϑ(0) ≈ 0       # -ω*A*sin(0) = 0

        # Test spring geometry calculation
        spring_coords = spring(0)
        @test length(spring_coords) == 3

        # Only test CairoMakie-dependent functions if available
        if @isdefined(CairoMakie) && isdefined(GeometricProblems.HarmonicOscillatorPlots, :plot_harmonic_oscillator)
            # Test phase space plot
            fig = plot_harmonic_oscillator(sol, ntime = 10)
            @test fig isa CairoMakie.Figure

            # Test spring animation functions (single frames)
            fig_single = plot_spring(0, sol.q[1])
            @test fig_single isa CairoMakie.Figure
        else
            @test true # Skip these tests if CairoMakie not available
        end

        # Test that spring generation works
        spring_coords = spring(0)
        @test length(spring_coords) == 3
        @test spring_coords[1] isa Vector
        @test spring_coords[2] isa Vector
        @test spring_coords[3] isa Vector
    end

    @testset "Function utilities" begin
        # Test coordinate functions
        @test xpos(0) ≈ 0
        @test ypos(0) ≈ 0
        @test zpos(0) ≈ 1.0  # A*cos(0) = 1.0
        @test ϑ(0) ≈ 0       # -ω*A*sin(0) = 0
    end

    @testset "Animation target directory creation" begin
        # Only test animation functions if CairoMakie is available
        if @isdefined(CairoMakie) && isdefined(GeometricProblems.HarmonicOscillatorPlots, :plot_spring_animation)
            test_dir = "test-harmonic-plots"

            # Test that we can create animations (but only do 2 frames for testing)
            plot_spring_animation(0:1, test_dir, true)

            # Check that files were created
            @test isdir(test_dir)
            @test isfile(joinpath(test_dir, "harmonic-oscillator-000.png"))
            @test isfile(joinpath(test_dir, "harmonic-oscillator-001.png"))

            # Clean up
            rm(test_dir, recursive = true, force = true)
        else
            @test true # Skip if CairoMakie not available
        end
    end

    # Note: We don't test the full animation generation in CI since it's time-consuming
    # and creates many files. The basic functionality test above verifies the core logic.

end
