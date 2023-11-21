using Test
import ExampleJuggler
using ExampleJuggler: cleanexamples, @testmodules, @testscripts

ExampleJuggler.verbose!(true)
#
# Include all Julia files in `testdir` whose name starts with `prefix`,
# Each file `prefixModName.jl` must contain a module named
# `prefixModName` which has a method test() returning true
# or false depending on success.
#
function run_tests_from_directory(testdir, prefix)
    @info "Directory $(testdir):"
    examples = filter(ex -> length(ex) >= length(prefix) && ex[1:length(prefix)] == prefix, basename.(readdir(testdir)))
    @info examples
    @testmodules(testdir, examples)
end

function run_all_tests(; run_notebooks = false)
    @testset "basictest" begin
        run_tests_from_directory(@__DIR__, "test_")
    end

    @testset "Development Examples" begin
        run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example0")
    end
    @testset "MultiD Examples" begin
        run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example5")
    end
    @testset "1D Examples" begin
        run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example1")
    end
    @testset "2D Examples" begin
        run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example2")
    end
    @testset "3D Examples" begin
        run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example3")
    end

    @testset "Misc Examples" begin
        run_tests_from_directory(joinpath(@__DIR__, "..", "examples"), "Example4")
    end

    if run_notebooks
        notebooks = ["nbproto.jl",
            "outflow.jl",
            "flux-reconstruction.jl",
            "interfaces1d.jl",
            "problemcase.jl",
            "nonlinear-solvers.jl",
            "api-update.jl",
        ]
        @testset "Notebooks" begin
            @testscripts(joinpath(@__DIR__, "..", "pluto-examples"), notebooks)
        end
    end
end

run_all_tests(; run_notebooks = true)
