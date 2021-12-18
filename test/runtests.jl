using Test
# Activate assembly loop allocation checking
# as default.
ENV["VORONOIFVM_CHECK_ALLOCS"]="true"

modname(fname)=splitext(basename(fname))[1]

#
# Include all Julia files in `testdir` whose name starts with `prefix`,
# Each file `prefixModName.jl` must contain a module named
# `prefixModName` which has a method test() returning true
# or false depending on success.
#
function run_tests_from_directory(testdir,prefix)
    println("Directory $(testdir):")
    examples=modname.(readdir(testdir))
    for example in examples
        if length(example)>=length(prefix) &&example[1:length(prefix)]==prefix
            println("  $(example):")
            path=joinpath(testdir,"$(example).jl")
            @eval begin
                include($path)
                # Compile + run test
                @test eval(Meta.parse("$($example).test()"))
                # Second run: pure execution time.
                @time eval(Meta.parse("$($example).test()"))
            end
        end
    end
end


function run_all_tests()

    ENV["VORONOIFVM_CHECK_ALLOCS"]="true"
    notebooks=["nbproto.jl",
               "api-update.jl",
               "flux-reconstruction.jl",
               "problemcase.jl"
               ]
    
    @time begin
        @testset "basictest" begin
            run_tests_from_directory(@__DIR__,"test_")
        end
        begin
            for notebook in notebooks
                @testset "$(notebook)" begin
                    include(joinpath(@__DIR__,"..","pluto-examples",notebook))
                end
            end
        end
        
        @testset "1D Examples" begin
            run_tests_from_directory(joinpath(@__DIR__,"..","examples"),"Example1")
        end
        @testset "2D Examples" begin
            run_tests_from_directory(joinpath(@__DIR__,"..","examples"),"Example2")
        end
        @testset "3D Examples" begin
            run_tests_from_directory(joinpath(@__DIR__,"..","examples"),"Example3")
        end
        @testset "Misc Examples" begin
            run_tests_from_directory(joinpath(@__DIR__,"..","examples"),"Example4")
        end
    end
end

run_all_tests()
