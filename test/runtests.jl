using Test

modname(fname)=splitext(basename(fname))[1]

function run_testdir(testdir,prefix)
    examples=modname.(readdir(testdir))
    for example in examples
        println("$(example):")
        if example[1:length(prefix)]==prefix
            path=joinpath(testdir,"$(example).jl")
            @eval begin
                include($path)
                print("  compile:")
                @time @test eval(Meta.parse("$($example).test()"))
                print("      run:")
                @time eval(Meta.parse("$($example).test()"))
            end
        end
    end
end


@time begin
    run_testdir(@__DIR__,"test_")
    run_testdir(joinpath(@__DIR__,"..","examples"),"Example")
end

