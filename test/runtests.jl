using Test
using Pluto
using Pkg
using VoronoiFVM


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
        if length(example)>=length(prefix) && example[1:length(prefix)]==prefix
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

    # de-markdown eventual cells with Pkg.develop and write
    # to pluto-tmp.jl
    # notebook=Pluto.load_notebook_nobackup(input)
    # pkgcellfortest=findfirst(c->occursin("Pkg.develop",c.code),notebook.cells)
    # if  pkgcellfortest!=nothing
    #     # de-markdown pkg cell
    #     notebook.cells[pkgcellfortest].code=replace(notebook.cells[pkgcellfortest].code,"md"=>"")
    #     notebook.cells[pkgcellfortest].code=replace(notebook.cells[pkgcellfortest].code,"\"\"\""=>"")
    #     notebook.cells[pkgcellfortest].code=replace(notebook.cells[pkgcellfortest].code,";"=>"")
    #     @info "Pkg cell: $(pkgcellfortest)\n$(notebook.cells[pkgcellfortest].code)"
    #     Pluto.save_notebook(notebook,"pluto-tmp.jl")
    #     input="pluto-tmp.jl"
    #     sleep(1)
    # end

"""
    run_notebook_in_current_environment(notebook)

Run Pluto notebook in current environment. This may be useful for testing
noteboos in CI. 
"""
function run_notebook_in_current_environment(notebookname)


    # Prevent Pluto from calling into registry update
    Pluto.PkgCompat._updated_registries_compat[]=true

    session = Pluto.ServerSession();
    session.options.server.disable_writing_notebook_files=true
    session.options.server.show_file_system=false
    session.options.server.launch_browser=false
    session.options.server.dismiss_update_notification=true
    session.options.evaluation.capture_stdout=false
    session.options.evaluation.workspace_use_distributed=false # this makes it fast

    wd=pwd()
    t= @elapsed notebook = Pluto.SessionActions.open(session, notebookname; run_async=false)
    @info "notebook executed in $(round(t,sigdigits=4)) seconds"
    cd(wd)
    errored=false
    for c in notebook.cells
        if c.errored
            errored=true
            @error "Error in  $(c.cell_id): $(c.output.body[:msg])\n $(c.code)"
        end
    end
    !errored
end

function run_all_tests(;run_notebooks=false)
    
    notebooks=["nbproto.jl",
               "flux-reconstruction.jl",
               "interfaces1d.jl",
               "problemcase.jl",
               "nonlinear-solvers.jl",
               "api-update.jl",
               ]


    @testset "basictest" begin
        run_tests_from_directory(@__DIR__,"test_")
    end
    @testset "Development Examples" begin
        run_tests_from_directory(joinpath(@__DIR__,"..","examples"),"Example0")
    end
    @testset "MultiD Examples" begin
         run_tests_from_directory(joinpath(@__DIR__,"..","examples"),"Example5")
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

    if run_notebooks && VERSION>v"1.8" && !(VERSION>v"1.9.99")
        notebookenv=joinpath(@__DIR__,"..","pluto-examples")
        Pkg.activate(notebookenv)
        Pkg.develop(path=joinpath(@__DIR__, ".."))
        Pkg.instantiate()
        ENV["PLUTO_PROJECT"]=notebookenv
        @testset "notebooks" begin
            for notebook in notebooks
                @info "notebook $(notebook):"
                @test run_notebook_in_current_environment(joinpath(@__DIR__,"..","pluto-examples",notebook))
                @info "notebook $(notebook) ok"
            end
        end
    end
    #    Pkg.activate(@__DIR__)
end

run_all_tests(;run_notebooks=true)
