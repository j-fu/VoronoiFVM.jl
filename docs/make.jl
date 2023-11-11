push!(LOAD_PATH, joinpath(@__DIR__, ".."))
using Documenter, VoronoiFVM, Literate, PlutoStaticHTML, ExtendableGrids, GridVisualize, Pkg
using LinearAlgebra

#
# Replace SOURCE_URL marker with github url of source
#
function replace_source_url(input, source_url)
    lines_in = collect(eachline(IOBuffer(input)))
    lines_out = IOBuffer()
    for line in lines_in
        println(lines_out, replace(line, "SOURCE_URL" => source_url))
    end
    return String(take!(lines_out))
end


function make_all(; with_examples = true, with_notebooks = true, example = nothing)

    generated_examples = []
    notebooks = []
    example_jl_dir = joinpath(@__DIR__, "..", "examples")
    example_md_dir = joinpath(@__DIR__, "src", "examples")
    notebook_md_dir = joinpath(@__DIR__, "src", "notebooks")

    with_examples && rm(example_md_dir, recursive = true, force = true)

    if with_notebooks
        rm(notebook_md_dir, recursive = true, force = true)
        mkdir(notebook_md_dir)
    end

    if with_notebooks
        thisdir = pwd()
        notebookdir = joinpath(@__DIR__, "..", "pluto-examples")
        Pkg.activate(notebookdir)
        Pkg.develop(path = joinpath(@__DIR__, ".."))
        Pkg.instantiate()
        Pkg.activate(thisdir)
        ENV["PLUTO_PROJECT"] = notebookdir
        notebooks =
            [
                "Outflow boundary conditions" => "outflow.jl",
                "Obtaining vector fields" => "flux-reconstruction.jl",
                "Internal interfaces (1D)" => "interfaces1d.jl",
                "A case for caution" => "problemcase.jl",
                "Nonlinear solver control" => "nonlinear-solvers.jl",
                "API Updates" => "api-update.jl",
            ]e
        oopts = OutputOptions(; append_build_context = true)
        output_format = documenter_output
        bopts = BuildOptions(notebookdir; output_format)

        notebookjl = last.(notebooks)
        notebookmd = [split(notebook, ".")[1] * ".md" for notebook in notebookjl]
        build_notebooks(bopts, notebookjl, oopts)
        for nb in notebookmd
            mv(joinpath(notebookdir, nb), joinpath(notebook_md_dir, nb))
        end

        notebooks = first.(notebooks) .=> joinpath.("notebooks", notebookmd)
        pushfirst!(notebooks, "About the notebooks" => "notebooks.md")
        @show notebooks
    end


    if with_examples
        #
        # Generate Markdown pages from examples
        #
        if example == nothing
            example_sources = readdir(example_jl_dir)
        else
            example_sources = [example]
        end
        for example_source in example_sources
            base, ext = splitext(example_source)
            if ext == ".jl"
                source_url =
                    "https://github.com/j-fu/VoronoiFVM.jl/raw/master/examples/" *
                    example_source
                preprocess(buffer) = replace_source_url(buffer, source_url)
                Literate.markdown(
                    joinpath(@__DIR__, "..", "examples", example_source),
                    example_md_dir,
                    documenter = false,
                    info = false,
                    preprocess = preprocess,
                )
            end
        end


        generated_examples =
            vcat(["runexamples.md"], joinpath.("examples", readdir(example_md_dir)))
    end

    makedocs(
        sitename = "VoronoiFVM.jl",
        modules = [VoronoiFVM, VoronoiFVM.SolverStrategies],
        checkdocs = :all,
        clean = false,
        doctest = true,
        warnonly = true,
        authors = "J. Fuhrmann",
        repo = "https://github.com/j-fu/VoronoiFVM.jl",
        pages = [
            "Home" => "index.md",
            "changes.md",
            "method.md",
            "API Documentation" => [
                "system.md",
                "physics.md",
                "solutions.md",
                "solver.md",
                "post.md",
                "quantities.md",
                "misc.md",
                "internal.md",
                "allindex.md",
                "devel.md",
            ],
            "Tutorial Notebooks" => notebooks,
            "Examples" => generated_examples,
        ],
    )

    with_examples && rm(example_md_dir, recursive = true, force = true)
    with_notebooks && rm(notebook_md_dir, recursive = true, force = true)

    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/VoronoiFVM.jl.git")
    end
end

if isinteractive()
    make_all(with_examples = false, with_notebooks = false)
else
    make_all(with_examples = true, with_notebooks = true)
end

