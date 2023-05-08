using Documenter, VoronoiFVM, Literate, PlutoSliderServer, ExtendableGrids, GridVisualize
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




function make_all(; with_examples = true, run_notebooks = true)

    generated_examples = []
    notebooks = []
    example_jl_dir = joinpath(@__DIR__, "..", "examples")
    example_md_dir = joinpath(@__DIR__, "src", "examples")
    notebook_html_dir = joinpath(@__DIR__, "src", "nbhtml")


    if with_examples
        #
        # Run notebooks
        #
        notebooks = [
            "API Updates" => "api-update.jl",
            "Nonlinear solver control" => "nonlinear-solvers.jl",
            "Obtaining vector fields" => "flux-reconstruction.jl",
            "Internal interfaces (1D)" => "interfaces1d.jl",
            "A case for caution" => "problemcase.jl",
        ]

        notebookjl = last.(notebooks)
        notebookmd = []


        # function rendernotebook(name)
        #     base=split(name,".")[1]
        #     input=joinpath(@__DIR__,"..","pluto-examples",base*".jl")
        #     output=joinpath(@__DIR__,"src","nbhtml",base*".html")
        #     session = Pluto.ServerSession();
        #     html_contents=PlutoStaticHTML.notebook2html(input;session)
        #     write(output, html_contents)
        # end


        # for notebook in notebookjl
        #     @info "Converting $(notebook)"
        #     rendernotebook(notebook)
        # end


        # Use sliderserver to generate html
        if run_notebooks
            export_directory(
                joinpath(@__DIR__, "..", "pluto-examples"),
                notebook_paths = notebookjl,
                Export_output_dir = joinpath(notebook_html_dir),
                Export_offer_binder = false,
            )
            
            # generate frame markdown for each notebook
            for notebook in notebookjl
                base = split(notebook, ".")[1]
                mdstring = """
                           ##### [$(base).jl](@id $(base))
                           [Download](https://github.com/j-fu/VoronoiFVM.jl/blob/master/pluto-examples/$(notebook))
                                       this [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook.
    
                           ```@raw html
                           <iframe style="height:20000px" width="100%" src="../$(base).html"> </iframe>
                           ```
                           """
                mdname = base * ".md"
                push!(notebookmd, joinpath("nbhtml", mdname))
                io = open(joinpath(notebook_html_dir, mdname), "w")
                write(io, mdstring)
                close(io)
            end
            
            notebooks = first.(notebooks) .=> notebookmd
            pushfirst!(notebooks, "About the notebooks" => "notebooks.md")
        end

        
        #
        # Generate Markdown pages from examples
        #
        
        for example_source in readdir(example_jl_dir)
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
        modules = [VoronoiFVM,VoronoiFVM.SolverStrategies],
        checkdocs = :all,
        clean = false,
        doctest = true,
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
            ],
            "Tutorial Notebooks" => notebooks,
            "Examples" => generated_examples,
        ],
    )

    with_examples && rm(example_md_dir, recursive = true)
    run_notebooks && rm(notebook_html_dir, recursive = true)

    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/VoronoiFVM.jl.git")
    end
end

make_all(with_examples = true, run_notebooks = true)
#make_all(with_examples=true,run_notebooks=false)
