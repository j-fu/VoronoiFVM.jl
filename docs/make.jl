using Documenter, ExampleJuggler, VoronoiFVM
using ExtendableGrids, GridVisualize, LinearAlgebra, OrdinaryDiffEq, RecursiveArrayTools, SciMLBase

function make_all(; with_examples = true, with_notebooks = true, example = nothing)
    ExampleJuggler.verbose!(true)

    cleanexamples()
    notebookdir = joinpath(@__DIR__, "..", "pluto-examples")
    exampledir = joinpath(@__DIR__, "..", "examples")

    notebooks = [
        "OrdinaryDiffEq.jl nonlinear diffusion" =>   "ode-diffusion1d.jl",
        "Outflow boundary conditions" => "outflow.jl",
        "Obtaining vector fields" => "flux-reconstruction.jl",
        "Internal interfaces (1D)" => "interfaces1d.jl",
        "A case for caution" => "problemcase.jl",
        "Nonlinear solver control" => "nonlinear-solvers.jl",
        "API Updates" => "api-update.jl",
    ]
   notebook_examples = @docplutonotebooks(notebookdir, notebooks, iframe=false)
   notebook_examples = vcat(["About the notebooks" => "notebooks.md"], notebook_examples)

    modules = filter(ex -> splitext(ex)[2] == ".jl", basename.(readdir(exampledir)))
    module_examples = @docmodules(exampledir, modules)
    module_examples = vcat(["About the examples" => "runexamples.md"], module_examples)

    makedocs(; sitename = "VoronoiFVM.jl",
             modules = [VoronoiFVM, VoronoiFVM.SolverStrategies],
             checkdocs = :all,
             clean = false,
             doctest = false,
             warnonly = true,
             authors = "J. Fuhrmann",
             repo = "https://github.com/j-fu/VoronoiFVM.jl",
             format = Documenter.HTML(; size_threshold_ignore = last.(notebook_examples),
                                      mathengine = MathJax3()),
            pages = [
                 "Home" => "index.md",
                 "changes.md",
                 "method.md",
                 "API Documentation" => [
                     "system.md",
                     "physics.md",
                     "solutions.md",
                     "solver.md",
                     "odesolver.md",
                     "post.md",
                     "quantities.md",
                     "misc.md",
                     "internal.md",
                     "allindex.md",
                     "devel.md",
                 ],
                "Tutorial Notebooks" => notebook_examples,
                 "Examples" => module_examples,
             ])

    cleanexamples()

    if !isinteractive()
        deploydocs(; repo = "github.com/j-fu/VoronoiFVM.jl.git")
    end
end

if isinteractive()
    make_all(; with_examples = false, with_notebooks = false)
else
    make_all(; with_examples = true, with_notebooks = true)
end
