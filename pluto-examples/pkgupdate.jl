#!/bin/sh
# -*-Julia-*-
#=
exec julia --startup-file=no --project=. "$0" "$@"
=#

using Pkg
Pkg.activate(@__DIR__)

using Pluto

notebooks = [
    "nbproto.jl",
    "api-update.jl",
    "flux-reconstruction.jl",
    "problemcase.jl",
    "nonlinear-solvers.jl",
]


tmp = mktempdir()
tmpjl = joinpath(tmp, "tmp.jl")

for notebook in notebooks
    println("Updating packages in $(notebook):")
    cp(joinpath(@__DIR__, notebook), tmpjl, force = true)
    Pluto.update_notebook_environment(tmpjl)
    cp(tmpjl, joinpath(@__DIR__, notebook), force = true)
    println("Updating of  $(notebook) done\n")
    Pkg.activate(@__DIR__)
end

dirs = ["test", "pluto-examples", "docs"]
for dir in dirs
    println("updating $(dir) environment")
    Pkg.activate(joinpath(@__DIR__, "..", dir))
    Pkg.status()
    Pkg.update()
    Pkg.status()
    Pkg.activate(@__DIR__)
end
