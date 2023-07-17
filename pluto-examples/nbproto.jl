### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ b285aca3-dee5-4b77-9276-537563e8643b
begin
    if initialized
        using VoronoiFVM
        using ExtendableGrids
        using Test
        using PlutoUI
        using GridVisualize
        using PlutoVista
        GridVisualize.default_plotter!(PlutoVista)
    end
end;

# ╔═╡ 4ed0c302-26e4-468a-a40d-0e6406f802d0
md"""
# Intro
"""

# ╔═╡ 7a104243-d3b9-421a-b494-5607c494b106
TableOfContents(; aside = false)

# ╔═╡ c8eda836-d719-4412-895e-c3a24fec21ec
scalarplot(sin.(0:0.1:10), size = (500, 200))

# ╔═╡ 3eef08af-f6ba-4874-82c0-65ff53e7f7da
@test true

# ╔═╡ fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
html"""<hr><hr><hr>"""

# ╔═╡ f50c6497-cba3-491a-bedd-5f94f88f76fb
md"""
# Appendix: Tests & Development
"""

# ╔═╡ ad899a81-baab-4433-8b7f-1e5c3b18dae6
md"""
The next cell activates a development environment if the notebook is loaded from a checked out VoronoiFVM.jl
and the environment variable `PLUTO_DEVELOP` is set, e.g. during continuous integration tests.
Otherwise, Pluto's built-in package manager is used.
"""

# ╔═╡ e00d0175-866e-4f0f-8121-49e7bbda6fb6
begin
    import Pkg as _Pkg
    if isfile(joinpath(@__DIR__, "..", "src", "VoronoiFVM.jl")) && haskey(ENV,"PLUTO_DEVELOP")
        _Pkg.activate(@__DIR__)
        _Pkg.instantiate()
        _Pkg.develop(path=joinpath(@__DIR__, ".."))
        using Revise
    end
    initialized = true
end;

# ╔═╡ bdbe6513-70b1-4d97-a79c-71534caad2b7


# ╔═╡ Cell order:
# ╠═b285aca3-dee5-4b77-9276-537563e8643b
# ╟─4ed0c302-26e4-468a-a40d-0e6406f802d0
# ╟─7a104243-d3b9-421a-b494-5607c494b106
# ╠═c8eda836-d719-4412-895e-c3a24fec21ec
# ╠═3eef08af-f6ba-4874-82c0-65ff53e7f7da
# ╟─fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
# ╟─f50c6497-cba3-491a-bedd-5f94f88f76fb
# ╟─ad899a81-baab-4433-8b7f-1e5c3b18dae6
# ╠═e00d0175-866e-4f0f-8121-49e7bbda6fb6
# ╟─bdbe6513-70b1-4d97-a79c-71534caad2b7
