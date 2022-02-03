### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ e00d0175-866e-4f0f-8121-49e7bbda6fb6
begin
    using Pkg
    inpluto=isdefined(Main,:PlutoRunner)
    developing=false	
    if inpluto && isfile(joinpath(@__DIR__,"..","src","VoronoiFVM.jl"))
	# We try to outsmart Pluto's cell parser here.
	# This activates an environment in VoronoiFVM/pluto-examples
	eval(:(Pkg.activate(joinpath(@__DIR__))))
    eval(:(Pkg.instantiate()))
	# use Revise if we develop VoronoiFVM
	using Revise
	# This activates the checked out version of VoronoiFVM.jl for development
	eval(:(Pkg.develop(path=joinpath(@__DIR__,".."))))
	developing=true
    end
end;

# ╔═╡ b285aca3-dee5-4b77-9276-537563e8643b
begin 
    using VoronoiFVM
    using ExtendableGrids
    using Test
    if inpluto
 	using PlutoUI
	using LinearAlgebra
	using GridVisualize
	using PlutoVista
	GridVisualize.default_plotter!(PlutoVista)
    end
end;

# ╔═╡ 4ed0c302-26e4-468a-a40d-0e6406f802d0
md"""
# Nonlinear solver control
"""

# ╔═╡ eb9ea477-6122-4774-b4a4-04dd7346e2b6
md"""
Generally, nonlinear systems in this package  are solved using Newton's method.  In many cases, the default settings provided by this package work well. However, the convergence of Newton's method is only guaranteed with initial values s7ufficiently close to the exact solution. This notebook describes how change the default settings for the solution of nonlinear problems with VoronoiFVM.jl. 
"""

# ╔═╡ 7a104243-d3b9-421a-b494-5607c494b106
inpluto && TableOfContents(;aside=false)

# ╔═╡ 9843b65c-6ca8-4ef8-a896-2c8cec4bff7c
md"""
Define a nonlinear Poisson equation to have an example. Let ``Ω=(0,10)`` and define
```math
\begin{aligned}
-Δ u + e^u-e^{-u} & = 0 & \text{in}\; Ω \\
	u(0)&=100\\
    u(10)&=0
\end{aligned}
```
"""

# ╔═╡ a70b8d85-66aa-4da8-8157-dd0244e3e4f6
X=0:0.001:1

# ╔═╡ c8eda836-d719-4412-895e-c3a24fec21ec
flux(y,u,edge)=y[1]=u[1,1]-u[1,2];

# ╔═╡ c09f5dfc-fc47-4952-8051-54731ec2b00b
function reaction(y,u,node)
	eplus=exp(u[1])
	eminus=1/eplus
	y[1]=eplus-eminus
end

# ╔═╡ eab04557-5084-4174-b275-b4d4399238e5
function bc(y,u,node)
	boundary_dirichlet!(y,u,node; region=1,value=100)
	boundary_dirichlet!(y,u,node; region=2,value=0.0)
end;

# ╔═╡ 316112fd-6553-494a-8e4a-65b34829891d
system=VoronoiFVM.System(X;
						flux=flux,
						reaction=reaction,
						bcondition=bc,
						species=1);

# ╔═╡ 42e54ff9-fc11-4a31-ba10-32f8d817d50c
md"""
## Solution using default settings
"""

# ╔═╡ 050ed807-1bca-4015-85f3-d7347ecb7e6b
begin
	sol=solve(system; log=true)
	hist=history(system)
end;

# ╔═╡ b9bb8020-5470-4964-818c-7f9b3bf2a8b4
inpluto && scalarplot(system,sol,resolution=(500,200),
        xlabel="x",ylable="y", title="solution")

# ╔═╡ b3124c06-1f40-46f5-abee-9c2e8e538162
md"""
With `log=true`, the `solve` method in addition to the solution records the solution
history which after finished solution can be obtatined as `history(system)`.
"""

# ╔═╡ 973db266-eb91-46e8-a917-9beeeb2c1ea7
md"""
The history can be plotted:
"""

# ╔═╡ 20e925f3-43fa-4db1-a656-79cf9c1c3303
plothistory(h)=scalarplot(1:length(h),h,resolution=(500,200),
             yscale=:log,
	        xlabel="step",
            ylabel="||δu||_∞",
            title= "Maximum norm of Newton update");

# ╔═╡ ebdc2c82-f72e-4e35-a63f-4ba5154e294f
inpluto && plothistory(hist)

# ╔═╡ 3d49aafd-79a6-4e2f-b422-24a5be7aa87a
md"""
History can be summarized:
"""

# ╔═╡ a217a308-5569-41e4-9d9d-418217017030
summary(hist)

# ╔═╡ f951b78c-a3fb-432e-bf6e-b956049d6a0d
md"""
History can be explored in detail:
"""

# ╔═╡ fcd7beb6-51b2-4ca8-a184-43ba6b5d2c1a
details(hist)

# ╔═╡ baed6e43-747b-4557-95c3-d4805f12b3a1
md"""
With default solver settings, for this particular problem, Newton's method needs $(length(history(system))) iteration steps.
"""

# ╔═╡ ccef0590-d5f8-4ee2-bb7a-d48ccfbd4d99
check(sol)= isapprox(sum(sol),2554.7106586964906,rtol=1.0e-12)

# ╔═╡ c0432a54-85ec-4478-bd75-f5b43770a117
@test check(sol)

# ╔═╡ a4c1a2f5-dddb-4a45-bc03-e168cbd7d569
md"""
## Damping 
"""

# ╔═╡ 38539474-af65-4f0f-9aa1-2292f4f6331c
md"""
Try to use a damped version of Newton method. The damping scheme is rather simple: an initial damping value `damp_initial` is increased by a growth factor `damp_growth` in each iteration until it reaches 1.
"""

# ╔═╡ d961d026-0b55-46c2-8555-8ef0763d8016
begin
  sol1=solve(system,log=true,inival=1,damp=0.15,damp_grow=1.5) 
  hist1=history(system)
end

# ╔═╡ e66d20f0-4987-471b-82ee-0c56160f9b01
inpluto && plothistory(history(system))

# ╔═╡ 35971019-fa07-4033-aebf-7872030a0cef
details(hist1)

# ╔═╡ 63ce84fc-e81b-4768-8122-36bfbd789727
summary(hist1)

# ╔═╡ bf3fe305-eecc-413e-956c-9737b9160f83
md"""
We see that the number of iterations decreased significantly.
"""

# ╔═╡ c8227ea2-2189-438f-bd02-e9d803031830
@test check(sol1)

# ╔═╡ f8cf5cdb-647c-4eb8-8bde-b12844f72b24
md"""
## Embedding
"""

# ╔═╡ eea852d0-937e-4af0-8688-f941e5d31697
md"""
Another possibility is the embedding (or homotopy) via a parameter: start with solving a simple problem and increase the level of complexity by increasing the parameter until the full problem is solved. This process is controlled by the parameters 
- `Δp`: initial parameter step size
- `Δp_min`: minimal parameter step size
- `Δp_max`: maximum parameter step size
- `Δp_grow`: maximum growth factor
- `Δu_opt`: optimal difference of solutions between two embedding steps

After successful solution of a parameter, the new parameter step size is calculated as
`Δp_new=min(Δp_max, Δp_grow, Δp*Δu_opt/(|u-u_old|+1.0e-14))` and adjusted to the end of the parameter interval. 

If the solution is unsuccessful, the parameter stepsize is halved and solution is retried, until the minimum step size is reached. 
"""

# ╔═╡ a71cbcd4-310e-47a8-94f9-1159995a7711
function pbc(y,u,node)
	boundary_dirichlet!(y,u,node; region=1,value=100*embedparam(node))
	boundary_dirichlet!(y,u,node; region=2,value=0)
end;

# ╔═╡ 89435c65-0520-4430-8727-9d013df6182d
system2=VoronoiFVM.System(X;
							flux=flux,
							reaction= function(y,u,node) 
								         reaction(y,u,node)
								         y[1]=y[1]*embedparam(node)
							          end,
							bcondition=pbc,
							species=1);

# ╔═╡ cb382145-c4f1-4222-aed7-32fa1e3bd7e4
begin sol2=solve(system2, inival=0,
					log=true,
					embed=(0,1),
					Δp=0.1,
	                max_lureuse=0,
					Δp_grow=1.2,
					Δu_opt=15)
	history2=history(system2)
end

# ╔═╡ a0b2aaf5-f7b1-40eb-ac4e-9790a8bbf09d
summary(history2)

# ╔═╡ 0d2311aa-f79a-4a44-bac4-59d6c5457ca5
inpluto && plothistory(vcat(history2...))

# ╔═╡ 75ab714d-0251-42ef-a415-ba6bed4c688f
@test check(sol2[end])

# ╔═╡ fe22e2d8-fd70-49e2-8baf-d5ec3faead24
md"""
For this particular problem, embedding uses less overall Newton steps than the default settings, but the damped method is faster.
"""

# ╔═╡ 2cb7f0f4-7635-462c-89f8-93d18e7f24fb
md"""
## Solver control
"""

# ╔═╡ 5f13831d-b73f-44fb-a870-16261f926ed5
md"""
Here we show the docsctring of `SolverControl` (formerly `NewtonControl`). This is a struct which can be passed to the `solve` method. Alternatively, as shown in this notebook, keyword arguments named like its entries can be passed directly to the `solve` method.
"""

# ╔═╡ 038d096a-1339-403c-aa9c-3112442d622d
@doc SolverControl

# ╔═╡ fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
html"""<hr><hr><hr>"""

# ╔═╡ f50c6497-cba3-491a-bedd-5f94f88f76fb
md"""
# Appendix: Tests & Development
"""

# ╔═╡ ad899a81-baab-4433-8b7f-1e5c3b18dae6
md"""
This notebook is also run during the automatic unit tests. In this case, all interactive elements and visualizations should be deactivated.
For this purposes, the next cell detects if the notebook is running under Pluto
and sets the `inpluto` flag accordingly.

Furthermore, the cell activates a development environment if the notebook is loaded from a checked out VoronoiFVM.jl. Otherwise, Pluto's built-in package manager is used.
"""

# ╔═╡ bdbe6513-70b1-4d97-a79c-71534caad2b7
if developing 
	md""" Developing VoronoiFVM at  $(pathof(VoronoiFVM))"""
else
	md""" Loaded VoronoiFVM from  $(pathof(VoronoiFVM))"""
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
ExtendableGrids = "~0.8.11"
GridVisualize = "~0.4.3"
PlutoUI = "~0.7.23"
PlutoVista = "~0.8.12"
Revise = "~3.2.0"
VoronoiFVM = "~0.14.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "265b06e2b1f6a216e0e8f183d28e4d354eab3220"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Cassette]]
git-tree-sha1 = "6ce3cd755d4130d43bab24ea5181e77b89b51839"
uuid = "7057c7e9-c182-5462-911a-8362d720325c"
version = "0.3.9"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "9bc5dac3c8b6706b58ad5ce24cffd9861f07c94f"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.9.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "a0fcc1bb3c9ceaf07e1d0529c9806ce94be6adf9"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.9"

[[ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "Test"]
git-tree-sha1 = "fbb0efd29f2ba5e25eeaf73b76257acfc1a28630"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.8.11"

[[ExtendableSparse]]
deps = ["DocStringExtensions", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "793bd32bb280668e80c476ce4a3d0f171c8122d5"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "0.6.6"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "2db648b6712831ecb333eae76dbfd1c156ca13bb"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.2"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2b72a5624e289ee18256111657663721d59c143e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.24"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "92243c07e786ea3458532e199eb3feee0e7e08eb"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.4.1"

[[GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "eef34bda67d8ea865d7467a141b2c7d8eccd5592"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "0.4.4"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "5335c4c9a30b4b823d1776d2db09882cbfac9f1e"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.16"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e273807f38074f033d94207a201e6e827d8417db"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.8.21"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "491a883c4fef1103077a7f648961adbf9c8dd933"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.1.2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "f755f36b19a5116bb580de457cda0c140153f283"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.6"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d7fa6237da8004be601e19bd6666083056649918"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.3"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "fed057115644d04fba7f4d768faeeeff6ad11a60"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.27"

[[PlutoVista]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "GridVisualize", "HypertextLiteral", "UUIDs"]
git-tree-sha1 = "2435d1d3e02db324414f268f30999b5c06a0d10f"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "0.8.12"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "39aa6d41a8166be535d8b40562f89387131dc4ff"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.21.0"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "8f82019e525f4d5c669692772a6f4b0a58b06a6a"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.2.0"

[[Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "e55f4c73ec827f96cd52db0bc6916a3891c726b5"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.2.1"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "f33a0f6721b270cdf417f0c986e93c973e5913c8"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.19.1"

[[SparsityDetection]]
deps = ["Cassette", "DocStringExtensions", "LinearAlgebra", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "9e182a311d169cb9fe0c6501aa252983215fe692"
uuid = "684fba80-ace3-11e9-3d08-3bc7ed6f96df"
version = "0.3.4"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[VoronoiFVM]]
deps = ["DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "IterativeSolvers", "JLD2", "LinearAlgebra", "Printf", "RecursiveArrayTools", "SparseArrays", "SparseDiffTools", "SparsityDetection", "StaticArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "0cf39adeb43883e0d33bde1c5a2afdbd34cbf330"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.14.0"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─4ed0c302-26e4-468a-a40d-0e6406f802d0
# ╟─eb9ea477-6122-4774-b4a4-04dd7346e2b6
# ╟─7a104243-d3b9-421a-b494-5607c494b106
# ╠═b285aca3-dee5-4b77-9276-537563e8643b
# ╟─9843b65c-6ca8-4ef8-a896-2c8cec4bff7c
# ╠═a70b8d85-66aa-4da8-8157-dd0244e3e4f6
# ╠═c8eda836-d719-4412-895e-c3a24fec21ec
# ╠═c09f5dfc-fc47-4952-8051-54731ec2b00b
# ╠═eab04557-5084-4174-b275-b4d4399238e5
# ╠═316112fd-6553-494a-8e4a-65b34829891d
# ╟─42e54ff9-fc11-4a31-ba10-32f8d817d50c
# ╠═050ed807-1bca-4015-85f3-d7347ecb7e6b
# ╟─b9bb8020-5470-4964-818c-7f9b3bf2a8b4
# ╟─b3124c06-1f40-46f5-abee-9c2e8e538162
# ╟─973db266-eb91-46e8-a917-9beeeb2c1ea7
# ╠═20e925f3-43fa-4db1-a656-79cf9c1c3303
# ╠═ebdc2c82-f72e-4e35-a63f-4ba5154e294f
# ╟─3d49aafd-79a6-4e2f-b422-24a5be7aa87a
# ╠═a217a308-5569-41e4-9d9d-418217017030
# ╟─f951b78c-a3fb-432e-bf6e-b956049d6a0d
# ╟─fcd7beb6-51b2-4ca8-a184-43ba6b5d2c1a
# ╟─baed6e43-747b-4557-95c3-d4805f12b3a1
# ╠═ccef0590-d5f8-4ee2-bb7a-d48ccfbd4d99
# ╠═c0432a54-85ec-4478-bd75-f5b43770a117
# ╟─a4c1a2f5-dddb-4a45-bc03-e168cbd7d569
# ╟─38539474-af65-4f0f-9aa1-2292f4f6331c
# ╠═d961d026-0b55-46c2-8555-8ef0763d8016
# ╟─e66d20f0-4987-471b-82ee-0c56160f9b01
# ╠═35971019-fa07-4033-aebf-7872030a0cef
# ╠═63ce84fc-e81b-4768-8122-36bfbd789727
# ╟─bf3fe305-eecc-413e-956c-9737b9160f83
# ╠═c8227ea2-2189-438f-bd02-e9d803031830
# ╟─f8cf5cdb-647c-4eb8-8bde-b12844f72b24
# ╟─eea852d0-937e-4af0-8688-f941e5d31697
# ╠═a71cbcd4-310e-47a8-94f9-1159995a7711
# ╠═89435c65-0520-4430-8727-9d013df6182d
# ╠═cb382145-c4f1-4222-aed7-32fa1e3bd7e4
# ╠═a0b2aaf5-f7b1-40eb-ac4e-9790a8bbf09d
# ╟─0d2311aa-f79a-4a44-bac4-59d6c5457ca5
# ╠═75ab714d-0251-42ef-a415-ba6bed4c688f
# ╟─fe22e2d8-fd70-49e2-8baf-d5ec3faead24
# ╟─2cb7f0f4-7635-462c-89f8-93d18e7f24fb
# ╟─5f13831d-b73f-44fb-a870-16261f926ed5
# ╠═038d096a-1339-403c-aa9c-3112442d622d
# ╟─fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
# ╟─f50c6497-cba3-491a-bedd-5f94f88f76fb
# ╟─ad899a81-baab-4433-8b7f-1e5c3b18dae6
# ╠═e00d0175-866e-4f0f-8121-49e7bbda6fb6
# ╟─bdbe6513-70b1-4d97-a79c-71534caad2b7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
