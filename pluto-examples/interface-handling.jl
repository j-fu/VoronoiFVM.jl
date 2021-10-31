### A Pluto.jl notebook ###
# v0.17.0

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
    ENV["LC_NUMERIC"]="C" # prevent pyplot from messing up string2float conversion
	using Pkg
	using Triangulate
	using PlutoUI
	using PyPlot
	using Printf
	using SimplexGridFactory
	using GridVisualize
	using ExtendableGrids
	using VoronoiFVM
	GridVisualize.default_plotter!(PyPlot)
	PyPlot.svg(true)
end

# ╔═╡ cc49b312-7e21-49aa-b42d-32c184919aa2
TableOfContents(title="",depth=4)

# ╔═╡ b61dce01-f84f-4e46-be76-8561bf0dddaa
md"""
# Handling of interfaces in VoronoiFVM

Let ``\Omega=(-0.5,1.5)\times (0,2)=\Omega_1\cap \Omega_2`` be subdivided into two
subregions. and  ``\Gamma_{12}`` be the interface between them. Let ``Γ_1=(-0.5)\times (0,2)``,
        ``Γ_2=(-1.5)\times (0,2)``, and ``Γ_N=\partial\Omega \setminus (Γ_1\cap Γ_2)``.
We assume that it is triangulated with a Delaunay triangulation wich conforms to interior
and exterior boundaries. The later property enssures, that subdomain boundaries are aligned
with triangle edges, and  that for $i=1,2$ the circumcenters
of triangles from ``\Omega_i`` are situated  in ``\bar\Omega_i``
"""

# ╔═╡ f37e9de2-a681-45d8-9a0a-f53613f19319
md"""
Here we create such a domain using the Triangle mesh generator of J. Shewchuk via Julia's Triangulate.jl  package.
"""

# ╔═╡ b912c3ff-42aa-4484-8b84-f011c80382e9
begin
    Γ_1=1
    Γ_2=2
    Γ_12=3
    Γ_N=4
    Ω_1=1
    Ω_2=2
end;

# ╔═╡ 8b806c22-3f34-11eb-0e03-4d3fa6461629
function tworegiongrid(;minangle=20)
    triin=Triangulate.TriangulateIO()
    triin.pointlist=Cdouble[-0.5 0.0 ;0.2 0.0; 1.5 0.0 ; 1.5  2.0 ; 0.8 2.0; -0.5 2.0]'
    triin.segmentlist=Cint[1 2 ; 2 3 ; 3 4 ; 4 5 ; 5 6 ; 6 1 ; 2 5]'
    triin.segmentmarkerlist=Cint[Γ_N, Γ_N, Γ_2, Γ_N, Γ_N, Γ_1, Γ_12]
    angle=@sprintf("%.15f",minangle)
    triin.regionlist=Cdouble[0.2 0.2 Ω_1 0.18;
		             0.8 0.2 Ω_2 0.18]'
    (triout, vorout)=triangulate("paAq$(angle)DQv", triin)	
end

# ╔═╡ 40795d22-6091-4498-bd07-f716fb07682e
function plot_with_numbers(triout, vorout)
    PyPlot.clf()
    Triangulate.plot_triangulateio(PyPlot,triout,voronoi=vorout)
    PyPlot.scatter(triout.pointlist[1,:],triout.pointlist[2,:],color=:red,alpha=1)
    dxy=[-0.03,0.01]
    for ipoint=1:size(triout.pointlist,2)
	PyPlot.text((dxy.+triout.pointlist[:,ipoint])...,"$(ipoint)",ha=:right,color=:red)
    end
    for ipoint=1:size(vorout.pointlist,2)
	PyPlot.text((dxy.+vorout.pointlist[:,ipoint])...,"$(ipoint)",ha=:right,color=:green)		
    end
    PyPlot.gcf()
end

# ╔═╡ d55fdbfc-3f34-11eb-3f35-c3bcdb156288
triout, vorout=tworegiongrid()

# ╔═╡ 36ef21e0-078c-4287-a0fe-f942cd14dc08
plot_with_numbers(triout, vorout)

# ╔═╡ d6276724-6f55-426c-aba9-fb89d2c8166e
md"""
We see here collocation points in red, and triangle circumcenters aka Voronoi nodes in green.
"""

# ╔═╡ b4b335ef-de2b-4563-8f0f-cde4a04ac115
md"""
## Heterogeneous diffusion

### Problem formulations

```math
\begin{aligned}
-\nabla\cdot D\nabla u &=0 & \text{in}\quad Ω\\
                     D\nabla u\cdot \vec n + αu &=g_1& \text{on}\quad Γ_1 \\
                     D\nabla u\cdot \vec n + αu &=g_2& \text{on}\quad Γ_2\\
                     D\nabla u\cdot \vec n &=0& \text{on}\quad Γ_N
\end{aligned}
```

with ``D|_{\Omega_i}= D_i>0`` and ``α>0`` In fact, this is a sloppy version of a more precise formulation
for i=1,2: 
```math
\begin{aligned}
-\nabla\cdot D_i\nabla u_i &=0 & \text{in}\quad Ω_i\\
                     D\nabla u\cdot \vec n_i+ α u_i &= g_i& \text{on}\quad Γ_i \\
                     D_i\nabla u_i\cdot \vec n_i &=0& \text{on}\quad Γ_N\cap \partial Ω_i
\end{aligned}
```math
with the interface condition (in the sense of traces of functions)
```math
\begin{aligned}
    u_1 &= u_2  & \text{on} \quad  Γ_{12}\\
 D_1\nabla u_1\cdot \vec n_1 &=-D_2\nabla u_2\cdot \vec n_2 & \text{on} \quad Γ_{12}
\end{aligned}
```

The interface condition consists in the continuity of the solution and the continuity
of fluxes.

The weak formulation of this problem in ``V=H^1(\Omega)`` is: find ``u∈V`` such that
```math
    \sum_{i=1}^2 \int_{Ω_i} D_i\nabla u \nabla v dω + \sum_{i=1}^2 \int_{Γ_i} \alpha uv dγ =  \sum_{i=1}^2 \int_{Γ_i} \alpha g_iv dγ \quad ∀ v ∈ V
```

This implicitely contains the interface condition (continuity through ``u∈V`` and flux continuity due to 
Green's theorem (TODO: work this out).


A P1 finite element discretization is straightforward if the interface is aligned with triangle edges.
"""

# ╔═╡ d95afcb5-68f9-4d2e-9111-06fc38fe949a
md"""
### Finite volume approach

Regard $0=-\int_\omega{K} \nabla D\nabla u d\omega

#### Interior control volume, e.g. ``ω_{13}``
Let ``\sigma_{KL}=\partial\omega_K\cap \partial\omega_L``

```math
0= \sum_{N_K} \int_{\sigma_{KL}} D\nabla u \cdot \vec n_{KL} \approx  \sum_{N_K} D(u_K-u_L)\frac{|\sigma_{KL}|}{|h_{KL}|}
```

#### Boundary  control volume, e.g. ``ω_{8}``
Let ``\gamma_{Ki}=\partial\omega_K\cap Γ_i`` (i=1,2,N)

```math
\begin{aligned}
0=&\sum_{N_K} \int_{\sigma_{KL}} D\nabla u \cdot \vec n_{KL} + \sum_i  \int_{\gamma_{Ki}} D\nabla u \cdot \vec n_{i} \\
  \approx&  \sum_{N_K} D(u_K-u_L)\frac{|\sigma_{KL}|}{|h_{KL}|} +  |\gamma_{K_i}| (\alpha u_K-g_i)
\end{aligned}
```

#### Interface  control volume, e.g. ``ω_{7}``
Here, we come to the main "trick" (TODO: notations could be improved)
Let ``\omega_{Km}= \omega_K\cap \Omega_{m}`` for ``m=1,2`` and 
``\sigma_{KL}=\partial\omega_K\cap \partial\omega_L\cap \Omega_m``

Let Let ``\gamma_{K12}=\omega_K\cap Γ_{12}`` be the part of the interface intersecting ``\omega_K``


Then we have
```math
\begin{aligned}
    0  =& -\int_{\omega_K} -\nabla\cdot D\nabla u = -\int_{\omega_{K1}} -\nabla\cdot D_1\nabla u -\int_{\omega_{K2}} -\nabla\cdot D_2\nabla u\\
       =&  \sum_{N_K} \int_{\sigma_{KL1}} D_1\nabla u \cdot \vec n_{KL} + \int_{\gamma_{K12}} D_1\nabla u \cdot \vec n_1\\
        & + \sum_{N_K} \int_{\sigma_{KL2}} D_2\nabla u \cdot \vec n_{KL} + \int_{\gamma_{K12}} D_2\nabla u \cdot \vec n_2 \\
       = &  \sum_{N_K} \int_{\sigma_{KL1}} D_1\nabla u \cdot \vec n_{KL} + \sum_{N_K} \int_{\sigma_{KL2}} D_2\nabla u \cdot \vec n_{KL} \\
    \approx & \sum_{N_K} D(u_K-u_L)\frac{|\sigma_{KL1}|}{|h_{KL}|} + \sum_{N_K} D(u_K-u_L)\frac{|\sigma_{KL2}|}{|h_{KL}|}
\end{aligned}
```
The development used the flux continuity at the interface. Thus, interface control volumes are split into their
respective intersections with the adjacent domains, and integrals are performed separately. 

In practice, we attach the domain information to the triangles and perform the integration for each triangle separately.



"""

# ╔═╡ f9886dc9-5b1b-4c27-ac0d-5e65805e45bb
md"""
### Numerical solution with VoronoiFVM
"""

# ╔═╡ add94eb5-7859-474b-87fc-a4e4ab00e691
heterogrid=simplexgrid(triout.pointlist, 
	                  triout.trianglelist,
	Int32.(triout.triangleattributelist[1,:]),
	triout.segmentlist,
	triout.segmentmarkerlist)

# ╔═╡ 427504f5-6c74-40d5-a4b9-9d2d7bd4f540
ispec=1

# ╔═╡ d908722f-a92d-49de-a3c2-aae33cc5a2e4
D=[1,100]

# ╔═╡ 6f108820-4d06-4fa3-b6ad-8f7454494670
flux(y,u,edge)= y[ispec]=D[edge.region]*(u[ispec,1]-u[ispec,2])

# ╔═╡ fca05141-71cb-4ecc-83b9-f50836dda100
md"""
Here we apply the Robin and Neumann boundary condtions in a `breaction`
"""

# ╔═╡ a8a69c7c-810a-4ea5-9a72-e0a2eea49782
g=[0,1,0,0]

# ╔═╡ bdc68399-19eb-48b5-ba69-5745a2541a17
α=[100,100,0,0]

# ╔═╡ ba273c41-f622-4f1a-9e1c-48b4143a6d43
breaction(y,u,bnode)=y[ispec]=α[bnode.region]*u[ispec]-g[bnode.region]

# ╔═╡ 92dadccc-44e9-49e5-ac94-a224dd8b61aa
begin
	heterophysics=VoronoiFVM.Physics(flux=flux,breaction=breaction)	 
    heterosys=VoronoiFVM.System(heterogrid,heterophysics)
	enable_species!(heterosys,ispec,[1,2])
end

# ╔═╡ 40dc1b8d-9f11-4ae2-b566-793c3fea99fc
solution=solve(unknowns(heterosys,inival=0),heterosys)

# ╔═╡ ead2425e-0109-4798-b79c-cbc8f76c1414
scalarplot(heterogrid,solution[1,:],resolution=(300,300))

# ╔═╡ d5946f89-d8d7-4196-90a6-8e34723265c3
md"""
## Bulk reactions

In this case we have to evaluate integrals over reaction terms, which 
are split into the respective contributions from each subdomain which may
differ in the coefficients or in the physics.

TODO: elaborate
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
SimplexGridFactory = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
Triangulate = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
ExtendableGrids = "~0.8.7"
GridVisualize = "~0.3.9"
PlutoUI = "~0.7.16"
PyPlot = "~2.10.0"
SimplexGridFactory = "~0.5.9"
Triangulate = "~2.1.0"
VoronoiFVM = "~0.13.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-rc2"
manifest_format = "2.0"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "d9352737cef8525944bf9ef34392d756321cbd54"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.38"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Cassette]]
git-tree-sha1 = "6ce3cd755d4130d43bab24ea5181e77b89b51839"
uuid = "7057c7e9-c182-5462-911a-8362d720325c"
version = "0.3.9"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "3533f5a691e60601fe60c90d8bc47a27aa2907ec"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "dce3e3fea680869eaa0b774b2e8343e9ff442313"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.40.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "7220bc21c33e990c14f4a9a319b1d242ebc5b269"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "a0fcc1bb3c9ceaf07e1d0529c9806ce94be6adf9"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.9"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "Test"]
git-tree-sha1 = "1e8e50f054057f23e908fbd6935766dca6293cc2"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.8.7"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "0341e41e45e6c5e46be89c33543082a0a867707c"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "0.6.5"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "2db648b6712831ecb333eae76dbfd1c156ca13bb"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.2"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "8756f9935b7ccc9064c6eef0bff0ad643df733a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.7"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "63777916efbcb0ab6173d09a658fb7f2783de485"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.21"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "Requires", "StaticArrays"]
git-tree-sha1 = "925ba2f11df005d894b113292d32fca9afe3f8c8"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "0.3.9"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "5efcf53d798efede8fee5b2c8b09284be359bf24"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.2"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "f0c6489b12d28fb4c2103073ec7452f3423bd308"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1a8c6237e78b714e901e406c096fc8a65528af7d"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.1"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["DataStructures", "FileIO", "MacroTools", "Mmap", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "46b7834ec8165c541b0b5d1c8ba63ec940723ffb"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.15"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "6193c3815f13ba1b78a51ce391db8be016ae9214"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.4"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d911b6a12ba974dabe2291c6d450094a7226b372"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlutoUI]]
deps = ["Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "4c8a7d080daca18545c56f1cac28710c362478f3"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.16"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "4ba3651d33ef76e24fef6a598b63ffd1c5e1cd17"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.5"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "14c1b795b9d764e1784713941e787e1384268103"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.10.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[deps.RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "c944fa4adbb47be43376359811c0a14757bdc8a8"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.20.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimplexGridFactory]]
deps = ["DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GridVisualize", "LinearAlgebra", "Printf", "Test"]
git-tree-sha1 = "af52ec74a4b6cfcc5b6d60d259099fa0596de2c1"
uuid = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
version = "0.5.9"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "LightGraphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "999444c932d3bd32e7cb5a349e219bec9cfa2bd5"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.17.1"

[[deps.SparsityDetection]]
deps = ["Cassette", "DocStringExtensions", "LinearAlgebra", "SparseArrays", "SpecialFunctions"]
git-tree-sha1 = "9e182a311d169cb9fe0c6501aa252983215fe692"
uuid = "684fba80-ace3-11e9-3d08-3bc7ed6f96df"
version = "0.3.4"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "e7bc80dc93f50857a5d1e3c8121495852f407e6a"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "fed34d0e71b91734bf0a7e10eb1bb05296ddbcd0"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.Triangle_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bfdd9ef1004eb9d407af935a6f36a4e0af711369"
uuid = "5639c1d2-226c-5e70-8d55-b3095415a16a"
version = "1.6.1+0"

[[deps.Triangulate]]
deps = ["DocStringExtensions", "Libdl", "Printf", "Test", "Triangle_jll"]
git-tree-sha1 = "2b4f716b192c0c615d96d541ee029e85666388cb"
uuid = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
version = "2.1.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "e575cf85535c7c3292b4d89d89cc29e8c3098e47"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.1"

[[deps.VertexSafeGraphs]]
deps = ["LightGraphs"]
git-tree-sha1 = "b9b450c99a3ca1cc1c6836f560d8d887bcbe356e"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.1.2"

[[deps.VoronoiFVM]]
deps = ["DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "IterativeSolvers", "JLD2", "LinearAlgebra", "Printf", "RecursiveArrayTools", "SparseArrays", "SparseDiffTools", "SparsityDetection", "StaticArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "65f1d72aa575d3d6c348fcfc0d0d905fdbd7f5c3"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.13.2"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─cc49b312-7e21-49aa-b42d-32c184919aa2
# ╟─b61dce01-f84f-4e46-be76-8561bf0dddaa
# ╟─f37e9de2-a681-45d8-9a0a-f53613f19319
# ╠═b912c3ff-42aa-4484-8b84-f011c80382e9
# ╠═8b806c22-3f34-11eb-0e03-4d3fa6461629
# ╠═40795d22-6091-4498-bd07-f716fb07682e
# ╠═d55fdbfc-3f34-11eb-3f35-c3bcdb156288
# ╠═36ef21e0-078c-4287-a0fe-f942cd14dc08
# ╟─d6276724-6f55-426c-aba9-fb89d2c8166e
# ╟─b4b335ef-de2b-4563-8f0f-cde4a04ac115
# ╟─d95afcb5-68f9-4d2e-9111-06fc38fe949a
# ╟─f9886dc9-5b1b-4c27-ac0d-5e65805e45bb
# ╠═add94eb5-7859-474b-87fc-a4e4ab00e691
# ╠═427504f5-6c74-40d5-a4b9-9d2d7bd4f540
# ╠═d908722f-a92d-49de-a3c2-aae33cc5a2e4
# ╠═6f108820-4d06-4fa3-b6ad-8f7454494670
# ╟─fca05141-71cb-4ecc-83b9-f50836dda100
# ╠═a8a69c7c-810a-4ea5-9a72-e0a2eea49782
# ╠═bdc68399-19eb-48b5-ba69-5745a2541a17
# ╠═ba273c41-f622-4f1a-9e1c-48b4143a6d43
# ╠═92dadccc-44e9-49e5-ac94-a224dd8b61aa
# ╠═40dc1b8d-9f11-4ae2-b566-793c3fea99fc
# ╠═ead2425e-0109-4798-b79c-cbc8f76c1414
# ╟─d5946f89-d8d7-4196-90a6-8e34723265c3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
