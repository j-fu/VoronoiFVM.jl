### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e00d0175-866e-4f0f-8121-49e7bbda6fb6
begin
    import Pkg as _Pkg
    developing = false
    if isfile(joinpath(@__DIR__, "..", "src", "VoronoiFVM.jl"))
        _Pkg.activate(@__DIR__)
        _Pkg.instantiate()
        using Revise
        developing = true
    end
    initialized = true
end;

# ╔═╡ b285aca3-dee5-4b77-9276-537563e8643b
begin
    if initialized
        using VoronoiFVM
        using ExtendableGrids
        using Test
        using PlutoUI
        using GridVisualize
        using PlutoVista
		using LinearSolve
		using ILUZero
		using LinearAlgebra
        GridVisualize.default_plotter!(PlutoVista)
    end
end;

# ╔═╡ 4ed0c302-26e4-468a-a40d-0e6406f802d0
md"""
# API updates 
"""

# ╔═╡ 3e6b4ffa-7b33-4b94-9fd6-75b030d5a115
md"""
Here we describe some updates for the API of `VoronoiFVM.jl`. These have been implemented mostly on top of the existing API, whose functionality is not affected.
"""

# ╔═╡ 7a104243-d3b9-421a-b494-5607c494b106
TableOfContents(; aside = false,depth=5)

# ╔═╡ a2f1e6ba-80b2-4902-9c5d-2806a7fb16f6
md"""
## v0.19
"""

# ╔═╡ 56976316-ca4f-4d5c-b303-026edd8751c2
md"""
This is a breaking release. Implementations using default solver settings should continue to work (albeit possibly with deprecation and 
allocation warnings). Really breaking is control of iterative linear solvers and allocation checks.

"""

# ╔═╡ 740dd980-ca12-4f51-88b2-408930704952
md"""
### Solve now a method of CommonSolve.solve
"""

# ╔═╡ db216ebf-1869-4b0f-8dbd-2a8a244c3e4c
md"""
As a consequence, all VoronoiFVM.solve methods with signatures others than `solve(system; kwargs...)`  are now deprecated
"""

# ╔═╡ 9eae0139-379d-47ce-a8d2-932c68504317
n=100

# ╔═╡ 8a31ed03-9cd8-4532-9fb9-8c060f777693
begin
	h = 1.0 / convert(Float64, n)
    const eps = 1.0e-2
    function reaction(f, u, node)
        f[1] = u[1]^2
    end

    function flux(f, u, edge)
        f[1] = eps * (u[1, 1]^2 - u[1, 2]^2)
    end

    function source(f, node)
        x1 = node[1] - 0.5
        x2 = node[2] - 0.5
        f[1] = exp(-20.0 * (x1^2 + x2^2))
    end

    function storage(f, u, node)
        f[1] = u[1]
    end

    function bcondition(f, u, node)
        boundary_dirichlet!(
            f,
            u,
            node,
            species = 1,
            region = 2,
            value = ramp(node.time, dt = (0, 0.1), du = (0, 1)),
        )
        boundary_dirichlet!(
            f,
            u,
            node,
            species = 1,
            region = 4,
            value = ramp(node.time, dt = (0, 0.1), du = (0, 1)),
        )
    end

    sys0 =
        VoronoiFVM.System(0.0:h:1.0,0.0:h:1.0; reaction, flux, source, storage, bcondition, species = [1])
end

# ╔═╡ 4279ae2e-a948-4358-9037-0c6895ecb809
md"""
Deprecated call:
"""

# ╔═╡ 58cf3e44-ee53-42f7-9132-eacfd900dc3a
begin
	inival=unknowns(sys0;inival=0.1)
	sol00=unknowns(sys0)
	solve!(sol00,inival,sys0)
end

# ╔═╡ b499a195-f06d-400f-9407-b07c3212c095
md"""
Replace this by:
"""

# ╔═╡ f660fc6e-0ab5-4fae-918f-39bf0c2153e5
sol0=solve(sys0;inival=0.1)

# ╔═╡ 13f284e8-4fe9-42ce-873a-c65febc7d7df
md"""
#### Docstring of solve
"""

# ╔═╡ c53cb54f-6099-4d09-9797-3da1f5428586
(@doc solve).content[end]

# ╔═╡ 2cebb072-bef1-4bce-9cbb-6b43b877b28b
md"""
#### Docstring of SolverControl
"""

# ╔═╡ 7b0f2021-a2fd-4bb2-a23a-432f61a38a07
@doc SolverControl

# ╔═╡ 182a9a9f-3f2c-4df3-abad-08488bb9fb33
md"""
### Rely on LinearSolve.jl for linear system solution
"""

# ╔═╡ 381c464e-a0cc-46ea-b0ba-089250c31555
md"""
This provides easy access to a large variety of linear solvers:
"""

# ╔═╡ ba4b8b92-617a-40b4-b3c9-1367067027fa
md"""
#### LU factorization from UMFPACK
"""

# ╔═╡ 00634dc4-e9c7-4165-a9ac-f3ffc8007e76
umf_sol=solve(sys0; inival=0.1, method_linear=UMFPACKFactorization(),verbose=true)

# ╔═╡ d6c88c3e-6b8a-45e9-9344-71f0de3fff51
@test isapprox(umf_sol, sol0, atol=1.0e-7)

# ╔═╡ cb2d42c4-8d5d-465c-9376-f568588ab453
md"""
#### LU factorization from Sparspak.jl
"""

# ╔═╡ 5908ae1d-b3b5-4681-a0b8-080f052af40f
sppk_sol=solve(sys0; inival=0.1, method_linear=SparspakFactorization(),verbose=true)

# ╔═╡ 22c3cecc-12d2-4f7d-8a53-894a5ea513f0
@test isapprox(sppk_sol, sol0, atol=1.0e-7)

# ╔═╡ 3107c644-5f43-4322-aad2-04a2ab0d576e
md"""
#### Iterative solvers
"""

# ╔═╡ 7e94ff44-4807-41b0-875d-526390903942
md"""
##### BICGstab from Krylov.jl with diagonal preconditioner 
"""

# ╔═╡ b70524fc-b8b1-4dee-b77d-e3f8d6d2837b
krydiag_sol=solve(sys0; inival=0.1, 
method_linear=KrylovJL_BICGSTAB(),
precon_linear=A->Diagonal(diag(A)),verbose=true)

# ╔═╡ d21d3236-b3d7-4cf8-ab9d-b0be44c9970b
@test isapprox(krydiag_sol, sol0, atol=1.0e-6)

# ╔═╡ f935061c-71db-4aaf-864b-bb566e68643e
md"""
##### BICGstab from Krylov.jl with delayed factorization preconditoner
"""

# ╔═╡ 4c7f8bbf-40c4-45b9-a62e-99ffaae30af1
    krydel_sol = solve(
        sys0;
        inival = 0.1,
        method_linear = KrylovJL_BICGSTAB(),
        precon_linear = A -> VoronoiFVM.factorization(A, SparspakFactorization()),
    verbose="nlad"
	)


# ╔═╡ 4fa1c608-19b2-4eaa-8c0a-5881b373807c
@test isapprox(krydel_sol, sol0, atol=1.0e-6)

# ╔═╡ 023ae612-6bef-4b0d-8d04-5c0461efbc18
md"""
##### BICGstab from Krylov.jl with ilu0 preconditioner from ILUZero.jl
"""

# ╔═╡ 6895cdf9-8291-47ce-bd1d-4c5beec594ea
    kryilu0_sol = solve(
        sys0;
        inival = 0.5,
        method_linear = KrylovJL_BICGSTAB(),
        precon_linear = ILUZero.ilu0,
        verbose=true
    )

# ╔═╡ d341f60e-191d-4cf9-9df9-fbe25c84a7da
@test isapprox(kryilu0_sol, sol0, atol=1.0e-6)

# ╔═╡ 3bca80f4-ec88-41e8-829a-278071442d41
md"""
### New verbosity handling
"""

# ╔═╡ c2c3088a-6c6e-492e-957b-c5d20fec4244
md"""
- `verbose` can now be a Bool or a String of flag characters, allowing for control of different output categories. I would love to do this via  logging, but there is still a [long way to go](https://github.com/JuliaLang/julia/issues/33418) IMHO 
- Allocation check is active by default with warnings which can be muted by passing a `verbose` string without 'a'. This is now the only control in this respect. All `check_allocs` methods/kwargs, control via environment variables have been removed.
- Deprecation warnings can be switched off by passing a `verbose` string without 'd'.
- Improve iteration logging etc., allow for logging of linear iterations ('l' flag character)


"""

# ╔═╡ ba6b0e35-75ce-4540-861c-7654fd4dee63
md"""
The following example gives some information in this respect:
"""

# ╔═╡ d1440945-54ac-4c12-9f34-0db1c7dfac11
D=0.1

# ╔═╡ 2ff77259-7a81-4c54-a291-cbc20ee56c5d
function xflux(f, u, edge)
        f[1] = D * (u[1, 1]^2 - u[1, 2]^2)
end

# ╔═╡ b0a845d7-95d7-4212-9620-e1948698c596
xsys=VoronoiFVM.System(0:0.001:1,flux=xflux,species=[1])

# ╔═╡ 2a85a42c-1a79-425f-8887-71e2944cb0f3
solve(xsys,inival=0.1,times=[0,1]);

# ╔═╡ 1370a5dc-c434-423f-a0fa-b25f0f2878f9
md"""
If we find these warnings annoying, we can switch them off:
"""

# ╔═╡ 62f605a0-6d5c-4ce2-a131-9ab7b6188d23
solve(xsys,inival=0.1,times=[0,1],verbose="");

# ╔═╡ 2b58318c-8e21-4f9a-97df-ff4af50c94b2
md"""
Or we get some more logging:
"""

# ╔═╡ 05381762-d8f8-46d2-8eb2-68275458787a
solve(xsys,inival=0.1,times=[0,1],verbose="en");

# ╔═╡ f752f390-c98e-4832-b186-f484ebe5a4cb
md"""
But we can also look for the reasons of the allocations. Here, global values should be declared as constants.
"""

# ╔═╡ 0e8bbe4c-66d8-4545-8196-c7d8d9e30bfd
const D1=0.1

# ╔═╡ c40f1954-4fb7-48b2-ab4c-cdf6459b7383
function xflux1(f, u, edge)
        f[1] = D1* (u[1, 1]^2 - u[1, 2]^2)
end

# ╔═╡ cec5152b-0597-4ea4-8387-b08e1e4ffcde
xsys1=VoronoiFVM.System(0:0.001:1,flux=xflux1,species=[1])

# ╔═╡ 76c6ee48-8061-4ba0-b61e-0e6b68ad6435
solve(xsys1,inival=0.1,times=[0,1]);

# ╔═╡ d4ee4693-8ecd-4916-a722-79f54eb99d42
md"""
## v0.14
"""

# ╔═╡ 8a4f336c-2016-453e-9a9f-beac66533fc0
md"""
### `VoronoiFVM.System` constructor
"""

# ╔═╡ 78e7c000-3a83-446a-b577-3a1809c664d2
md"""
### Implicit creation of physics

The `VoronoiFVM.Physics` struct almost never was used outside of the constructor of `VoronoiFVM.System`. Now it is possible to specify the flux functions directly in the system constructor. By default, it is als possible to set a list of species which are attached to all interior and boundary regions of the grid.
"""

# ╔═╡ c59876bd-0cb1-4157-9ba4-bdbede151a44
grid1 = simplexgrid(0:0.1:1);

# ╔═╡ 90bbf212-c6c8-44f0-8132-4a98f094750e
function multispecies_flux(y, u, edge)
    for i = 1:edge.nspec
        y[i] = u[i, 1] - u[i, 2]
    end
end

# ╔═╡ adff41d1-9398-4a66-9a8e-e03809973fa6
function test_reaction(y, u, node)
    y[1] = u[1]
    y[2] = -u[1]
end

# ╔═╡ 5e6d83ab-65c7-4f33-b0a8-29cd5717b4d6
begin
    system1 = VoronoiFVM.System(
        grid1,
        flux = multispecies_flux,
        reaction = test_reaction,
        species = [1, 2],
    )
    boundary_dirichlet!(system1, species = 1, region = 1, value = 1)
    boundary_dirichlet!(system1, species = 2, region = 2, value = 0)
end;

# ╔═╡ f11c03a3-7899-42fd-a2da-a257715815dc
sol1 = solve(system1);

# ╔═╡ 2754d4c8-bbc1-4283-8156-c660c33cd62d
let
    vis = GridVisualizer(resolution = (500, 300), legend = :rt)
    scalarplot!(vis, grid1, sol1[1, :], color = :red, label = "species1")
    scalarplot!(vis, grid1, sol1[2, :], color = :green, label = "species2", clear = false)
    reveal(vis)
end

# ╔═╡ fbd75cf1-64e4-4f07-b54f-f90626f3f6ba
@test isapprox(sum(sol1), 11.323894375033476, rtol = 1.0e-14)

# ╔═╡ b5f7e133-500d-4a27-8f78-11ff3582599c
md"""
### Boundary conditions as part of physics

This makes the API more consistent and opens an easy possibility to have
space and time dependent boundary conditions. One can specify them either in `breaction` or the synomymous `bcondition`.
"""

# ╔═╡ ec188c81-3374-4eed-9b7e-e22350886df2
function bcond2(y, u, bnode)
    boundary_neumann!(y, u, bnode, species = 1, region = 1, value = sin(bnode.time))
    boundary_dirichlet!(y, u, bnode, species = 2, region = 2, value = 0)
end;

# ╔═╡ c86e8a0f-299f-42ab-96f8-0cd62d50f196
system2 = VoronoiFVM.System(
    grid1,
    flux = multispecies_flux,
    reaction = test_reaction,
    species = [1, 2],
    bcondition = bcond2,
    check_allocs = false,
);

# ╔═╡ b3d936fe-69ab-4013-b787-2f0b5410638a
sol2 = solve(system2, times = (0, 10), Δt_max = 0.01);

# ╔═╡ 17749697-d5d8-4629-a625-e96590a5f0ac
vis2 = GridVisualizer(resolution = (500, 300), limits = (-2, 2), legend = :rt)

# ╔═╡ 0c916da5-2d6e-42df-ac4b-4a062f931ccd
md"""
time: $(@bind t2 Slider(0:0.01:10; default=5,show_value=true))
"""

# ╔═╡ 783618f8-2470-4c7c-afc1-9800586625c1
let
    s = sol2(t2)
    scalarplot!(vis2, grid1, s[1, :], color = :red, label = "species1")
    scalarplot!(
        vis2,
        grid1,
        s[2, :],
        color = :green,
        label = "species2",
        clear = false,
        title = "time=$(t2)",
    )
    reveal(vis2)
end

# ╔═╡ 4cbea340-9c02-4e69-8f5e-62bf45312bdd
@test isapprox(sum(sol2) / length(sol2), 2.4908109508494247, rtol = 1.0e-14)

# ╔═╡ 1c18b5a0-cca6-46a1-bb9f-b3d65b8043c5
md"""
### Implicit creation of grid
"""

# ╔═╡ 47280b56-e5ec-4345-b4a1-7c3c92536b2e
md"""
By passing data for grid creation (one  to thre abstract vectors) instead a grid, a tensor product grid is implicitely created.
This example also demonstrates position dependent boundary values.
"""

# ╔═╡ a71086fa-4ec6-4842-a4e1-6a6b60441fc2
function bcond3(y, u, bnode)
    boundary_dirichlet!(y, u, bnode, region = 4, value = bnode[2])
    boundary_dirichlet!(y, u, bnode, region = 2, value = -bnode[2])
end;

# ╔═╡ a514231a-e465-4f05-ba4c-b20aa968d96f
system3 = VoronoiFVM.System(
    -1:0.1:1,
    -1:0.1:1,
    flux = multispecies_flux,
    bcondition = bcond3,
    species = 1,
);

# ╔═╡ d55f615c-d586-4ef7-adf9-5faf052b75ac
sol3 = solve(system3);


# ╔═╡ c17e5104-4d3a-4d54-81c1-d7253245a8bb
@test isapprox(sum(sol3), 0.0, atol = 1.0e-14)

# ╔═╡ 087ea16d-742e-4398-acf5-37248af1b5b4
md"""
### GridVisualize API extended to System
Instead of a grid, a system can be passed to `gridplot` and `scalarplot`.
"""

# ╔═╡ 34d465a5-7cc5-4348-b9ba-6d9381bb3a87
scalarplot(system3, sol3, resolution = (300, 300), levels = 10, colormap = :hot)

# ╔═╡ ac48f5bd-fd1e-4aa7-a2c9-90f0f427143c
md"""
### Parameters of `solve`
"""

# ╔═╡ 69974c02-57e6-4eb5-acf4-b2d480fbd67d
md"""
The `solve` API has been simplified and made more Julian. All entries of `VoronoiFVM.NewtonControl` can be now passed as keyword arguments to `solve`.

Another new keyword argument is `inival` which allows to pass an initial value which by default is initialized to zero. Therefore we now can write `solve(system)` as we already have seen above.
"""

# ╔═╡ 1e12afcf-cf46-4672-9434-44fa8af95ef7
reaction4(y, u, bnode) = y[1] = -bnode[1]^2 + u[1]^4;

# ╔═╡ 938ef63c-58c4-41a0-b3dd-4eb76987a4d7
bc4(args...) = boundary_dirichlet!(args..., value = 0);

# ╔═╡ fe424654-f070-46a9-850a-738b1d4aca8f
system4 = VoronoiFVM.System(
    -10:0.1:10,
    species = [1],
    reaction = reaction4,
    flux = multispecies_flux,
    bcondition = bc4,
);

# ╔═╡ 37fc8816-5ccd-436e-8335-ebb1218d8a35
sol4 = solve(system4, log = true, damp_initial = 0.001, damp_growth = 3);

# ╔═╡ 6a256a29-f15f-4d82-8e84-7ceacb786715
scalarplot(
    system4,
    sol4,
    resolution = (500, 300),
    xlabel = "x",
    ylabel = "u",
    title = "solution",
)

# ╔═╡ 5c2a3836-dc81-4950-88e5-7f603514b1c0
@test isapprox(sum(sol4), 418.58515700568535, rtol = 1.0e-14)

# ╔═╡ fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
html"""<hr><hr><hr>"""

# ╔═╡ f50c6497-cba3-491a-bedd-5f94f88f76fb
md"""
## Appendix: Tests & Development
"""

# ╔═╡ ad899a81-baab-4433-8b7f-1e5c3b18dae6
md"""
This notebook is also run during the automatic unit tests.
Furthermore, the cell activates a development environment if the notebook is loaded from a checked out VoronoiFVM.jl. Otherwise, Pluto's built-in package manager is used.
"""

# ╔═╡ bdbe6513-70b1-4d97-a79c-71534caad2b7
if developing
    @info "Developing VoronoiFVM at  $(pathof(VoronoiFVM))"
else
    @info "Loaded VoronoiFVM from  $(pathof(VoronoiFVM))"
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
ILUZero = "88f59080-6952-5380-9ea5-54057fb9a43f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PlutoVista = "646e1f28-b900-46d7-9d87-d554eb38a413"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
ExtendableGrids = "~0.9.16"
GridVisualize = "~0.6.3"
ILUZero = "~0.2.0"
LinearSolve = "~1.34.1"
PlutoUI = "~0.7.49"
PlutoVista = "~0.8.16"
Revise = "~3.5.0"
VoronoiFVM = "~0.19.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "ff64551fd0cf7d13f1b70fbefdc9cdb552b0e0ee"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "df23d15b1090a3332a09a7a51da45bd9f0a07f92"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.8"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0310e08cb19f5da31d08341c6120c047598f5b9c"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "SnoopPrecompile", "Static"]
git-tree-sha1 = "dedc16cbdd1d32bead4617d27572f582216ccf23"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.25"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e5f08b5689b1aad068e01751889f2f615c7db36d"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.29"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "3d1a9a01976971063b3930d1aed1d9c4af0817f8"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.7"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "f12dc65aef03d0a49650b20b2fdaf184928fd886"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "93c8ba53d8d26e124a5a8d4ec914c3a16e6a0970"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.3"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4aff5fa660eb95c2e0deb6bcdabe4d9a96bc4667"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.8.18"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b277f5d9737fde7b5dce7f43faf49c80f542806c"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.17.10"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "0c5f81f47bbbcf4aea7b2959135713459170798b"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.5"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "2c144ddb46b552f72d7eafe7cc2f50746e41ea21"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "d61300b9895f129f4bd684b2aff97cf319b6c493"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.11"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "0e5c14c3bb8a61b3d53b2c0620570c332c8d0663"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.2.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.Configurations]]
deps = ["ExproniconLite", "OrderedCollections", "TOML"]
git-tree-sha1 = "62a7c76dbad02fdfdaa53608104edf760938c4ca"
uuid = "5218b696-f38b-4ac9-8b61-a12ec717816d"
version = "0.17.4"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "c5b6685d53f933c11404a3ae9822afe30d522494"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.12.2"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "74911ad88921455c6afcad1eefa12bd7b1724631"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.80"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "988e2db482abeb69efc76ae8b6eba2e93805ee70"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.15"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "d0fa82f39c2a5cdb3ee385ad52bc05c42cb4b9f0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.5"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "e1c40d78de68e9a2be565f0202693a158ec9ad85"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.11"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.ExproniconLite]]
deps = ["Pkg", "TOML"]
git-tree-sha1 = "c2eb763acf6e13e75595e0737a07a0bec0ce2147"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.7.11"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "Test", "WriteVTK"]
git-tree-sha1 = "310b903a560b7b18f63486ff93da1ded9cae1f15"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.9.16"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "ILUZero", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "Sparspak", "SuiteSparse", "Test"]
git-tree-sha1 = "075dfc8c0049b676a6af6ce2d7e9e84f56371808"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "0.9.6"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FastLapackInterface]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "7fbaf9f73cd4c8561702ea9b16acf3f99d913fe4"
uuid = "29a986be-02c6-4525-aec4-84b980013641"
version = "1.2.8"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "Setfield", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "04ed1f0029b6b3af88343e439b995141cb0d0b8d"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.17.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "a69dd6db8a809f78846ff259298678f0d6212180"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.34"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "a5e6e7f12607e90d71b09e6ce2c965e41b337968"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.1"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.FuzzyCompletions]]
deps = ["REPL"]
git-tree-sha1 = "e16dd964b4dfaebcded16b2af32f05e235b354be"
uuid = "fb4132e2-a121-4a70-b8a1-d5b831dcdcc2"
version = "0.5.1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "57f7cde02d7a53c9d1d28443b9f11ac5fbe7ebc9"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.3"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "e07a1b98ed72e3cdd02c6ceaab94b8a606faca40"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.2.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ba2d094a88b6b287bd25cfa86f301e7693ffae2f"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.4"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "GridVisualizeTools", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "15c3d40efa9c581f16fd534b42171ac2f93f2cd0"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "0.6.3"

[[deps.GridVisualizeTools]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "StaticArraysCore"]
git-tree-sha1 = "5964fd3e4080af45bfdbdaff75567759fd0367bd"
uuid = "5573ae12-3b76-41d9-b48c-81d0b6e61cc5"
version = "0.2.1"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "47f0f03eddecd7ad59c42b1dd46d5f42916aff63"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.11"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "734fd90dd2f920a2f1921d5388dcebe805b262dc"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.14"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.ILUZero]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "b007cfc7f9bee9a958992d2301e9c5b63f332a90"
uuid = "88f59080-6952-5380-9ea5-54057fb9a43f"
version = "0.2.0"

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
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "3f91cd3f56ea48d4d2a75c2a65455c5fc74fa347"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "c3244ef42b7d4508c638339df1bdbf4353e144db"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.30"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "72ab280d921e8a013a83e64709f66bc3e854b2ed"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.20"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "764164ed65c30738750965d55652db9c94c59bfe"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.4.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "dd90aacbfb622f898a97c2a4411ac49101ebab8a"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.9.0"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "0a92979c14dfa71adbf892f0cd073e34b7189197"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.13.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "0ad6f0c51ce004dadc24a28a0dfecfb568e52242"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.13"

[[deps.LazilyInitializedFields]]
git-tree-sha1 = "410fe4739a4b092f2ffe36fcb0dcc3ab12648ce1"
uuid = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
version = "1.2.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterfaceCore", "DocStringExtensions", "FastLapackInterface", "GPUArraysCore", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "Preferences", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "SnoopPrecompile", "SparseArrays", "Sparspak", "SuiteSparse", "UnPack"]
git-tree-sha1 = "960da8a80f9882fb52a5a199e944d3b86f0d2b94"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.34.1"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "45b288af6956e67e621c5cbb2d75a261ab58300b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.20"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SnoopPrecompile", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "f127dff5c29692069c8ada8267864d252d417821"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.149"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MsgPack]]
deps = ["Serialization"]
git-tree-sha1 = "a8cbf066b54d793b9a48c5daa5d586cf2b5bd43d"
uuid = "99f44e22-a591-53d1-9472-aa23ef4bd671"
version = "1.1.0"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "aa532179d4a643d4bd9f328589ca01fa20a0d197"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.1.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "8175fc2b118a3755113c8e68084dc1a9e63c61ee"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.Pluto]]
deps = ["Base64", "Configurations", "Dates", "Distributed", "FileWatching", "FuzzyCompletions", "HTTP", "HypertextLiteral", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "MsgPack", "Pkg", "PrecompileSignatures", "REPL", "RegistryInstances", "RelocatableFolders", "SnoopPrecompile", "Sockets", "TOML", "Tables", "URIs", "UUIDs"]
git-tree-sha1 = "0ee5bd226e5b95e2232229f7c4a97309ccd8158b"
uuid = "c3e4b0f8-55cb-11ea-2926-15256bba5781"
version = "0.19.22"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

[[deps.PlutoVista]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "GridVisualizeTools", "HypertextLiteral", "Pluto", "UUIDs"]
git-tree-sha1 = "5af654ba1660641b3b80614a7be7eacae4c49875"
uuid = "646e1f28-b900-46d7-9d87-d554eb38a413"
version = "0.8.16"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "e8e0fabcff4df8686c4267503887202a783d498e"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.2"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "240d7170f5ffdb285f9427b92333c3463bf65bf6"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.1"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff"]
git-tree-sha1 = "758f3283aba57c53960c8e1900b4c724bf24ba74"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.8"

[[deps.PrecompileSignatures]]
git-tree-sha1 = "18ef344185f25ee9d51d80e179f8dad33dc48eb1"
uuid = "91cefc8d-f054-46dc-8f8c-26e11d7c5411"
version = "3.0.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "de191bc385072cc6c7ed3ffdc1caeed3f22c74d4"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.7.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "ZygoteRules"]
git-tree-sha1 = "f311e004143b4dc7c5492a2e9b9a1d945058fd8c"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.36.0"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "SnoopPrecompile", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "315b2c85818eea6ad1b6b84fd4ecb40cd4146665"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.17"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegistryInstances]]
deps = ["LazilyInitializedFields", "Pkg", "TOML", "Tar"]
git-tree-sha1 = "ffd19052caf598b8653b99404058fce14828be51"
uuid = "2792f1a3-b283-48e8-9a74-f99dce5104f3"
version = "0.1.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "fd5dba2f01743555d8435f7c96437b29eae81a17"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "50314d2ef65fce648975a8e80ae6d8409ebbf835"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.5"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "dd4195d308df24f33fb10dde7c22103ba88887fa"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "cda0aece8080e992f6370491b08ef3909d1c04e7"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.38"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "RuntimeGeneratedFunctions", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "c4e3c4bb8544b2f0bc537baa44925b6b39a8c383"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.81.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "4245283bee733122a9cb4545748d64e0c63337c0"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.30.0"

[[deps.Sparspak]]
deps = ["Libdl", "LinearAlgebra", "Logging", "OffsetArrays", "Printf", "SparseArrays", "Test"]
git-tree-sha1 = "f9ca0222e4c0b7b1520298a0880e2f7101938831"
uuid = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
version = "0.3.7"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "c35b107b61e7f34fa3f124026f2a9be97dea9e1c"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.3"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "45190b743cdc6f761da1e079bb15ff103a89069c"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.4.6"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+0"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "6b764c160547240d868be4e961a5037f47ad7379"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.1"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "348ad5af9c916b6e1641c74378fac8bb49236688"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.1"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "eb34c1e20b225c1de5adeeff3a085c9e985df532"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.0.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "c97f60dd4f2331e1a495527f80d242501d2f9865"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.1"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f2fd3f288dfc6f507b0c3a2eb3bac009251e548b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.22"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "31eedbc0b6d07c08a700e26d31298ac27ef330eb"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.19"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "4c59c2df8d2676c4691a39fa70495a6db0c5d290"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.58"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.VoronoiFVM]]
deps = ["BandedMatrices", "CommonSolve", "DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "InteractiveUtils", "JLD2", "LinearAlgebra", "LinearSolve", "Printf", "Random", "RecursiveArrayTools", "SnoopPrecompile", "SparseArrays", "SparseDiffTools", "StaticArrays", "Statistics", "SuiteSparse", "Symbolics", "Test"]
git-tree-sha1 = "39f1ffcca22b8c3357789d0f2dc63803f5299fb5"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.19.0"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams"]
git-tree-sha1 = "f1b5fce1c4849bc49212c5f471284abf11c57eec"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.17.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─4ed0c302-26e4-468a-a40d-0e6406f802d0
# ╟─3e6b4ffa-7b33-4b94-9fd6-75b030d5a115
# ╠═7a104243-d3b9-421a-b494-5607c494b106
# ╠═b285aca3-dee5-4b77-9276-537563e8643b
# ╟─a2f1e6ba-80b2-4902-9c5d-2806a7fb16f6
# ╟─56976316-ca4f-4d5c-b303-026edd8751c2
# ╟─740dd980-ca12-4f51-88b2-408930704952
# ╟─db216ebf-1869-4b0f-8dbd-2a8a244c3e4c
# ╠═9eae0139-379d-47ce-a8d2-932c68504317
# ╠═8a31ed03-9cd8-4532-9fb9-8c060f777693
# ╟─4279ae2e-a948-4358-9037-0c6895ecb809
# ╠═58cf3e44-ee53-42f7-9132-eacfd900dc3a
# ╟─b499a195-f06d-400f-9407-b07c3212c095
# ╠═f660fc6e-0ab5-4fae-918f-39bf0c2153e5
# ╟─13f284e8-4fe9-42ce-873a-c65febc7d7df
# ╟─c53cb54f-6099-4d09-9797-3da1f5428586
# ╟─2cebb072-bef1-4bce-9cbb-6b43b877b28b
# ╟─7b0f2021-a2fd-4bb2-a23a-432f61a38a07
# ╟─182a9a9f-3f2c-4df3-abad-08488bb9fb33
# ╟─381c464e-a0cc-46ea-b0ba-089250c31555
# ╟─ba4b8b92-617a-40b4-b3c9-1367067027fa
# ╠═00634dc4-e9c7-4165-a9ac-f3ffc8007e76
# ╠═d6c88c3e-6b8a-45e9-9344-71f0de3fff51
# ╟─cb2d42c4-8d5d-465c-9376-f568588ab453
# ╠═5908ae1d-b3b5-4681-a0b8-080f052af40f
# ╠═22c3cecc-12d2-4f7d-8a53-894a5ea513f0
# ╟─3107c644-5f43-4322-aad2-04a2ab0d576e
# ╟─7e94ff44-4807-41b0-875d-526390903942
# ╠═b70524fc-b8b1-4dee-b77d-e3f8d6d2837b
# ╠═d21d3236-b3d7-4cf8-ab9d-b0be44c9970b
# ╟─f935061c-71db-4aaf-864b-bb566e68643e
# ╠═4c7f8bbf-40c4-45b9-a62e-99ffaae30af1
# ╠═4fa1c608-19b2-4eaa-8c0a-5881b373807c
# ╟─023ae612-6bef-4b0d-8d04-5c0461efbc18
# ╠═6895cdf9-8291-47ce-bd1d-4c5beec594ea
# ╠═d341f60e-191d-4cf9-9df9-fbe25c84a7da
# ╟─3bca80f4-ec88-41e8-829a-278071442d41
# ╟─c2c3088a-6c6e-492e-957b-c5d20fec4244
# ╟─ba6b0e35-75ce-4540-861c-7654fd4dee63
# ╠═d1440945-54ac-4c12-9f34-0db1c7dfac11
# ╠═2ff77259-7a81-4c54-a291-cbc20ee56c5d
# ╠═b0a845d7-95d7-4212-9620-e1948698c596
# ╠═2a85a42c-1a79-425f-8887-71e2944cb0f3
# ╟─1370a5dc-c434-423f-a0fa-b25f0f2878f9
# ╠═62f605a0-6d5c-4ce2-a131-9ab7b6188d23
# ╟─2b58318c-8e21-4f9a-97df-ff4af50c94b2
# ╠═05381762-d8f8-46d2-8eb2-68275458787a
# ╟─f752f390-c98e-4832-b186-f484ebe5a4cb
# ╠═0e8bbe4c-66d8-4545-8196-c7d8d9e30bfd
# ╠═c40f1954-4fb7-48b2-ab4c-cdf6459b7383
# ╠═cec5152b-0597-4ea4-8387-b08e1e4ffcde
# ╠═76c6ee48-8061-4ba0-b61e-0e6b68ad6435
# ╟─d4ee4693-8ecd-4916-a722-79f54eb99d42
# ╟─8a4f336c-2016-453e-9a9f-beac66533fc0
# ╟─78e7c000-3a83-446a-b577-3a1809c664d2
# ╠═c59876bd-0cb1-4157-9ba4-bdbede151a44
# ╠═90bbf212-c6c8-44f0-8132-4a98f094750e
# ╠═adff41d1-9398-4a66-9a8e-e03809973fa6
# ╠═5e6d83ab-65c7-4f33-b0a8-29cd5717b4d6
# ╠═f11c03a3-7899-42fd-a2da-a257715815dc
# ╟─2754d4c8-bbc1-4283-8156-c660c33cd62d
# ╠═fbd75cf1-64e4-4f07-b54f-f90626f3f6ba
# ╟─b5f7e133-500d-4a27-8f78-11ff3582599c
# ╠═ec188c81-3374-4eed-9b7e-e22350886df2
# ╠═c86e8a0f-299f-42ab-96f8-0cd62d50f196
# ╠═b3d936fe-69ab-4013-b787-2f0b5410638a
# ╟─17749697-d5d8-4629-a625-e96590a5f0ac
# ╟─0c916da5-2d6e-42df-ac4b-4a062f931ccd
# ╟─783618f8-2470-4c7c-afc1-9800586625c1
# ╠═4cbea340-9c02-4e69-8f5e-62bf45312bdd
# ╟─1c18b5a0-cca6-46a1-bb9f-b3d65b8043c5
# ╟─47280b56-e5ec-4345-b4a1-7c3c92536b2e
# ╠═a71086fa-4ec6-4842-a4e1-6a6b60441fc2
# ╠═a514231a-e465-4f05-ba4c-b20aa968d96f
# ╠═d55f615c-d586-4ef7-adf9-5faf052b75ac
# ╠═c17e5104-4d3a-4d54-81c1-d7253245a8bb
# ╟─087ea16d-742e-4398-acf5-37248af1b5b4
# ╠═34d465a5-7cc5-4348-b9ba-6d9381bb3a87
# ╟─ac48f5bd-fd1e-4aa7-a2c9-90f0f427143c
# ╟─69974c02-57e6-4eb5-acf4-b2d480fbd67d
# ╠═1e12afcf-cf46-4672-9434-44fa8af95ef7
# ╠═938ef63c-58c4-41a0-b3dd-4eb76987a4d7
# ╠═fe424654-f070-46a9-850a-738b1d4aca8f
# ╠═37fc8816-5ccd-436e-8335-ebb1218d8a35
# ╟─6a256a29-f15f-4d82-8e84-7ceacb786715
# ╠═5c2a3836-dc81-4950-88e5-7f603514b1c0
# ╟─fc0245fe-1bf2-45a3-aa7c-9cce8d7eef37
# ╟─f50c6497-cba3-491a-bedd-5f94f88f76fb
# ╟─ad899a81-baab-4433-8b7f-1e5c3b18dae6
# ╠═e00d0175-866e-4f0f-8121-49e7bbda6fb6
# ╠═bdbe6513-70b1-4d97-a79c-71534caad2b7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
