### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 4b75cde6-db35-11eb-1b95-bba4e05edb2f
begin
    using Pkg
    Pkg.activate(mktempdir())
    Pkg.add("Revise")
	develop=false
    using Revise
    if haskey(ENV,"JULIA_PLUTO_DEVELOP") 
        Pkg.develop(name="VoronoiFVM")
   	    develop=true
    else
        Pkg.add(url="https://github.com/j-fu/VoronoiFVM.jl.git", rev="discontinouos_quantities")
    end
    Pkg.add(["ExtendableGrids","GridVisualize","PlutoVista","PlutoUI"])
    using VoronoiFVM
    using ExtendableGrids
    using GridVisualize
    using PlutoVista
    using PlutoUI
    GridVisualize.default_plotter!(PlutoVista)
	develop
end

# ╔═╡ 3fa189e4-9e1c-470c-bf26-15b631945d2d
md"""
# A discussion of various interface jump conditions


## Two subdomains, current API
For a simple stationary diffusion equation with an interior interface, we discuss possible interface conditions between two subdomains.

Let ``\Omega=\Omega_1\cup\Omega_2`` where ``\Omega_1=(-1,0)`` and ``\Omega_2=(0,1)``.
Let ``\Gamma_1={-1}``,``\Gamma_2={1}``  and ``\Gamma_3={0}``.


Regard the following problem:


``\begin{aligned}
     -\Delta u_1 &= 0 & \text{in}\quad \Omega_1\\ 
     -\Delta u_2 &= 0 & \text{in}\quad \Omega_2\\ 
\end{aligned}
``

with exterior boundary conditions

``u_1|_{\Gamma_1} = g_1`` and ``u_2|_{\Gamma_2} = g_2`` 



For the interior boundary (interface) conditons we set 


``\nabla u_1|_{\Gamma_3}+f_1(u_1,u_2)=0``


``-\nabla u_2|_{\Gamma_3}+f_2(u_1,u_2)=0``

where ``f_1``, ``f_2`` are discussed later.
"""

# ╔═╡ 327af4a8-74cc-4834-ab19-d3a1d6873982
TableOfContents()

# ╔═╡ d5a0ee0d-959d-476b-b3c5-79b741059992
md"""
### Set up
"""

# ╔═╡ f03ff283-c989-4b1a-b73e-2e616054e3db
md"""
Create a grid with two subdomins and an interface in the center.
"""

# ╔═╡ 670c78c1-d0be-4362-975b-2c944620681f
nref=2

# ╔═╡ 79193d53-9dfa-47b7-aed2-c7eb43769b5f
begin
	hmax=0.2/2.0^nref
    hmin=0.05/2.0^nref
    X1=geomspace(-1.0,0.0, hmax,hmin)
    X2=geomspace(0.0,1.0, hmin,hmax)
    X=glue(X1,X2)
    grid=VoronoiFVM.Grid(X)
	
    bfacemask!(grid, [0.0],[0.0],3)
    ## Material 1 left of 0
    cellmask!(grid, [-1.0],[0.0],1)
    ## Material 2 right of 0
    cellmask!(grid, [0.0],[1.0],2)
end;

# ╔═╡ 4cb07222-587b-4a74-a444-43f5433d5b03
gridplot(grid,legend=:rt,resolution=(600,200))

# ╔═╡ 02ec3c0b-6e68-462b-84df-931370cbdcac
md"""
For later use (plotting) extract the two subgrids from the grid
"""

# ╔═╡ 592429a1-108c-4e9e-8961-497c2c31f319
subgrid1=subgrid(grid,[1]);

# ╔═╡ 1fa5d9b3-0558-4558-a5f1-f34f70c8d9a0
subgrid2=subgrid(grid,[2]);

# ╔═╡ 5b539ad2-46d6-43b0-8b3e-f8a7e1ae0a6d
md"""
Define the diffusion flux
"""

# ╔═╡ 6aabfbe1-de7d-49ba-8144-6d364b21b34f
    function flux!(f,u,edge)
        if edge.region==1
            f[1]=u[1,1]-u[1,2]
        end
        if edge.region==2
            f[2]=u[2,1]-u[2,2]
        end
    end

# ╔═╡ dbd27f10-9fd9-450d-9f78-89ea738d605b
md"""
Specify  the outer boundary values.
"""

# ╔═╡ e5b432b2-875e-49a5-8e78-68f7f47f06c3
g_1=1.0

# ╔═╡ 8ae32292-df3e-4190-83a5-ec5ba529299e
g_2=0.1

# ╔═╡ 677853a7-0c43-4d3c-bc74-799535f95aeb
md"""
Create the system. We pass the interface condition function as a parameter.
"""

# ╔═╡ 139058a9-44a0-43d5-a377-4fc72927fa28
function make_system(breaction)
	physics=VoronoiFVM.Physics(
        flux=flux!,
        breaction=breaction
    )

    ## Create system
    sys=VoronoiFVM.System(grid,physics,unknown_storage=:sparse)

    ##  put potential into both regions
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[2])

    ## Set boundary conditions
    for ispec=1:2
        boundary_dirichlet!(sys,ispec,1,g_1)
        boundary_dirichlet!(sys,ispec,2,g_2)
    end
	sys
end

# ╔═╡ de251b03-dd3a-4a44-9440-b7e654c32dac
md"""
Stationary solution with zero initial value
"""

# ╔═╡ 5b3a49d6-30ce-4d6e-94b4-a76645d7d8ce
function mysolve(sys)
	U=solve(unknowns(sys,inival=0),sys)	
	U1=view(U[1,:],subgrid1)
    U2=view(U[2,:],subgrid2)	
	U1,U2
end

# ╔═╡ 0735c061-68e1-429f-80f1-d8410989a91d
md"""
Plot the results
"""

# ╔═╡ 467dc381-3b3d-4de7-a7f9-bfc51300832b
function plot(U1,U2;title="")
	vis=GridVisualizer(resolution=(600,300))   
    scalarplot!(vis,subgrid1,U1,clear=false,show=false, color=:green, label="u1")
    scalarplot!(vis,subgrid2,U2,clear=false,show=true, color=:blue, label="u2",legend=:rt,title=title,flimits=(-0.5,1.5))
end

# ╔═╡ fa4fcc0d-1d3a-45a2-8857-50536bbe39cc
md"""
### No interface reaction

This means we set ``f_1(u_1,u_2)=0`` and ``f_2(u_1,u_2)=0``. 
"""

# ╔═╡ 8f210696-fcf4-47bc-a5a2-c561ad7efcbd
function noreaction(f,u,node) end

# ╔═╡ 57e8515e-3be1-4478-af98-430501438ee7
system1=make_system(noreaction);

# ╔═╡ 56136cd1-0c01-449d-9297-68924ac99ee7
plot(mysolve(system1)...)

# ╔═╡ aad305a9-aac6-4aff-9f8e-08d6a2f756c8
md"""
### Mass action law reaction ``u_1 \leftrightharpoons u_2``

This is a rather general ansatz where we assume a backward-forward reaction at the interface with reaction constants ``k_1`` and ``k_2``, respectively.

According to the mass action law, this tranlates to a reaction rate

``r(u_1,u_2)=k_1u_1 - k_2u_2``

and correspondingly

``f_1(u_1,u_2)=r``

``f_2(u_1,u_2)=-r`` 

Note, that ``f_i`` is monotonically increasing in ``u_i`` and monotonically decreasing in the respective other argument, leading to an M-Property of the overall discretization matrix.


Note that the "no reaction" case is just a special case where ``k_1,k_2=0``.
"""

# ╔═╡ d027ff24-3ad1-4528-b5cf-10814caf30db
k1=1; k2=0.1;

# ╔═╡ 1328b4bf-2d64-4b02-a910-1995da8be28b
function mal_reaction(f,u,node)
        if node.region==3
            react=k1*u[1]-k2*u[2]
            f[1]= react
            f[2]= -react
		
        end
end

# ╔═╡ 610a0761-1c23-415d-a187-f7d93a1b7637
system2=make_system(mal_reaction)

# ╔═╡ 87edce1f-df6d-4cd8-bce5-24fb666cd6b5
begin
	k1,k2
	U1,U2=mysolve(system2)
    plot(U1,U2;title="k1=$(k1), k2=$(k2)")
end

# ╔═╡ b82fc6b2-eee1-4a91-a115-61b86621f686
md"""
### Penalty enforcing continuity


Setting ``k_1,k_2`` to a large number leads to another special case of the above reaction - similar to the penalty method to implement the Dirichlet boundary conditions, this lets the reaction equation dominate, which in this case forces
``u_1-u_2=0`` at the interface, and thus continuity.
"""

# ╔═╡ 9eaea813-2628-47d0-9d36-54c367689142
function penalty_reaction(f,u,node)
        if node.region==3
            react=1.0e10*(u[1]-u[2])
            f[1]= react
            f[2]= -react
        end
end

# ╔═╡ 817738c0-f1a3-4779-9075-7ea051a81e73
system3=make_system(penalty_reaction);

# ╔═╡ 80019ef1-bf41-4a55-9262-613a2d20be1f
plot(mysolve(system3)...)

# ╔═╡ 31e00855-8906-4be5-8e69-e2d8d9539e04
md"""
### Interface recombination

Here, we implement an annihilation reaction ``u_1 + u_2 \to \emptyset``
Acoording to the mass action law, this is implemented via

``r(u_1,u_2)=k_r u_1 u_2``

``f_1(u_1,u_2)=r``

``f_2(u_1,u_2)=r``



"""

# ╔═╡ f2490f99-04ca-4f42-af2a-53adae51ca68
k_r=1.0

# ╔═╡ 39a0db1b-3a4e-4108-b43f-d4e578c92608
function recombination(f,u,node)
        if node.region==3
            react=k_r*(u[1]*u[2])
            f[1]= react
            f[2]= react
        end
end;

# ╔═╡ 644149fb-2264-42bd-92c9-193ab07c08f6
system4=make_system(recombination);

# ╔═╡ b479402f-ef00-4425-8b0f-45f2dae74d80
plot(mysolve(system4)...)

# ╔═╡ ed068b51-92af-48d5-9230-debc178ec827
md"""
### Thin  conductive interface layer

Let us assume that the interface is of thickness $d$ which is however small with respect to ``\Omega`` that we want to derive an interface condition from the assumption of an exact continuous solution within the interface.

So let ``\Omega_I=(x_l,x_r)`` be  the interface region where
we have ``-\Delta u_I=0`` with values ``u_l``, ``u_r`` at the boundaries. 

Then we have for the flux  in the interface region, ``q_I=\nabla u = \frac1{d}(u_r - u_l)``

Continuity of fluxes then gives ``f_1=q_I`` and ``f_2=-q_I``.

Continuity of ``u`` gives ``u_{1,I}=u_l, u_{2,I}=u_r``
This gives

``r=q_I=\frac{1}{d}(u_1-u_{2})``

``f_1(u_1,v_1)=r``

``f_2(u_1,v_1)=-r``

and therefore another special case of the mass action law condition.
"""

# ╔═╡ a2d919a5-a395-40fb-8f93-db742f8a77c2
d=0.5

# ╔═╡ 58d8831b-ad66-4f77-a33a-933c15c46a52
function thinlayer(f,u,node)
        if node.region==3
            react=(u[1]-u[2])/d
            f[1]= react
            f[2]= -react
        end
end

# ╔═╡ 8c0b4ab5-09da-4d8f-b001-5e15f823423c
system5=make_system(thinlayer);

# ╔═╡ d3d99b9c-ad18-4a04-bb3e-f17dd542f9f3
plot(mysolve(system5)...)

# ╔═╡ eb9abf2e-372c-4f79-afbe-772b90eff9ad
md"""
## An API layer for discontinuous species handling

From this discussion it seems that discontinuous interface conditions can be formulated in a rather general way via linear or nonlinear robin boundary conditions for each of the adjacent discontinuous species. Technically, it is necessary to be able to access the adjacent bulk data.
"""

# ╔═╡ 9f3ae7b5-51b3-48bc-b4db-b7236ba30682
md"""
Here, we propose an API layer on top  of the species handling of VoronoiFVM.
For a start, we call these "Meta species" "quantities". There may be a different name, however.
"""

# ╔═╡ 2da8a5b1-b168-41d9-baa8-d65a4ef5c4c0
md"""
We define a grid with N subregions
"""

# ╔═╡ d44407de-8c9c-42fa-b1a2-ae02b826eccc
N=5

# ╔═╡ ae268316-c058-4db8-9b71-57b0d9425274
begin
	XX=collect(0:0.1:1)
	xcoord=XX
	for i=1:N-1
		xcoord=glue(xcoord,XX.+i)
	end
	grid2=simplexgrid(xcoord)
	for i=1:N
		cellmask!(grid2,[i-1],[i],i)
	end	
	for i=1:N-1
		bfacemask!(grid2,[i],[i],i+2)
	end
end

# ╔═╡ b53b9d28-4c25-4fb8-a3e4-599b0e385121
gridplot(grid2,legend=:lt,resolution=(600,200))

# ╔═╡ e7ce7fd4-cfa4-4cc6-84a2-7e20ed2f4e5c
md"""
To work with quantities, we first introduce a new constructor call without the "physics" parameter:
"""

# ╔═╡ 29f36902-e355-4b02-b7b0-c4db12c47d33
system6=VoronoiFVM.System(grid2,unknown_storage=:sparse)

# ╔═╡ 673e9320-ea30-4831-ad85-ba7936293ee2
md"""
First, we introduce a continuous quantity which we name "cspec". Note that the "species number" can be assigned automatically if not given explicitely.
"""

# ╔═╡ f35f419a-94dd-4051-a533-4b1ec9a4c9ec
cspec=ContinuousQuantity(system6,1:N;ispec=1)

# ╔═╡ 9661e4fc-55e1-4c2c-a3ad-515cdac3b514
md"""
A discontinuous quantity can be introduced as well. by default, each reagion gets a new species number. This can be overwritten by the user.
"""

# ╔═╡ 90298676-fda7-4168-8a40-7ff53e7c761b
dspec=DiscontinuousQuantity(system6,1:N; regionspec=[2+i%2 for i=1:N])

# ╔═╡ cebabf33-e769-47bd-b6f1-ddf525fea895
md"""
For both quantities, we define simple diffusion fluxes:
"""

# ╔═╡ 719f206a-5b9f-4d78-8778-1d89edb2bc4d
	function flux2(f,u,edge)
		f[dspec]=u[dspec,1]-u[dspec,2]
		f[cspec]=u[cspec,1]-u[cspec,2]
	end

# ╔═╡ 1d7f442f-c057-4379-8a40-c6ce3646ad5c
md"""
Define a thin layer inteface condition for `dspec` and an interface source for `cspec`.
"""

# ╔═╡ da41b22e-114d-4eee-81d0-73e6f3b45242
md"""
Add physics to the system, set dirichlet bc at both ends, and extract subgrids
for plotting (until there will be a plotting API for this...)
"""

# ╔═╡ 0c3ca093-e011-409c-be58-e617e02f5d4c
function plot2(U,subgrids,system6)
	dvws=VoronoiFVM.views(U,dspec,subgrids,system6)
	cvws=VoronoiFVM.views(U,cspec,subgrids,system6)
	vis=GridVisualizer(resolution=(600,300))
	for i=1:length(dvws)
		scalarplot!(vis,subgrids[i],dvws[i],flimits=(-0.5,1.5),clear=false,color=:red)
		scalarplot!(vis,subgrids[i],cvws[i],flimits=(-0.5,1.5),clear=false,color=:green)
 	end
	reveal(vis)
end

# ╔═╡ b8cd6ad1-d323-4888-bbd1-5deba5a5870d
d1=1

# ╔═╡ 441a39a0-a7de-47db-8539-12dee30b8312
q1=0.2

# ╔═╡ d6e1c6c7-060d-4c2f-8054-d8f33f54bd55
function breaction2(f,u,node)
	       if node.region>2
            react=(u[dspec,1]-u[dspec,2])/d1
            f[dspec,1]= react
            f[dspec,2]= -react
	
		    f[cspec]=-q1
        end
end

# ╔═╡ 59c83a22-a4cc-4b51-a1cc-5eb39588eacd
begin
	physics!(system6,VoronoiFVM.Physics(
        flux=flux2,
        breaction=breaction2
    ))
	

    ## Set boundary conditions
    boundary_dirichlet!(system6,dspec,1,g_1)
    boundary_dirichlet!(system6,dspec,2,g_2)
    boundary_dirichlet!(system6,cspec,1,0)
    boundary_dirichlet!(system6,cspec,2,0)
	subgrids=VoronoiFVM.subgrids(dspec,system6)
end;

# ╔═╡ 5d8aca85-f12d-4a6e-84e4-781d45b65742
plot2(solve(unknowns(system6,inival=0),system6),subgrids,system6)

# ╔═╡ 6203d95b-13d2-48a6-a69f-39e33e2edcdb
md"""

## Open problems
- Testfunctions: this will be a another problem to be handeled
- Alternative is "glueing together" systems
- Forwarddiff sees the full species matrix. A slight remedy is the possibility to assign species numbers by the user.
"""

# ╔═╡ Cell order:
# ╟─3fa189e4-9e1c-470c-bf26-15b631945d2d
# ╟─327af4a8-74cc-4834-ab19-d3a1d6873982
# ╟─d5a0ee0d-959d-476b-b3c5-79b741059992
# ╠═4b75cde6-db35-11eb-1b95-bba4e05edb2f
# ╟─f03ff283-c989-4b1a-b73e-2e616054e3db
# ╠═670c78c1-d0be-4362-975b-2c944620681f
# ╠═79193d53-9dfa-47b7-aed2-c7eb43769b5f
# ╠═4cb07222-587b-4a74-a444-43f5433d5b03
# ╟─02ec3c0b-6e68-462b-84df-931370cbdcac
# ╠═592429a1-108c-4e9e-8961-497c2c31f319
# ╠═1fa5d9b3-0558-4558-a5f1-f34f70c8d9a0
# ╟─5b539ad2-46d6-43b0-8b3e-f8a7e1ae0a6d
# ╠═6aabfbe1-de7d-49ba-8144-6d364b21b34f
# ╟─dbd27f10-9fd9-450d-9f78-89ea738d605b
# ╠═e5b432b2-875e-49a5-8e78-68f7f47f06c3
# ╠═8ae32292-df3e-4190-83a5-ec5ba529299e
# ╟─677853a7-0c43-4d3c-bc74-799535f95aeb
# ╠═139058a9-44a0-43d5-a377-4fc72927fa28
# ╟─de251b03-dd3a-4a44-9440-b7e654c32dac
# ╠═5b3a49d6-30ce-4d6e-94b4-a76645d7d8ce
# ╟─0735c061-68e1-429f-80f1-d8410989a91d
# ╠═467dc381-3b3d-4de7-a7f9-bfc51300832b
# ╟─fa4fcc0d-1d3a-45a2-8857-50536bbe39cc
# ╠═8f210696-fcf4-47bc-a5a2-c561ad7efcbd
# ╠═57e8515e-3be1-4478-af98-430501438ee7
# ╠═56136cd1-0c01-449d-9297-68924ac99ee7
# ╟─aad305a9-aac6-4aff-9f8e-08d6a2f756c8
# ╠═1328b4bf-2d64-4b02-a910-1995da8be28b
# ╠═610a0761-1c23-415d-a187-f7d93a1b7637
# ╠═d027ff24-3ad1-4528-b5cf-10814caf30db
# ╟─87edce1f-df6d-4cd8-bce5-24fb666cd6b5
# ╟─b82fc6b2-eee1-4a91-a115-61b86621f686
# ╠═9eaea813-2628-47d0-9d36-54c367689142
# ╠═817738c0-f1a3-4779-9075-7ea051a81e73
# ╠═80019ef1-bf41-4a55-9262-613a2d20be1f
# ╟─31e00855-8906-4be5-8e69-e2d8d9539e04
# ╠═39a0db1b-3a4e-4108-b43f-d4e578c92608
# ╠═644149fb-2264-42bd-92c9-193ab07c08f6
# ╠═f2490f99-04ca-4f42-af2a-53adae51ca68
# ╠═b479402f-ef00-4425-8b0f-45f2dae74d80
# ╟─ed068b51-92af-48d5-9230-debc178ec827
# ╠═58d8831b-ad66-4f77-a33a-933c15c46a52
# ╠═8c0b4ab5-09da-4d8f-b001-5e15f823423c
# ╠═d3d99b9c-ad18-4a04-bb3e-f17dd542f9f3
# ╠═a2d919a5-a395-40fb-8f93-db742f8a77c2
# ╟─eb9abf2e-372c-4f79-afbe-772b90eff9ad
# ╟─9f3ae7b5-51b3-48bc-b4db-b7236ba30682
# ╟─2da8a5b1-b168-41d9-baa8-d65a4ef5c4c0
# ╠═d44407de-8c9c-42fa-b1a2-ae02b826eccc
# ╠═ae268316-c058-4db8-9b71-57b0d9425274
# ╠═b53b9d28-4c25-4fb8-a3e4-599b0e385121
# ╟─e7ce7fd4-cfa4-4cc6-84a2-7e20ed2f4e5c
# ╠═29f36902-e355-4b02-b7b0-c4db12c47d33
# ╟─673e9320-ea30-4831-ad85-ba7936293ee2
# ╠═f35f419a-94dd-4051-a533-4b1ec9a4c9ec
# ╟─9661e4fc-55e1-4c2c-a3ad-515cdac3b514
# ╠═90298676-fda7-4168-8a40-7ff53e7c761b
# ╟─cebabf33-e769-47bd-b6f1-ddf525fea895
# ╠═719f206a-5b9f-4d78-8778-1d89edb2bc4d
# ╟─1d7f442f-c057-4379-8a40-c6ce3646ad5c
# ╠═d6e1c6c7-060d-4c2f-8054-d8f33f54bd55
# ╟─da41b22e-114d-4eee-81d0-73e6f3b45242
# ╠═59c83a22-a4cc-4b51-a1cc-5eb39588eacd
# ╠═0c3ca093-e011-409c-be58-e617e02f5d4c
# ╠═5d8aca85-f12d-4a6e-84e4-781d45b65742
# ╠═b8cd6ad1-d323-4888-bbd1-5deba5a5870d
# ╠═441a39a0-a7de-47db-8539-12dee30b8312
# ╟─6203d95b-13d2-48a6-a69f-39e33e2edcdb
