### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 60941eaa-1aea-11eb-1277-97b991548781
begin 
    ENV["LC_NUMERIC"]="C" # prevent pyplot from messing up string2float conversion
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
TableOfContents(title="",depth=4,aside=false)

# ╔═╡ b61dce01-f84f-4e46-be76-8561bf0dddaa
md"""
# Handling of interfaces in VoronoiFVM

Let ``\Omega=(-0.5,1.5)\times (0,2)`` be subdivided into two
subregions, i.e. ``\Omega = \Omega_1\cap \Omega_2``. Let  ``\Gamma_{12}`` be the interface between them. Further, let ``Γ_1=\{-0.5\}\times (0,2)``,
        ``Γ_2=\{1.5\}\times (0,2)``, and ``Γ_N=\partial\Omega \setminus (Γ_1\cup Γ_2)``.
We assume that ``\Omega`` is triangulated with a Delaunay triangulation which conforms to interior
and exterior boundaries. The later property ensures, that subdomain boundaries are aligned
with triangle edges, and  that for $i=1,2$ the circumcenters
of triangles from ``\Omega_i`` are situated  in ``\bar\Omega_i``.
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
-\nabla\cdot D\nabla u &=0, & \text{in}\quad Ω,\\
                     D\nabla u\cdot \vec n + αu &=g_1,& \text{on}\quad Γ_1, \\
                     D\nabla u\cdot \vec n + αu &=g_2,& \text{on}\quad Γ_2,\\
                     D\nabla u\cdot \vec n &=0,& \text{on}\quad Γ_N,
\end{aligned}
```

with ``D|_{\Omega_i}= D_i>0`` and ``α>0``. In fact, this is a sloppy version of a more precise formulation.

Let ``u_i`` be solely defined on ``\Omega_i, i = 1, 2`` and satisfy 
```math
\begin{aligned}
-\nabla\cdot D_i\nabla u_i &=0, & \text{in}\quad Ω_i,\\
                     D_i\nabla u_i\cdot \vec n_i+ α u_i &= g_i,& \text{on}\quad Γ_i, \\
                     D_i\nabla u_i\cdot \vec n_i &=0,& \text{on}\quad Γ_N\cap \partial Ω_i.
\end{aligned}
```
This means, two seperate problems on ``\Omega_1, \Omega_2`` are considered. With the following interface conditions (in the sense of traces of functions) ``u_1, u_2`` are linked
```math
\begin{aligned}
    u_1 &= u_2,  & \text{on} \quad  Γ_{12},\\
 D_1\nabla u_1\cdot \vec n_1 &=-D_2\nabla u_2\cdot \vec n_2, & \text{on} \quad Γ_{12}.
\end{aligned}
```
From the later conditions we can directly follow the continuity of a solution ``u`` with ``u_{|\Omega_1} = u_1, u_{|\Omega_2} = u_2`` and the continuity of respective fluxes.

**Weak Formulation**

Let ``V=H^1(\Omega)`` and denote by ``v \in V`` a respective test function. Integrating over the domain, splitting the problem onto the subdomains, multiplying by the testfunction and using Greens first identity yields
```math
\begin{aligned}
  0 =&   \sum_{i=1}^2 \Bigl( \int_{Ω_i} D_i\nabla u \nabla v dx - \int_{Γ_i} (D_i\nabla u_i )\cdot \vec n_i v dγ  - \int_{Γ_N} (D_i\nabla u_i )\cdot \vec n_i v dγ\Bigr) \\
~
&  - \int_{Γ_{12}} \Bigl((D_1\nabla u_1)\cdot \vec n_1 + (D_2\nabla u_2 )\cdot \vec n_2 \Bigr)v dγ \quad ∀ v ∈ V
\end{aligned}
```
The integral over ``\Gamma_N`` vanishes due to the homogenous Neumann boundary conditions, the integral over ``Γ_{12}`` likewise vanishes due to the interface conditions. Inserting the interface conditions on ``\Gamma_i, i = 1, 2`` brings us the following weak formulation:


Find ``u∈V`` such that
```math
    \sum_{i=1}^2 \int_{Ω_i} D_i\nabla u \nabla v dx + \sum_{i=1}^2 \int_{Γ_i} \alpha uv dγ =  \sum_{i=1}^2 \int_{Γ_i}  g_iv dγ \quad ∀ v ∈ V.
```


A P1 finite element discretization is straightforward if the interface is aligned with triangle edges.
"""

# ╔═╡ d95afcb5-68f9-4d2e-9111-06fc38fe949a
md"""
### Finite volume approach
We introduce the following notation: Let $\omega_K$ be a respective control volume. We denote by ``\sigma_{KL} = \partial \omega_K \cap \partial\omega_L`` the face between two adjacent control volumes. To each control volume we assign a collocation point ``x_K \in \omega_K``. The set of all neighbouring control volumes of ``\omega_K`` is given by ``\mathcal{N}(\omega_K)``. The Euclidean distance between two collocation points is given by ``h_{KL} = \Vert x_L - x_K \Vert_2``. Note that ``h_{KL}`` and ``\sigma_{KL}`` are orthogonal to each other.

In the following we consider the integral over a control volume ``\omega_K``
```math
0=-\int_{\omega_{K}} \nabla D\nabla u dx = -\int_{\partial\omega_{K}} (D\nabla u) \cdot \vec n_K d\gamma
```
and the respective discrete version of it.


#### Interior control volume, e.g. ``ω_{13}``
We can further rewrite the previous integral equation as a sum over shared faces, i.e.

```math
0= \sum_{\omega_L \in \mathcal{N}(\omega_K)} \int_{\sigma_{KL}} - D\nabla u \cdot \vec n_{KL} d\gamma \approx  \sum_{\omega_L \in \mathcal{N}(\omega_K)} \frac{|\sigma_{KL}|}{h_{KL}} D(u_K-u_L).
```

#### Boundary  control volume, e.g. ``ω_{8}``
We denote by ``\gamma_{K,j}=\partial\omega_K\cap Γ_j, j = 1, 2`` the faces of ``\omega_K`` which are located at the outer boundary. Thus, the reformulation of the integrated flux over a boundary control volume ``\omega_K`` for a fixed ``j`` reads

```math
\begin{aligned}
0=&\sum_{\substack{\omega_L \in \mathcal{N}(\omega_K)\\ \partial\omega_L \cap \Gamma_j = \emptyset}} \int_{\sigma_{KL}} -  D\nabla u \cdot \vec n_{KL} d\gamma + \int_{\gamma_{K, j}} - D\nabla u \cdot \vec n_{j}d\gamma  \\
  \approx&  \sum_{\substack{\omega_L \in \mathcal{N}(\omega_K)\\ \partial\omega_L \cap \Gamma_j = \emptyset}} \frac{|\sigma_{KL}|}{h_{KL}} D(u_K-u_L) +  |\gamma_{K, j}| (\alpha u_K-g_j)
\end{aligned}
```

#### Interface  control volume, e.g. ``ω_{7}``
This is the case we are interested to and where the magic happens. The control colume ``\omega_K`` is part of both domains ``\Omega_1, \Omega_2``. We can define ``\omega_{K,j} = \omega_K \cap \Omega_j, j = 1,2 `` as the part of ``\omega_K`` which is solely defined on one subdomain. Thus, we likewise introduce ``\sigma_{KL}^j = \partial \omega_{K,j} \cap \partial \omega_L`` as the respective face between a neighboring control volume and the sub-control volume. Lastly, we use the following notation for the face of ``\omega_{K,j}`` intersecting with the inner boundary ``Γ_{12}``, which is likewise the part of the interface intersecting ``\omega_K``: 
``\gamma_{K,12}=\omega_K\cap Γ_{12}``.


Integrating our PDE over this control volume ``\omega_K`` yields
```math
\begin{aligned}
    0  =& -\int_{\omega_K} -\nabla\cdot D\nabla u dx\\
    =& -\int_{\omega_{K,1}} -\nabla\cdot D_1\nabla u dx -\int_{\omega_{K,2}} -\nabla\cdot D_2\nabla udx\\
       =&  \sum_{\omega_L \in \mathcal{N}(\omega_{K,1})} \int_{\sigma_{KL}^1} - D_1\nabla u \cdot \vec n_{KL} d\gamma + \sum_{\omega_L \in \mathcal{N}(\omega_{K,2})} \int_{\sigma_{KL}^2} -D_2\nabla u \cdot \vec n_{KL} d\gamma \\
        & - \int_{\gamma_{K,12}}  D_1\nabla u \cdot \vec n_1 + D_2\nabla u \cdot \vec n_2 d\gamma \\
    \approx & \sum_{\omega_L \in \mathcal{N}(\omega_{K,1})} \frac{|\sigma_{KL}^1|}{h_{KL}} D_1(u_K-u_L) + \sum_{\omega_L \in \mathcal{N}(\omega_{K,2})} \frac{|\sigma_{KL}^2|}{h_{KL}} D_2(u_K-u_L)
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

# ╔═╡ b16061ae-8a08-4649-af8e-58ddf054996d
md"""
# Handling of a boundary flux

In this section, we will study the effect of the boundary flux. We consider the following test problem 
```math
- \nabla D \nabla u + k u = c
```
on two domains which shall in the limiting case coincide.
"""

# ╔═╡ 3d0eebe4-5eaa-4339-b6f6-bb5acd573033
md"""
The necessary parameters for the problem are defined here, likewise as the respective physics functions.
"""

# ╔═╡ d7910827-238a-4b71-ad44-7ce43ff1c0db
begin
	dA_bulk       = 2.0
	dB_bulk       = 2.0
	dbB_interface = 5.0
	k             = 6.0
	c             = 1.0
end

# ╔═╡ 7b51cf5a-fa12-4756-a016-ed00cf90b1ba
md"""
## Domain A: Introducing a small interface layer

We introduce the height of this thin layer by following variable as a variable for the mesh refinement
"""



# ╔═╡ 3b429e0c-c07a-40df-ab0b-1f71d515ece3
n          = 3

# ╔═╡ c770801b-d8f8-4095-98ff-e546f661e141
thin_layer = 0.1 # just set 0.01, then you see agreement

# ╔═╡ 7d902e19-4ced-4685-ac5b-129ce101c9fe
begin
	function fluxA!(f, u, edge)
        if edge.region == 3
            f[1] =  (dA_bulk + dbB_interface) * (u[1,1] - u[1,2])
        else
            f[1] = dA_bulk *  (u[1,1] - u[1,2])
        end
    end
 
    function fluxB!(f, u, edge)
        f[1] =  dB_bulk *  (u[1,1] - u[1,2])
    end
         
    function reaction!(f, u, node)      
        f[1] = k *  u[1]
    end
         
    function source!(f, node)
        f[1] = c
    end
 
    function bfluxB!(f, u, bedge)
        if bedge.region == 3
            f[1] =  thin_layer *  dbB_interface *  (u[1,1] - u[1,2])
        end
    end

    function bfluxA!(f, u, bedge)
        f[1] = 0.0
    end
end

# ╔═╡ fd59fa72-bc40-41ef-bdac-485b171da03e
begin
    h_pdopingA         = 1.0
	h_intrinsic1A      = 1.0/2
	h_thinlayerA       = thin_layer
	h_intrinsic2A      = 1.0/2 - thin_layer
    h_ndopingA         = 1.0
    h_totalA           = h_pdopingA + h_intrinsic1A + h_thinlayerA + h_intrinsic2A + h_ndopingA
    lengthA            = 1.0 

    regionAcceptorA    = 1
    regionIntrinsic1A  = 2
    regionThinLayerA   = 3
    regionIntrinsic2A  = 4
    regionDonorA       = 5
    regionsA           = [regionAcceptorA, regionIntrinsic1A, regionThinLayerA, regionIntrinsic2A, regionDonorA]

    # boundary region numbers
    bregionAcceptorA   = 1
    bregionDonorA      = 2
    
    coord_pA       = collect(range(0.0,                                    stop = h_pdopingA,                                          length = n))
    coord_i_1A     = collect(range(h_pdopingA,                             stop = h_pdopingA+h_intrinsic1A,                            length = n))
    coord_i_thinLA = collect(range(h_pdopingA+h_intrinsic1A,               stop = h_pdopingA+h_intrinsic1A+h_thinlayerA,               length = n))
    coord_i_2A     = collect(range(h_pdopingA+h_intrinsic1A+h_thinlayerA,  stop = h_pdopingA+h_intrinsic1A+h_thinlayerA+h_intrinsic2A, length = n))
    coord_nA       = collect(range(h_totalA-h_ndopingA,                    stop = h_totalA,                                            length = n))

    coordA         = glue(coord_pA, coord_i_1A, tol=1.0e-30)
    coordA         = glue(coordA,   coord_i_thinLA, tol=1.0e-30)
    coordA         = glue(coordA,   coord_i_2A, tol=1.0e-30) 
    coordA         = glue(coordA,   coord_nA, tol=1.0e-30)
    coord_lengthA  = collect(range(0.0, stop = lengthA, length=3*n))

    gridA          = simplexgrid(coordA, coord_lengthA)
    numberOfNodesA = size(gridA[Coordinates])[2]
 
    println("number of nodes in grid A are: ", numberOfNodesA)

    # specify inner regions
    cellmask!(gridA,[0.0, 0.0],                                  [h_pdopingA, lengthA],                                         regionAcceptorA,   tol=1.0e-30)
    cellmask!(gridA,[h_pdopingA, 0.0],                           [h_pdopingA+h_intrinsic1A, lengthA],                           regionIntrinsic1A, tol=1.0e-30)
    cellmask!(gridA,[h_pdopingA+h_intrinsic1A, 0.0],             [h_pdopingA+h_intrinsic1A+h_thinlayerA, lengthA],              regionThinLayerA,  tol=1.0e-30)
    cellmask!(gridA,[h_pdopingA+h_intrinsic1A+h_thinlayerA, 0.0],[h_pdopingA+h_intrinsic1A+h_thinlayerA+h_intrinsic2A, lengthA],regionIntrinsic2A, tol=1.0e-30) 
    cellmask!(gridA,[h_totalA-h_ndopingA, 0.0],                  [h_totalA, lengthA],                                           regionDonorA,      tol=1.0e-30)
 
    # specifiy outer regions
    # metal interfaces
    bfacemask!(gridA, [0.0, lengthA],  [h_pdopingA, lengthA], bregionAcceptorA, tol=1.0e-30) # BregionNumber = 1
    bfacemask!(gridA, [h_totalA, 0.0], [h_totalA, lengthA],   bregionDonorA,    tol=1.0e-30) # BregionNumber = 2

    bfacemask!(gridA, [h_pdopingA+1.0/2, 0.0], [h_pdopingA+1.0/2, lengthA], 3,  tol=1.0e-30)

    # no flux interfaces [xmin, ymin], [xmax, ymax]
    # @Jürgen: If not defining this, then one of these boundaries is set to 1??
    bfacemask!(gridA, [0.0, 0.0],            [0.0, lengthA],      4,  tol=1.0e-30) 
    bfacemask!(gridA, [h_pdopingA, lengthA], [h_totalA, lengthA], 4,  tol=1.0e-30) 
    bfacemask!(gridA, [0.0, 0.0],            [h_totalA, 0.0],     4,  tol=1.0e-30) 

    gridplot(gridA, Plotter=PyPlot, legend=:rt)
end

# ╔═╡ 0875ce0d-a22b-4e08-9b58-57673c986a67
md"""
## Domain B: Defining an interior interface

"""

# ╔═╡ a9d85817-df56-48ea-8aee-f700dc6310c8
begin
	h_pdopingB       = 1.0
    h_intrinsicB     = 1.0
    h_ndopingB       = 1.0
    h_totalB         = h_pdopingB + h_intrinsicB + h_ndopingB
    lengthB          = 1.0

    # region numbers
    regionAcceptorB  = 1
    regionIntrinsicB = 2
    regionDonorB     = 3
    regionsB         = [regionAcceptorB, regionIntrinsicB, regionDonorB]

    # boundary region numbers
    bregionAcceptorB = 1
    bregionDonorB    = 2
    binnerInterfaceB = 3

    coord_pB         = collect(range(0.0,                       stop = h_pdopingB,                length = n))
    coord_i1B        = collect(range(h_pdopingB,                stop = h_pdopingB+h_intrinsicB/2, length = n))
    coord_i2B        = collect(range(h_pdopingB+h_intrinsicB/2, stop = h_pdopingB+h_intrinsicB,   length = n))
    coord_nB         = collect(range(h_pdopingB+h_intrinsicB,   stop = h_totalB,                  length = n))

    coordB           = glue(coord_pB,  coord_i1B)
    coordB           = glue(coordB, coord_i2B)
    coordB           = glue(coordB,    coord_nB)
    coord_lengthB    = collect(range(0.0, stop = lengthB, length=4*n))

    gridB            = simplexgrid(coordB, coord_lengthB)
    numberOfNodesB   = size(gridB[Coordinates])[2]
 
    println("number of nodes in grid B are: ", numberOfNodesB)

    # specify inner regions
    cellmask!(gridB, [0.0, 0.0],                      [h_pdopingB, lengthB],              regionAcceptorB,  tol=1.0e-18)
    cellmask!(gridB, [h_pdopingB, 0.0],               [h_pdopingB+h_intrinsicB, lengthB], regionIntrinsicB, tol=1.0e-18) 
    cellmask!(gridB, [h_pdopingB+h_intrinsicB, 0.0], [h_totalB, lengthB],                regionDonorB,     tol=1.0e-18)
 
    # specifiy outer regions
    # metal interfaces
    bfacemask!(gridB, [0.0, lengthB],  [h_pdopingB, lengthB], bregionAcceptorB, tol=1.0e-30) 
    bfacemask!(gridB, [h_totalB, 0.0], [h_totalB, lengthB],   bregionDonorB,    tol=1.0e-30) 

    bfacemask!(gridB, [h_pdopingB+h_intrinsicB/2, 0.0], [h_pdopingB+h_intrinsicB/2, lengthB], binnerInterfaceB,  tol=1.0e-30)

    # no flux interfaces [xmin, ymin], [xmax, ymax]
    # @Jürgen: If not defining this, then one of these boundaries is set to 1??
    bfacemask!(gridB, [0.0, 0.0],            [0.0, lengthB],      4,  tol=1.0e-30) 
    bfacemask!(gridB, [h_pdopingB, lengthB], [h_totalB, lengthB], 4,  tol=1.0e-30) 
    bfacemask!(gridB, [0.0, 0.0],            [h_totalB, 0.0],     4,  tol=1.0e-30) 

    gridplot(gridB, Plotter=PyPlot, legend=:rt)
end

# ╔═╡ b7a40324-dea5-4fd8-9d7a-1031dceee305
md"""
## Solving the two test problems

"""

# ╔═╡ 905fd6fe-4858-4d85-aeba-7e57439531f5
begin
	ispec_flux       = 1
 
    sysA = VoronoiFVM.System(gridA, VoronoiFVM.Physics(flux=fluxA!,  reaction=reaction!, source=source!, bflux=bfluxA!))
    sysB = VoronoiFVM.System(gridB, VoronoiFVM.Physics(flux=fluxB!,  reaction=reaction!, source=source!, bflux=bfluxB!))
 
    # enable species in all regions 
    enable_species!(sysA, ispec_flux, regionsA)
    enable_species!(sysB, ispec_flux, regionsB)
 
    # boundary conditions
    boundary_dirichlet!(sysA, ispec_flux, bregionAcceptorA, 1.0)
    boundary_dirichlet!(sysA, ispec_flux, bregionDonorA,    0.0)
    ####
    boundary_dirichlet!(sysB, ispec_flux, bregionAcceptorB, 1.0)
    boundary_dirichlet!(sysB, ispec_flux, bregionDonorB,    0.0)
 
    inivalA  = unknowns(sysA);    inivalA .= 0.0;    solA     = unknowns(sysA)
    inivalB  = unknowns(sysB);    inivalB .= 0.0;    solB     = unknowns(sysB)
 
    ## Create solver control info
    control   = VoronoiFVM.NewtonControl()
    
    ## Stationary solution of both problems
    solve!(solA, inivalA, sysA, control=control)
    solve!(solB, inivalB, sysB, control=control)
end

# ╔═╡ 80cfa231-23b1-48a3-a5ad-ee0e122da746
begin
	p2 = GridVisualizer(;Plotter=PyPlot, dim=2, layout=(1, 1), clear=false, resolution=(800, 500))
    # this is for variable transformation, since we consider right outer boundary and want to transform to x-axis.
    function tran32!(a,b)
        a[1] = b[2]
    end

    # note that if adjusting active_boundary to 3 or 4, then transform needs to be deleted.
    bgridA      = subgrid(gridA, [3], boundary = true, transform = tran32!)
    bgridB      = subgrid(gridB, [binnerInterfaceB], boundary = true, transform = tran32!)

    sol_boundA = view(solA[ispec, :], bgridA)
    sol_boundB = view(solB[ispec, :], bgridB)


    scalarplot!(p2[1,1], bgridA, sol_boundA, show = true, clear=true, flimits=(0.1,0.4), cellwise = true, color=:green, label = "Problem A ")
    scalarplot!(p2[1,1], bgridB, sol_boundB, show = true, clear=false, cellwise = true, color=:red, label = "Problem B ", legend=:lt)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ExtendableGrids = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
GridVisualize = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
PyPlot = "d330b81b-6aea-500a-939a-2ce795aea3ee"
SimplexGridFactory = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
Triangulate = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
VoronoiFVM = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"

[compat]
ExtendableGrids = "~0.9.10"
GridVisualize = "~0.5.1"
PlutoUI = "~0.7.39"
PyPlot = "~2.10.0"
SimplexGridFactory = "~0.5.16"
Triangulate = "~2.1.3"
VoronoiFVM = "~0.16.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-rc1"
manifest_format = "2.0"
project_hash = "84d95e7b5fc2144cf83436eefeacd49837b8313d"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "dd2f52bc149ff35158827471453e2e4f1a2685a6"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.26.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "d956c0606a3bc1112a1f99a8b2309b79558d9921"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.17"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d618d3cf75e8ed5064670e939289698ecf426c7f"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.12"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "d7dc30474e73173a990eca86af76cae8790fa9f2"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "b15a6bc52594f5e4a3b825858d1089618871bf9d"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.36"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "7297381ccb5df764549818d9a7d57e45f1057d30"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.18.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "332a332c97c7071600984b3c31d9067e1a4e6e25"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

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

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "28d605d9a0ac17118fe2c5e9ce0fbb76c3ceb120"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "0ec161f87bf4ab164ff96dfacf4be8ffff2375fd"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.62"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "c6c5eb6c4fde80d1a90e5a7e05cf2adfb14e3706"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.10"

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
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.ElasticArrays]]
deps = ["Adapt"]
git-tree-sha1 = "a0fcc1bb3c9ceaf07e1d0529c9806ce94be6adf9"
uuid = "fdbdab4c-e67f-52f5-8c3f-e7b388dad3d4"
version = "1.2.9"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.ExtendableGrids]]
deps = ["AbstractTrees", "Dates", "DocStringExtensions", "ElasticArrays", "InteractiveUtils", "LinearAlgebra", "Printf", "Random", "SparseArrays", "StaticArrays", "Test", "WriteVTK"]
git-tree-sha1 = "1dcc8c9a6a480bae1974b74f5a0d4af8dbfada0e"
uuid = "cfc395e8-590f-11e8-1f13-43a2532b2fa8"
version = "0.9.10"

[[deps.ExtendableSparse]]
deps = ["DocStringExtensions", "LinearAlgebra", "Printf", "Requires", "SparseArrays", "SuiteSparse", "Test"]
git-tree-sha1 = "eb3393e4de326349a4b5bccd9b17ed1029a2d0ca"
uuid = "95c220a8-a1cf-11e9-0c77-dbfce5f500b3"
version = "0.6.7"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ee13c773ce60d9e95a6c6ea134f25605dce2eda3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.13.0"

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
git-tree-sha1 = "2f18915445b248731ec5db4e4a17e451020bf21e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.30"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArrays]]
deps = ["Adapt", "LLVM", "LinearAlgebra", "Printf", "Random", "Serialization", "Statistics"]
git-tree-sha1 = "c783e8883028bf26fb05ed4022c450ef44edd875"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.3.2"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "4888af84657011a65afc7a564918d281612f983a"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.0"

[[deps.GridVisualize]]
deps = ["ColorSchemes", "Colors", "DocStringExtensions", "ElasticArrays", "ExtendableGrids", "GeometryBasics", "HypertextLiteral", "LinearAlgebra", "Observables", "OrderedCollections", "PkgVersion", "Printf", "StaticArrays"]
git-tree-sha1 = "5d845bccf5d690879f4f5f01c7112e428b1fa543"
uuid = "5eed8a63-0fb0-45eb-886d-8d5a387d12b8"
version = "0.5.1"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "7ac3333b82b85f753dd22bed79b85cabcd5e7317"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.7"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "cb7099a0109939f16a4d3b572ba8396b1f6c7c31"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.10"

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

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Statistics"]
git-tree-sha1 = "ad841eddfb05f6d9be0bff1fa48dcae32f134a2d"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.6.2"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

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
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

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

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "e7e9184b0bf0158ac4e4aa9daf00041b5909bf1a"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.14.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg", "TOML"]
git-tree-sha1 = "771bfe376249626d3ca12bcd58ba243d3f961576"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.16+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ArrayInterfaceStaticArrays", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "a63da17ff71f41a1f818e0e1d3c02a32cf4c51f7"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.10.2"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.81.0+0"

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
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "e129d9391168c677cd4800f5c0abb1ed8cb3794f"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

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
version = "2.28.0+0"

[[deps.MeshIO]]
deps = ["ColorTypes", "FileIO", "GeometryBasics", "Printf"]
git-tree-sha1 = "8be09d84a2d597c7c0c34d7d604c039c9763e48c"
uuid = "7269a6da-0436-5bbc-96c2-40638cbb6118"
version = "0.4.10"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "a160e323d3684889e6026914576f1f4288de131d"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.4"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "4e675d6e9ec02061800d6cfb695812becbd03cdf"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.4"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "fe29afdef3d0c4a8286128d4e45cc50621b1e43d"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.4.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

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
git-tree-sha1 = "7f4869861f8dac4990d6808b66b57e5a425cfd99"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.13"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "747f4261ebe38a2bc6abf0850ea8c6d9027ccd07"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "1fc929f47d7c151c839c5fc1375929766fb8edcc"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.93.1"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "14c1b795b9d764e1784713941e787e1384268103"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.10.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

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
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArrays", "LinearAlgebra", "RecipesBase", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "c8bb13a16838ce37f94149c356c5664562b46548"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.29.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "ac248d767048e681843ab674b18e483b05bedc09"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.41.2"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimplexGridFactory]]
deps = ["DocStringExtensions", "ElasticArrays", "ExtendableGrids", "FileIO", "GridVisualize", "LinearAlgebra", "MeshIO", "Printf", "Test"]
git-tree-sha1 = "7fc289795d2eb1ccd8c7917fa346289c95639f04"
uuid = "57bfcd06-606e-45d6-baf4-4ba06da0efd5"
version = "0.5.16"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "32025c052719c6353f22f7c6de7d7b97b7cd2c88"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.24.0"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "5d2c08cef80c7a3a8ba9ca023031a85c263012c5"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2bbd9f2e40afd197a1379aef05e0d85dba649951"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2c11d7290036fe7aac9038ff312d3b3a2a5bf89e"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.4.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "9abba8f8fb8458e9adf07c8a2377a070674a24f1"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.8"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "8bd05747214d9f085ab13d94d2d298111aefe39b"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.9"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "3229d0c2dde1669f239430efece4c540a23bca01"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.7.0"

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
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "d223de97c948636a4f34d1f84d92fd7602dc555b"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.10"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "464d64b2510a25e6efe410e7edab14fffdc333df"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.20"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c76399a3bbe6f5a88faa33c8f8a65aa631d95013"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.73"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Triangle_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bfdd9ef1004eb9d407af935a6f36a4e0af711369"
uuid = "5639c1d2-226c-5e70-8d55-b3095415a16a"
version = "1.6.1+0"

[[deps.Triangulate]]
deps = ["DocStringExtensions", "Libdl", "Printf", "Test", "Triangle_jll"]
git-tree-sha1 = "796a9c0b02a3414af6065098bb7cf0e88dfa450e"
uuid = "f7e6ffb2-c36d-4f8f-a77e-16e897189344"
version = "2.1.3"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.VoronoiFVM]]
deps = ["DiffResults", "DocStringExtensions", "ExtendableGrids", "ExtendableSparse", "ForwardDiff", "GridVisualize", "IterativeSolvers", "JLD2", "LinearAlgebra", "Parameters", "Printf", "RecursiveArrayTools", "Requires", "SparseArrays", "SparseDiffTools", "StaticArrays", "Statistics", "SuiteSparse", "Symbolics", "Test"]
git-tree-sha1 = "f3a434d7621f21688a273f774031cef4bc3baed3"
uuid = "82b139dc-5afc-11e9-35da-9b9bdfd336f3"
version = "0.16.4"

[[deps.WriteVTK]]
deps = ["Base64", "CodecZlib", "FillArrays", "LightXML", "TranscodingStreams"]
git-tree-sha1 = "bff2f6b5ff1e60d89ae2deba51500ce80014f8f6"
uuid = "64499a7a-5c06-52f2-abe2-ccb03c286192"
version = "1.14.2"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

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
version = "5.1.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.41.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═60941eaa-1aea-11eb-1277-97b991548781
# ╟─cc49b312-7e21-49aa-b42d-32c184919aa2
# ╠═b61dce01-f84f-4e46-be76-8561bf0dddaa
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
# ╠═b16061ae-8a08-4649-af8e-58ddf054996d
# ╟─3d0eebe4-5eaa-4339-b6f6-bb5acd573033
# ╠═d7910827-238a-4b71-ad44-7ce43ff1c0db
# ╠═7d902e19-4ced-4685-ac5b-129ce101c9fe
# ╠═7b51cf5a-fa12-4756-a016-ed00cf90b1ba
# ╠═3b429e0c-c07a-40df-ab0b-1f71d515ece3
# ╠═c770801b-d8f8-4095-98ff-e546f661e141
# ╠═fd59fa72-bc40-41ef-bdac-485b171da03e
# ╠═0875ce0d-a22b-4e08-9b58-57673c986a67
# ╠═a9d85817-df56-48ea-8aee-f700dc6310c8
# ╠═b7a40324-dea5-4fd8-9d7a-1031dceee305
# ╠═905fd6fe-4858-4d85-aeba-7e57439531f5
# ╠═80cfa231-23b1-48a3-a5ad-ee0e122da746
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
