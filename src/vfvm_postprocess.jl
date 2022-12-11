struct SolutionIntegral{T}<:AbstractMatrix{T} 
    value::Matrix{T}
end

Base.size(I::SolutionIntegral)=size(I.value)
Base.getindex(I::SolutionIntegral,ispec::Integer,ireg)=I.value[ispec,ireg]
Base.setindex!(I::SolutionIntegral,v,ispec::Integer,ireg)=I.value[ispec,ireg]=v



################################################################
"""
````
integrate(system,F,U; boundary=false)    
````

Integrate node function (same signature as reaction or storage)
 `F` of  solution vector region-wise over domain or boundary.
The result is  `nspec x nregion` vector.
"""
function integrate(system::AbstractSystem{Tv,Tc,Ti,Tm},F::Function,U::AbstractMatrix{Tu}; boundary=false) where {Tu,Tv,Tc,Ti,Tm}
    grid=system.grid
    data=system.physics.data
    nspecies=num_species(system)
    res=zeros(Tu,nspecies)

    if boundary
        bnode=BNode(system)
        bnodeparams=(bnode,)
        if isdata(data)
            bnodeparams=(bnode,data,)
        end
#!!!        bnode.time=time
#!!!        bnode.embedparam=embedparam
        
        geom=grid[BFaceGeometries][1]
        bfaceregions=grid[BFaceRegions]
        nbfaceregions=maximum(bfaceregions)
        integral=zeros(Tu,nspecies,nbfaceregions)
        
        for ibface=1:num_bfaces(grid)
            for inode=1:num_nodes(geom)
                _fill!(bnode,inode,ibface)
                res.=zero(Tv)
                @views F(rhs(bnode,res),unknowns(bnode,U[:,bnode.index]),bnodeparams...)
                for ispec=1:nspecies
                    if isdof(system, ispec, bnode.index)
                        integral[ispec,bnode.region]+=system.bfacenodefactors[inode,ibface]*res[ispec]
                    end
                end
            end
        end
    else
        node=Node(system)
        nodeparams=(node,)
        if isdata(data)
            nodeparams=(node,data,)
        end
#!!!        node.time=time
#!!!        node.embedparam=embedparam
    
        
        geom=grid[CellGeometries][1]
        cellnodes=grid[CellNodes]
        cellregions=grid[CellRegions]
        ncellregions=maximum(cellregions)
        integral=zeros(Tu,nspecies,ncellregions)
        
        
        for icell=1:num_cells(grid)
            for inode=1:num_nodes(geom)
                _fill!(node,inode,icell)
                res.=zero(Tv)
                @views F(rhs(node,res),unknowns(node,U[:,node.index]),nodeparams...)
                for ispec=1:nspecies
                    if isdof(system, ispec, node.index)
                        integral[ispec,node.region]+=system.cellnodefactors[inode,icell]*res[ispec]
                    end
                end
            end
        end
    end
    
    return SolutionIntegral(integral)
end



############################################################################
"""
$(SIGNATURES)

Reconstruction of  edge flux as  vector function  on the nodes  of the
triangulation.  The result  can be seen as a  piecewiesw linear vector
function in the FEM space spanned by the discretization nodes exhibiting
the flux density.  

The reconstruction is based on the  "magic formula"
R. Eymard, T. Gallouet, R. Herbin, IMA Journal of Numerical Analysis (2006)
26, 326−353, Lemma 2.4 .

The return value is a `dim x nspec x nnodes` vector. The flux of species i
can  e.g. plotted via GridVisualize.vectorplot.

Example:
```julia
    ispec=3
    vis=GridVisualizer(Plotter=Plotter)
    scalarplot!(vis,grid,solution[ispec,:],clear=true,colormap=:summer)
    vectorplot!(vis,grid,nf[:,ispec,:],clear=false)
    reveal(vis)
```

CAVEAT: there is a possible unsolved problem with the values at domain 
corners in the code. Please see any potential boundary artifacts as a manifestation
of this issue and report them.
"""
function nodeflux(system::AbstractSystem{Tv,Tc,Ti,Tm},U::AbstractArray{Tu,2}) where {Tu,Tv,Tc,Ti,Tm}
    grid=system.grid
    dim=dim_space(grid)
    nnodes=num_nodes(grid)
    nspecies=num_species(system)
    nodeflux=zeros(Tu,dim,nspecies,nnodes)
    edgeflux=zeros(Tu,nspecies)
    xsigma=grid[VoronoiFaceCenters]
    coord=grid[Coordinates]
    nodevol=zeros(Tv,nnodes)
    cellnodes=grid[CellNodes]
    physics=system.physics
    edge=Edge(system)

    # !!! TODO Parameter handling here
    UKL=Array{Tu,1}(undef,2*nspecies)
    flux_eval = ResEvaluator(system.physics,:flux,UKL,edge,nspecies)
    
    geom=grid[CellGeometries][1]

    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            _fill!(edge,iedge,icell)
            K=edge.node[1]
            L=edge.node[2]
            fac=system.celledgefactors[iedge,icell]
            @views UKL[1:nspecies].=U[:,edge.node[1]]
            @views UKL[nspecies+1:2*nspecies].=U[:,edge.node[2]]
            evaluate!(flux_eval,UKL)
            edgeflux = res(flux_eval)
            for ispec=1:nspecies
                if isdof(system, ispec,K) && isdof(system, ispec,L) 
                    @views nodeflux[:,ispec,K].+=fac*edgeflux[ispec]*(xsigma[:,edge.index]-coord[:,K])
                    @views nodeflux[:,ispec,L].-=fac*edgeflux[ispec]*(xsigma[:,edge.index]-coord[:,L])
                end
            end
        end

        for inode=1:num_nodes(geom)
            nodevol[cellnodes[inode,icell]]+=system.cellnodefactors[inode,icell]
        end
    end
    for inode=1:nnodes
        @views nodeflux[:,:,inode]/=nodevol[inode]
    end
    nodeflux
end



