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
    node=Node(system)
    bnode=BNode(system)
    edge=Edge(system)
    bedge=BEdge(system)
    @create_physics_wrappers(physics,node,bnode,edge,bedge)

    UKL=Array{Tu,1}(undef,2*nspecies)
    geom=grid[CellGeometries][1]

    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            _fill!(edge,iedge,icell)
            K=edge.node[1]
            L=edge.node[2]
            fac=system.celledgefactors[iedge,icell]
            @views UKL[1:nspecies].=U[:,edge.node[1]]
            @views UKL[nspecies+1:2*nspecies].=U[:,edge.node[2]]
            edgeflux.=zero(Tv)
            fluxwrap(edgeflux,UKL)
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
##############################################################################################



"""
    This shall evaluate at a given solution value to an array of type Tp (type of parameters)
    We expect that in this way, parameter derivatives are easily created.

    For a stationary problem, assume we have
    F(u,p)=0

    Then also  F_u(u,p) u + F_p(u,p) p = 0
    If we want u_p for further use, then 
    F_u(u,p) u_p + F_p(u,p) =0

    But actually, we always want to evaluate functionals of the solution.
    ---> need to have a better plan here...

    ... the current trial is pointless as it involves global Jacobin evaluation which we want to avoid.
    We would need to have different flux wrappers etc,

    function reactionwrap(Y,P)
         node.params=P
         reaction(Y,UK,node)
    end
    
    
"""
function _eval(system::AbstractSystem{Tv, Tc, Ti, Tm},
               U::AbstractMatrix{Tv}, # Actual solution iteration
               UOld::AbstractMatrix{Tv}, # Old timestep solution
               F::AbstractMatrix{Tp},# Right hand side: 
               time,
               tstep,# time step size. Inf means stationary solution
               embedparam,
               params::Vector{Tp};
               edge_cutoff=0.0
               ) where {Tv, Tc, Ti, Tm, Tp}
    

    _complete!(system) # needed here as well for test function system which does not use newton
    
    grid    = system.grid
    physics = system.physics

    grid    = system.grid
    physics = system.physics
    node    = Node(system,time,embedparam,params)
    bnode   = BNode(system,time,embedparam,params)
    edge    = Edge(system,time,embedparam,params)
    bedge   = BEdge(system,time,embedparam,params)



    
    @create_physics_wrappers(physics, node, bnode, edge, bedge)

    nspecies::Int = num_species(system)
    matrix        = system.matrix
    
    cellnodefactors::Array{Tv,2}  = system.cellnodefactors
    celledgefactors::Array{Tv,2}  = system.celledgefactors
    bfacenodefactors::Array{Tv,2} = system.bfacenodefactors
    bfaceedgefactors::Array{Tv,2} = system.bfaceedgefactors    
    
    # Reset rhs
    F   .= zero(Tp)
    
    # Arrays holding function results
    res_react= Array{Tp,1}(undef,nspecies)
    res_stor= Array{Tp,1}(undef,nspecies)
    oldstor= Array{Tp,1}(undef,nspecies)
    res_flux= Array{Tp,1}(undef,nspecies)

    
    # Arrays for gathering solution data
    UK    = Array{Tv,1}(undef,nspecies)
    UKOld = Array{Tv,1}(undef,nspecies)
    UKL   = Array{Tv,1}(undef,2*nspecies)

    
    # array holding source term
    src   = zeros(Tv,nspecies)
    bsrc  = zeros(Tv,nspecies)
    # arrays holding storage terms for old solution
    oldstor   = zeros(Tv, nspecies)
    res_react = zeros(Tv, nspecies)
    jac_react = zeros(Tv, nspecies, nspecies)
    oldbstor  = zeros(Tv, nspecies)

    
    
    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here for the
    # case of stationary problems.
    tstepinv = 1.0/tstep 

    
    boundary_factors::Array{Tv,2} = system.boundary_factors
    boundary_values::Array{Tv,2}  = system.boundary_values

    hasbc = !iszero(boundary_factors) || !iszero(boundary_values)

    bfaceregions::Vector{Ti} = grid[BFaceRegions]
    cellregions::Vector{Ti} = grid[CellRegions]

    nbfaces = num_bfaces(grid)
    ncells  = num_cells(grid)
    geom=grid[CellGeometries][1]
    bgeom   = grid[BFaceGeometries][1]

    
    nn::Int = num_nodes(geom)
    ne::Int = num_edges(geom)
    
    
    for icell=1:ncells
        ireg=cellregions[icell]
        for inode=1:nn
        _fill!(node,inode,icell)
            @views UK[1:nspecies]    .= U[:,node.index]
            @views UKOld[1:nspecies] .= UOld[:,node.index]

            if issource
                sourcewrap(src)
            end
            
            if isreaction
                res_react .= zero(Tv)
                reactionwrap(res_react,UK)
            end
            
            storagewrap(res_stor,UK)
            
            # Evaluate storage term for old timestep
            oldstor.=zero(Tv)
            storagewrap(oldstor,UKOld)
            K=node.index
            for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idof,K)
                if system.region_species[ispec,ireg]<=0
                    continue
                end
                _add(F,idof,cellnodefactors[inode,icell]*(res_react[ispec]-src[ispec] + (res_stor[ispec]-oldstor[ispec])*tstepinv))
            end
        end
        
        for iedge=1:ne
            if abs(celledgefactors[iedge,icell])<edge_cutoff
                continue
            end
            _fill!(edge,iedge,icell)
            
            #Set up argument for fluxwrap
            @views UKL[1:nspecies]            .= U[:, edge.node[1]]
            @views UKL[nspecies+1:2*nspecies] .= U[:, edge.node[2]]
            
            fluxwrap(res_flux,UKL)
            
            K=edge.node[1]
            L=edge.node[2]
            
            fac=celledgefactors[iedge,icell]
            for idofK=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idofK,K)
                idofL=dof(F,ispec,L)
                if idofL==0
                    continue
                end
                if system.region_species[ispec,ireg]<=0
                    continue
                end
                
                _add(F,idofK,fac*res[ispec])
                _add(F,idofL,-fac*res[ispec])
            end
        end
    end
    
    nbn::Int = num_nodes(bgeom)
    nbe::Int = num_edges(bgeom)
    
    # Assembly loop for boundary conditions
    for ibface=1:nbfaces
        ibreg=bfaceregions[ibface]
        
        # Loop over nodes of boundary face
        for ibnode=1:nbn
            # Fill bnode data shuttle with data from grid
            _fill!(bnode,ibnode,ibface)
            
            
            # Copy unknown values from solution into dense array
            @views UK[1:nspecies].=U[:,bnode.index]
            
            # Measure of boundary face part assembled to node
            bnode_factor::Tv=bfacenodefactors[ibnode,ibface]
            bnode.Dirichlet=Dirichlet/bnode_factor
            
            if isbsource
                bsourcewrap(bsrc)
            end
            
            # Global index of node
            K::Int=bnode.index
            
            
            if hasbc
                # Assemble "standard" boundary conditions: Robin or
                # Dirichlet
                # valid only for interior species, currently not checked
                for ispec=1:nspecies
                    idof=dof(F,ispec,K)
                    # If species is present, assemble the boundary condition
                    if idof>0
                        # Get user specified data
                        boundary_factor = boundary_factors[ispec,ibreg]
                        boundary_value  = boundary_values[ispec,ibreg]
                        
                        if boundary_factor==Dirichlet
                            # Dirichlet is encoded in the boundary condition factor
                            # Penalty method: scale   (u-boundary_value) by penalty=boundary_factor
                            
                            # Add penalty*boundry_value to right hand side
                            F[ispec,K]+=boundary_factor*(U[ispec,K]-boundary_value)
                        else
                            # Robin boundary condition
                            F[ispec,K]+=bnode_factor*(boundary_factor*U[ispec,K]-boundary_value)
                        end
                    end
                end
            end
            
            # Boundary reaction term has been given.
            # Valid for both boundary and interior species
            if isbreaction
                # Evaluate function + derivative
                Y.=zero(Tv)
                ForwardDiff.vector_mode_jacobian!(result_r,breactionwrap,Y,UK,cfg_br)
                res_breact=DiffResults.value(result_r)
                jac_breact=DiffResults.jacobian(result_r)
                
                # Assemble RHS + matrix
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    if !isdof(system,ispec,K)
                        continue
                    end
                    _add(F,idof,bnode_factor*(res_breact[ispec]-bsrc[ispec]))
                end
                
            end
            
            # Boundary reaction term has been given.
            # Valid only for boundary species, but system is currently not checked.
            if isbstorage
            
                # Fetch data for old timestep
                @views UKOld.=UOld[:,bnode.index]
                
                # Evaluate storage term for old timestep
                bstoragewrap(oldbstor,UKOld)
                
                # Evaluate & differentiate storage term for new time step value
                Y.=zero(Tv)
                ForwardDiff.vector_mode_jacobian!(result_s,bstoragewrap,Y,UK,cfg_bs)
                res_bstor=DiffResults.value(result_s)
                jac_bstor=DiffResults.jacobian(result_s)
                
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    if !isdof(system,ispec,K)
                        continue
                    end
                    
                    # Assemble finite difference in time for right hand side
                    _add(F,idof,bnode_factor*(res_bstor[ispec]-oldbstor[ispec])*tstepinv)
                end
            end # if isbstorage
            
        end # ibnode=1:nbn

        if isbflux
            for ibedge=1:nbe
                if abs(bfaceedgefactors[ibedge,ibface]) < edge_cutoff
                    continue
                end
                
                _fill!(bedge,ibedge,ibface)
                @views UKL[1:nspecies]            .= U[:, bedge.node[1]]
                @views UKL[nspecies+1:2*nspecies] .= U[:, bedge.node[2]]
                
                Y.=zero(Tv)
                ForwardDiff.vector_mode_jacobian!(result_bflx, bfluxwrap, Y, UKL, cfg_bflx)
                res = DiffResults.value(result_bflx)
                jac = DiffResults.jacobian(result_bflx)
                
                K   = bedge.node[1]
                L   = bedge.node[2]
            
                fac = bfaceedgefactors[ibedge, ibface]
                
                for idofK = _firstnodedof(F, K):_lastnodedof(F, K)
                    ispec =_spec(F, idofK, K)
                    if !isdof(system,ispec,K)
                        continue
                    end
                    
                    idofL = dof(F, ispec, L)
                    if idofL == 0
                        continue
                    end
                    _add(F, idofK,  fac*res[ispec])
                    _add(F, idofL, -fac*res[ispec])
                end
            end
        end
     end

    _eval_generic_operator(system,U,F)
    _eval_inactive_species(system,U,UOld,F)
end


