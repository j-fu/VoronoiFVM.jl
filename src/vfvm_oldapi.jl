#
# Methods for pre-0.7 API.
#

# Main assembly loop
function _eval_and_assemble_oldapi(this::AbstractSystem{Tv,Ti},
                            U::AbstractMatrix{Tv}, # Actual solution iteration
                            UOld::AbstractMatrix{Tv}, # Old timestep solution
                            F::AbstractMatrix{Tv},# Right hand side
                            tstep::Tv, # time step size. Inf means stationary solution
                            xstorage::FSTOR, # ensure type stability
                            xbstorage::FBSTOR, # ensure type stability
                            xsource::FSRC,  # ensure type stability
                           ) where {Tv,Ti,FSRC,FSTOR, FBSTOR}

    _complete!(this) # needed here as well for test function system which does not use newton
    
    grid=this.grid

    physics::Physics=this.physics
    data=physics.data
    node=Node{Tv,Ti}(this)
    bnode=BNode{Tv,Ti}(this)
    edge=Edge{Tv,Ti}(this)
    edge_cutoff=1.0e-12
    nspecies::Int32=num_species(this)
    
    matrix=this.matrix
    issource=(physics.source!=nofunc)
    isreaction=(physics.reaction!=nofunc)
    isbreaction=(physics.breaction!=nofunc)
    isbstorage=(physics.bstorage!=nofunc)

    
    function fluxwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.flux(y,u,edge,data)
    end


    function reactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        ## for ii in ..  uu[node.speclist[ii]]=u[ii]
        physics.reaction(y,u,node,data)
        ## for ii in .. y[ii]=y[node.speclist[ii]]
    end

    function storagewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        xstorage(y,u,node,data)
    end

    function sourcewrap(y)
        y.=0
        xsource(y,node,data)
    end

    @inline function breactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.breaction(y,u,bnode,data)
    end

    @inline function bstoragewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        xbstorage(y,u,bnode,data)
    end
    
    # Reset matrix + rhs
    nzv=nonzeros(matrix)
    nzv.=0.0
    F.=0.0

    # structs holding diff results for storage, reaction,  flux ...
    result_r=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_s=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_flx=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,2*nspecies))
    
    # Array holding function results
    Y=Array{Tv,1}(undef,nspecies)

    # Arrays for gathering solution data
    UK=Array{Tv,1}(undef,nspecies)
    UKOld=Array{Tv,1}(undef,nspecies)
    UKL=Array{Tv,1}(undef,2*nspecies)

    # array holding source term
    src=zeros(Tv,nspecies)

    # arrays holding storage terms for old solution
    oldstor=zeros(Tv,nspecies)
    res_react=zeros(Tv,nspecies)
    jac_react=zeros(Tv,nspecies,nspecies)
    oldbstor=zeros(Tv,nspecies)

    cfg_r=ForwardDiff.JacobianConfig(reactionwrap, Y, UK)
    cfg_s=ForwardDiff.JacobianConfig(storagewrap, Y, UK)
    cfg_br=ForwardDiff.JacobianConfig(breactionwrap, Y, UK)
    cfg_bs=ForwardDiff.JacobianConfig(bstoragewrap, Y, UK)
    cfg_flx=ForwardDiff.JacobianConfig(fluxwrap, Y, UKL)

    
    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here.
    tstepinv=1.0/tstep 

    
    # Arrays holding for factor data
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    bnode_factors=zeros(Tv,num_nodes_per_bface(grid))

    # Main cell loop
    for icell=1:num_cells(grid)
        # set up form factors
        cellfactors!(grid,icell,node_factors,edge_factors)
        # set up data for callbacks

        node.region=reg_cell(grid,icell)
        edge.region=reg_cell(grid,icell)

        for inode=1:num_nodes_per_cell(grid)
            _fill_oldapi!(node,grid, inode,icell)
            for ispec=1:nspecies
                # xx gather:
                # ii=0
                # for i region_spec.colptr[ireg]:region_spec.colptr[ireg+1]-1
                #    ispec=Fdof.rowval[idof]
                #    ii=ii+1
                #    node.speclist[ii]=ispec
                #    UK[ii]=U[ispec,K]
                #   UK[1:nspecies]=U[:,node.index]
                #   UKOld[1:nspecies]=UOld[:,node.index]
                # Evaluate source term
                UK[ispec]=U[ispec,node.index]
                UKOld[ispec]=UOld[ispec,node.index]
            end


            if issource
                sourcewrap(src)
            end

            if isreaction
                # Evaluate & differentiate reaction term if present
                ForwardDiff.jacobian!(result_r,reactionwrap,Y,UK,cfg_r)
                res_react=DiffResults.value(result_r)
                jac_react=DiffResults.jacobian(result_r)
            end

            # Evaluate & differentiate storage term
            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK,cfg_s)
            res_stor=DiffResults.value(result_s)
            jac_stor=DiffResults.jacobian(result_s)

            # Evaluate storage term for old timestep
            ### ALLOC ???
            storagewrap(oldstor,UKOld) 
            
            
            K=node.index

            for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idof,K)
                _add(F,idof,node_factors[inode]*(res_react[ispec]-src[ispec] + (res_stor[ispec]-oldstor[ispec])*tstepinv))
                for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                    jspec=_spec(F,jdof,K)
                    _addnz(matrix,idof,jdof,jac_react[ispec,jspec]+ jac_stor[ispec,jspec]*tstepinv,node_factors[inode])
                end
            end
        end
        
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end

            _fill_oldapi!(edge,grid,iedge,icell)

            #Set up argument for fluxwrap
            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.node[1]]
                UKL[ispec+nspecies]=U[ispec,edge.node[2]]
            end

            ForwardDiff.jacobian!(result_flx,fluxwrap,Y,UKL,cfg_flx)
            res=DiffResults.value(result_flx)
            jac=DiffResults.jacobian(result_flx)

            K=edge.node[1]
            L=edge.node[2]
            fac=edge_factors[iedge]
            for idofK=_firstnodedof(F,K):_lastnodedof(F,K)
                ispec=_spec(F,idofK,K)
                idofL=dof(F,ispec,L)
                if idofL==0
                    continue
                end
                _add(F,idofK,fac*res[ispec])
                _add(F,idofL,-fac*res[ispec])
                
                for jdofK=_firstnodedof(F,K):_lastnodedof(F,K)
                    jspec=_spec(F,jdofK,K)
                    jdofL=dof(F,jspec,L)
                    if jdofL==0
                        continue
                    end
                    
                    _addnz(matrix,idofK,jdofK,+jac[ispec,jspec         ],fac)
                    _addnz(matrix,idofL,jdofK,-jac[ispec,jspec         ],fac)
                    _addnz(matrix,idofK,jdofL,+jac[ispec,jspec+nspecies],fac)
                    _addnz(matrix,idofL,jdofL,-jac[ispec,jspec+nspecies],fac)
                    
                end
            end
        end
    end

    # Assembly loop for boundary conditions
    for ibface=1:num_bfaces(grid)

        # Calculate measure of boundary face contribution
        # to the corresponding nodes
        bfacefactors!(grid,ibface,bnode_factors)

        # Obtain boundary region number
        ibreg=grid.bfaceregions[ibface]

        # Loop over nodes of boundary face
        for ibnode=1:num_nodes_per_bface(grid)

            # Fill bnode data shuttle with data from grid
            _fill_oldapi!(bnode,grid,ibnode,ibface)

            # Copy unknown values from solution into dense array
            for ispec=1:nspecies
                UK[ispec]=U[ispec,bnode.index]
            end

            # Measure of boundary face part assembled to node
            bnode_factor=bnode_factors[ibnode]

            # Global index of node
            K=bnode.index
            
            # Assemble "standard" boundary conditions: Robin or
            # Dirichlet
            # valid only for interior species, currently not checked
            for ispec=1:nspecies
                idof=dof(F,ispec,K)
                # If species is present, assemble the boundary condition
                if idof>0
                    # Get user specified data
                    boundary_factor=this.boundary_factors[ispec,ibreg]
                    boundary_value=this.boundary_values[ispec,ibreg]
                    
                    if boundary_factor==Dirichlet
                        # Dirichlet is encoded in the boundary condition factor
                        # Penalty method: scale   (u-boundary_value) by penalty=boundary_factor
                        
                        # Add penalty*boundry_value to right hand side
                        F[ispec,K]+=boundary_factor*(U[ispec,K]-boundary_value)
                        
                        # Add penalty to matrix main diagonal (without bnode factor, so penalty
                        # is independent of h)
                        _addnz(matrix,idof,idof,boundary_factor,1)
                    else
                        # Robin boundary condition
                        F[ispec,K]+=bnode_factor*(boundary_factor*U[ispec,K]-boundary_value)
                        _addnz(matrix,idof,idof,boundary_factor, bnode_factor)
                    end
                end
            end

            # Boundary reaction term has been given.
            # Valid for both boundary and interior species
            if isbreaction

                # Evaluate function + derivative
                ForwardDiff.jacobian!(result_r,breactionwrap,Y,UK,cfg_br)
                res_breact=DiffResults.value(result_r)
                jac_breact=DiffResults.jacobian(result_r)
                
                # Assemble RHS + matrix
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    _add(F,idof,bnode_factor*res_breact[ispec])
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        _addnz(matrix,idof,jdof,jac_breact[ispec,jspec],bnode_factor)
                    end
                end

            end
            
            # Boundary reaction term has been given.
            # Valid only for boundary species, but this is currently not checked.
            if isbstorage
                
                # Fetch data for old timestep
                for ispec=1:nspecies
                    UKOld[ispec]=UOld[ispec,K]
                end
                # Evaluate storage term for old timestep
                bstoragewrap(oldbstor,UKOld)

                # Evaluate & differentiate storage term for new time step value
                ForwardDiff.jacobian!(result_s,bstoragewrap,Y,UK,cfg_bs)
                res_bstor=DiffResults.value(result_s)
                jac_bstor=DiffResults.jacobian(result_s)
                
                for idof=_firstnodedof(F,K):_lastnodedof(F,K)
                    ispec=_spec(F,idof,K)
                    # Assemble finite difference in time for right hand side
                    _add(F,idof,bnode_factor*(res_bstor[ispec]-oldbstor[ispec])*tstepinv)
                    # Assemble matrix.
                    for jdof=_firstnodedof(F,K):_lastnodedof(F,K)
                        jspec=_spec(F,jdof,K)
                        _addnz(matrix,idof,jdof,jac_bstor[ispec,jspec],bnode_factor*tstepinv)
                    end
                end
            end
            
        end
    end
    _eval_and_assemble_inactive_species(this,U,UOld,F)
end

# Test function integration
function _integrate_oldapi(this::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractMatrix{Tv}, Uold::AbstractMatrix{Tv}, tstep::Real) where {Tv,Ti}
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    stor=zeros(Tv,nspecies)
    storold=zeros(Tv,nspecies)
    tstepinv=1.0/tstep
    node=Node{Tv,Ti}(this)
    edge=Edge{Tv,Ti}(this)
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    edge_cutoff=1.0e-12
    UKL=Array{Tv,1}(undef,2*nspecies)

    
    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        edge.region=reg_cell(grid,icell)
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end
            _fill_oldapi!(edge,grid,iedge,icell)

            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.node[1]]
                UKL[ispec+nspecies]=U[ispec,edge.node[2]]
            end
            res.=0
            @views this.physics.flux(res,UKL,edge, this.physics.data)
            for ispec=1:nspecies
                if this.node_dof[ispec,edge.node[1]]==ispec && this.node_dof[ispec,edge.node[2]]==ispec
                    integral[ispec]+=edge_factors[iedge]*res[ispec]*(tf[edge.node[1]]-tf[edge.node[2]])
                end
            end
        end
        
        for inode=1:num_nodes_per_cell(grid)
            _fill_oldapi!(node,grid,inode,icell)
            @views begin
                res.=0
                stor.=0
                storold.=0
                this.physics.reaction(res,U[:,node.index],node,this.physics.data)
                this.physics.storage(stor,U[:,node.index],node,this.physics.data)
                this.physics.storage(storold,Uold[:,node.index],node,this.physics.data)
            end
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=node_factors[inode]*(res[ispec]+(stor[ispec]-storold[ispec])*tstepinv)*tf[node.index]
                end
            end
        end
    end
    return integral
end



function _integrate_stdy_oldapi(this::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractArray{Tu,2}) where {Tu,Tv,Ti}
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tu,nspecies)
    res=zeros(Tu,nspecies)
    stor=zeros(Tu,nspecies)
    node=Node{Tv,Ti}(this)
    edge=Edge{Tv,Ti}(this)
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    edge_cutoff=1.0e-12
    UKL=Array{Tu,1}(undef,2*nspecies)
    UK=Array{Tu,1}(undef,nspecies)
    
    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        edge.region=reg_cell(grid,icell)
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end
            _fill_oldapi!(edge,grid,iedge,icell)

            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.node[1]]
                UKL[ispec+nspecies]=U[ispec,edge.node[2]]
                res[ispec]=0.0
            end
            this.physics.flux(res,UKL,edge, this.physics.data)
            for ispec=1:nspecies
                if this.node_dof[ispec,edge.node[1]]==ispec && this.node_dof[ispec,edge.node[2]]==ispec
                    integral[ispec]+=edge_factors[iedge]*res[ispec]*(tf[edge.node[1]]-tf[edge.node[2]])
                end
            end
        end
        
        for inode=1:num_nodes_per_cell(grid)
            _fill_oldapi!(node,grid,inode,icell)
            for ispec=1:nspecies
                res[ispec]=0.0
                UK[ispec]=U[ispec,node.index]
            end
            this.physics.reaction(res,UK,node,this.physics.data)
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=node_factors[inode]*res[ispec]*tf[node.index]
                end
            end
        end
    end
    return integral
end


function _integrate_tran_oldapi(this::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractArray{Tu,2}) where {Tu,Tv,Ti}
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tu,nspecies)
    res=zeros(Tu,nspecies)
    stor=zeros(Tu,nspecies)
    node=Node{Tv,Ti}(this)
    edge=Edge{Tv,Ti}(this)
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    edge_cutoff=1.0e-12
    UK=Array{Tu,1}(undef,nspecies)
    
    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        for inode=1:num_nodes_per_cell(grid)
            _fill_oldapi!(node,grid,inode,icell)
            for ispec=1:nspecies
                UK[ispec]=U[ispec,node.index]
                res[ispec]=0.0
                stor[ispec]=0.0
            end
            this.physics.storage(stor,U[:,node.index],node,this.physics.data)
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=node_factors[inode]*stor[ispec]*tf[node.index]
                end
            end
        end
    end
    return integral
end



##################################################################
"""
$(TYPEDSIGNATURES)
   
Calculate the length of an edge. 
"""
function edgelength(edge::Edge)
    l=0.0
    for i=1:length(edge.coordK)
        d=edge.coordK[i]-edge.coordL[i]
        l=l+d*d
    end
    return sqrt(l)
end

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
@inline viewK(edge::Edge,u::AbstractArray)=@views u[1:edge.nspec]


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
@inline viewL(edge::Edge,u::AbstractArray)=@views u[edge.nspec+1:2*edge.nspec]

"""
$(TYPEDSIGNATURES)

Solution view on first edge node
"""
@inline viewK(nspec,u::AbstractArray)=@views u[1:nspec]


"""
$(TYPEDSIGNATURES)

Solution view on second edge node
"""
@inline viewL(nspec,u::AbstractArray)=@views u[nspec+1:2*nspec]

function _fill_oldapi!(node::Node,grid::Grid,inode,icell)
    _fill!(node,grid,inode,icell)
    for i=1:length(node.coord)
        node.coord[i]=grid.coord[i,node.index]
    end
end


function _fill_oldapi!(node::BNode,grid::Grid,ibnode,ibface)
    _fill!(node,grid,ibnode,ibface)
    for i=1:length(node.coord)
        node.coord[i]=grid.coord[i,node.index]
    end
end


function _fill_oldapi!(edge::Edge,grid::Grid,iedge,icell)
    _fill!(edge,grid,iedge,icell)
    for i=1:length(edge.coordK)
        edge.coordK[i]=grid.coord[i,edge.node[1]]
        edge.coordL[i]=grid.coord[i,edge.node[2]]
    end
end

data(item::AbstractGeometryItem)=item.physics_data
