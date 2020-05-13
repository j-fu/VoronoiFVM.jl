################################################
"""
$(TYPEDEF)

Data structure containing DenseSystem used to calculate
test functions for boundary flux calculations.


$(TYPEDFIELDS)
"""
mutable struct TestFunctionFactory{Tv}

    """
    Original system
    """
    system::AbstractSystem{Tv}

    """
    Test function system
    """
    tfsystem::DenseSystem{Tv}
end


################################################
"""
$(TYPEDSIGNATURES)

Constructor for TestFunctionFactory from System
"""
function TestFunctionFactory(system::AbstractSystem{Tv}) where Tv
    physics=Physics( 
        flux=function(f,u,edge)
        f[1]=u[1]-u[2]
        end,
        storage=function(f,u,node)
        f[1]=u[1]
        end
    )
    tfsystem=System(system.grid,physics,unknown_storage=:dense)
    enable_species!(tfsystem,1,[i for i=1:num_cellregions(system.grid)])
    return TestFunctionFactory{Tv}(system,tfsystem)
end



############################################################################
"""
$(TYPEDSIGNATURES)

Create testfunction which has Dirichlet zero boundary conditions  for boundary
regions in bc0 and Dirichlet one boundary conditions  for boundary
regions in bc1.
"""
function testfunction(factory::TestFunctionFactory{Tv}, bc0, bc1) where Tv

    u=unknowns(factory.tfsystem)
    f=unknowns(factory.tfsystem)
    u.=0
    f.=0
    
    factory.tfsystem.boundary_factors.=0
    factory.tfsystem.boundary_values.=0

    for i=1:length(bc1)
        factory.tfsystem.boundary_factors[1,bc1[i]]=Dirichlet
        factory.tfsystem.boundary_values[1,bc1[i]]=-1
    end

    for i=1:length(bc0)
        factory.tfsystem.boundary_factors[1,bc0[i]]=Dirichlet
        factory.tfsystem.boundary_values[1,bc0[i]]=0
    end

    eval_and_assemble(factory.tfsystem,u,u,f,Inf)

    _initialize!(u,factory.tfsystem)
    lufact=LinearAlgebra.lu(factory.tfsystem.matrix)
    ldiv!(vec(u),lufact,vec(f))
    return vec(u)
end


############################################################################

"""
$(SIGNATURES)

Calculate test function integral for transient solution.
"""
function integrate(this::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractMatrix{Tv}, Uold::AbstractMatrix{Tv}, tstep::Real) where {Tv,Ti}
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    stor=zeros(Tv,nspecies)
    storold=zeros(Tv,nspecies)
    tstepinv=1.0/tstep
    node=Node{Tv,Ti}(this)
    edge=Edge{Tv,Ti}(this)
    data=this.physics.data
    UKL=Array{Tv,2}(undef,nspecies,2)

    nodeparams=(node,)
    edgeparams=(edge,)

    if isdata(data)
        nodeparams=(node,data,)
        edgeparams=(edge,data,)
    end    

    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]

    if haskey(grid,CellEdges)
        cellx=grid[CellEdges]
        edgenodes=grid[EdgeNodes]
        has_celledges=true
    else
        cellx=grid[CellNodes]
        edgenodes=local_celledgenodes(geom)
        has_celledges=false
    end
    
    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            _fill!(edge,cellx,edgenodes,cellregions,iedge,icell, has_celledges)

            for ispec=1:nspecies
                UKL[ispec,1]=U[ispec,edge.node[1]]
                UKL[ispec,2]=U[ispec,edge.node[2]]
            end
            res.=0
            @views this.physics.flux(res,UKL,edgeparams...)
            for ispec=1:nspecies
                if this.node_dof[ispec,edge.node[1]]==ispec && this.node_dof[ispec,edge.node[2]]==ispec
                    integral[ispec]+=this.celledgefactors[iedge,icell]*res[ispec]*(tf[edge.node[1]]-tf[edge.node[2]])
                end
            end
        end
        
        for inode=1:num_nodes(geom)
            _fill!(node,cellnodes,cellregions,inode,icell)
            @views begin
                res.=0
                stor.=0
                storold.=0
                this.physics.reaction(res,U[:,node.index],nodeparams...)
                this.physics.storage(stor,U[:,node.index],nodeparams...)
                this.physics.storage(storold,Uold[:,node.index],nodeparams...)
            end
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=this.cellnodefactors[inode,icell]*(res[ispec]+(stor[ispec]-storold[ispec])*tstepinv)*tf[node.index]
                end
            end
        end
    end
    return integral
end

############################################################################
"""
$(SIGNATURES)

Calculate test function integral for steady state solution.
"""
integrate(this::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractMatrix{Tv}) where {Tv,Ti} =integrate(this,tf,U,U,Inf)



############################################################################
"""
$(SIGNATURES)

Steady state part of test function integral.
"""
function integrate_stdy(this::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractArray{Tu,2}) where {Tu,Tv,Ti}
    grid=this.grid
    data=this.physics.data
    nspecies=num_species(this)
    integral=zeros(Tu,nspecies)
    res=zeros(Tu,nspecies)
    stor=zeros(Tu,nspecies)
    node=Node{Tv,Ti}(this)
    edge=Edge{Tv,Ti}(this)
    UKL=Array{Tu,1}(undef,2*nspecies)
    UK=Array{Tu,1}(undef,nspecies)
    nodeparams=(node,)
    edgeparams=(edge,)
    if isdata(data)
        nodeparams=(node,data,)
        edgeparams=(edge,data,)
    end
    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]

    
    if haskey(grid,CellEdges)
        cellx=grid[CellEdges]
        edgenodes=grid[EdgeNodes]
        has_celledges=true
    else
        cellx=grid[CellNodes]
        edgenodes=local_celledgenodes(geom)
        has_celledges=false
    end
        

   
    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            _fill!(edge,cellx,edgenodes,cellregions,iedge,icell, has_celledges)

            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.node[1]]
                UKL[ispec+nspecies]=U[ispec,edge.node[2]]
                res[ispec]=0.0
            end
            this.physics.flux(res,UKL,edgeparams...)
            for ispec=1:nspecies
                if this.node_dof[ispec,edge.node[1]]==ispec && this.node_dof[ispec,edge.node[2]]==ispec
                    integral[ispec]+=this.celledgefactors[iedge,icell]*res[ispec]*(tf[edge.node[1]]-tf[edge.node[2]])
                end
            end
        end
        
        for inode=1:num_nodes(geom)
            _fill!(node,cellnodes,cellregions,inode,icell)
            for ispec=1:nspecies
                res[ispec]=0.0
                UK[ispec]=U[ispec,node.index]
            end
            this.physics.reaction(res,UK,nodeparams...)
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=this.cellnodefactors[inode,icell]*res[ispec]*tf[node.index]
                end
            end
        end
    end
    return integral
end

############################################################################
"""
$(SIGNATURES)

Calculate transient part of test function integral.
"""
function integrate_tran(this::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractArray{Tu,2}) where {Tu,Tv,Ti}
    grid=this.grid
    data=this.physics.data
    nspecies=num_species(this)
    integral=zeros(Tu,nspecies)
    res=zeros(Tu,nspecies)
    stor=zeros(Tu,nspecies)
    node=Node{Tv,Ti}(this)
    edge=Edge{Tv,Ti}(this)

    nodeparams=(node,)
    edgeparams=(edge,)
    if isdata(data)
        nodeparams=(node,data,)
        edgeparams=(edge,data,)
    end
    
    UK=Array{Tu,1}(undef,nspecies)
    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]

    
    for icell=1:num_cells(grid)
        for inode=1:num_nodes(geom)
            _fill!(node,cellnodes,cellregions,inode,icell)
            for ispec=1:nspecies
                UK[ispec]=U[ispec,node.index]
                res[ispec]=0.0
                stor[ispec]=0.0
            end
            this.physics.storage(stor,U[:,node.index],nodeparams...)
            for ispec=1:nspecies
                if this.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=this.cellnodefactors[inode,icell]*stor[ispec]*tf[node.index]
                end
            end
        end
    end
    return integral
end
