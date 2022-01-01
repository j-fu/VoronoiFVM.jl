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

    eval_and_assemble(factory.tfsystem,u,u,f,Inf,Inf,0.0,zeros(0))

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
function integrate(system::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractMatrix{Tv},
                   Uold::AbstractMatrix{Tv}, tstep::Real) where {Tv,Ti}
    grid=system.grid
    nspecies=num_species(system)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    src=zeros(Tv,nspecies)
    stor=zeros(Tv,nspecies)
    storold=zeros(Tv,nspecies)
    tstepinv=1.0/tstep

    physics=system.physics
    node=Node{Tv,Ti}(system)
    bnode=BNode{Tv,Ti}(system)
    edge=Edge{Tv,Ti}(system)
    bedge=Edge{Tv,Ti}(system)
    @create_physics_wrappers(physics,node,bnode,edge,bedge)

    UKL=Array{Tv,1}(undef,2*nspecies)
    UK=Array{Tv,1}(undef,nspecies)
    UKold=Array{Tv,1}(undef,nspecies)
    

    geom=grid[CellGeometries][1]
    
    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            _fill!(edge,iedge,icell)

            @views UKL[1:nspecies].=U[:,edge.node[1]]
            @views UKL[nspecies+1:2*nspecies].=U[:,edge.node[2]]

            res.=zero(Tv)
            fluxwrap(res,UKL)
            for ispec=1:nspecies
                if system.node_dof[ispec,edge.node[1]]==ispec && system.node_dof[ispec,edge.node[2]]==ispec
                    integral[ispec]+=system.celledgefactors[iedge,icell]*res[ispec]*(tf[edge.node[1]]-tf[edge.node[2]])
                end
            end
        end
        
        for inode=1:num_nodes(geom)
            _fill!(node,inode,icell)
            begin
                res.=zero(Tv)
                stor.=zero(Tv)
                storold.=zero(Tv)
                for ispec=1:nspecies
                    UK[ispec]=U[ispec,node.index]
                    UKold[ispec]=Uold[ispec,node.index]
                end
                reactionwrap(res,UK)
                sourcewrap(src)
                storagewrap(stor,UK)
                storagewrap(storold,UKold)
            end
            for ispec=1:nspecies
                if system.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=system.cellnodefactors[inode,icell]*(res[ispec]-src[ispec]+(stor[ispec]-storold[ispec])*tstepinv)*tf[node.index]
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
integrate(system::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractMatrix{Tv}) where {Tv,Ti} =integrate(system,tf,U,U,Inf)



############################################################################
"""
$(SIGNATURES)

Steady state part of test function integral.
"""
function integrate_stdy(system::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractArray{Tu,2}) where {Tu,Tv,Ti}
    grid=system.grid
    nspecies=num_species(system)
    integral=zeros(Tu,nspecies)
    res=zeros(Tu,nspecies)
    src=zeros(Tu,nspecies)
    stor=zeros(Tu,nspecies)

    physics=system.physics
    node=Node{Tv,Ti}(system)
    bnode=BNode{Tv,Ti}(system)
    edge=Edge{Tv,Ti}(system)
    bedge=BEdge{Tv,Ti}(system)
    @create_physics_wrappers(physics,node,bnode,edge,bedge)

    UKL=Array{Tu,1}(undef,2*nspecies)
    UK=Array{Tu,1}(undef,nspecies)
    geom=grid[CellGeometries][1]
   
    for icell=1:num_cells(grid)
        for iedge=1:num_edges(geom)
            _fill!(edge,iedge,icell)

            @views UKL[1:nspecies].=U[:,edge.node[1]]
            @views UKL[nspecies+1:2*nspecies].=U[:,edge.node[2]]
            res.=zero(Tv)
            
            fluxwrap(res,UKL)
            for ispec=1:nspecies
                if system.node_dof[ispec,edge.node[1]]==ispec && system.node_dof[ispec,edge.node[2]]==ispec
                    integral[ispec]+=system.celledgefactors[iedge,icell]*res[ispec]*(tf[edge.node[1]]-tf[edge.node[2]])
                end
            end
        end
        
        for inode=1:num_nodes(geom)
            _fill!(node,inode,icell)

            res.=zeros(Tv)
            src.=zeros(Tv)
            @views UK.=U[:,node.index]
            reactionwrap(res,UK)
            sourcewrap(src)
            for ispec=1:nspecies
                if system.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=system.cellnodefactors[inode,icell]*(res[ispec]-src[ispec])*tf[node.index]
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
function integrate_tran(system::AbstractSystem{Tv,Ti},tf::Vector{Tv},U::AbstractArray{Tu,2}) where {Tu,Tv,Ti}
    grid=system.grid
    nspecies=num_species(system)
    integral=zeros(Tu,nspecies)
    res=zeros(Tu,nspecies)
    stor=zeros(Tu,nspecies)


    physics=system.physics
    node=Node{Tv,Ti}(system)
    bnode=BNode{Tv,Ti}(system)
    edge=Edge{Tv,Ti}(system)
    bedge=BEdge{Tv,Ti}(system)
    @create_physics_wrappers(physics,node,bnode,edge,bedge)
    
    UK=Array{Tu,1}(undef,nspecies)
    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    
    for icell=1:num_cells(grid)
        for inode=1:num_nodes(geom)
            _fill!(node,inode,icell)

            res.=zeros(Tv)
            stor.=zeros(Tv)
            @views UK.=U[:,node.index]

            storagewrap(stor,U[:,node.index])
            for ispec=1:nspecies
                if system.node_dof[ispec,node.index]==ispec
                    integral[ispec]+=system.cellnodefactors[inode,icell]*stor[ispec]*tf[node.index]
                end
            end
        end
    end
    return integral
end



# function checkdelaunay(grid)
#     nreg=num_cellregions(grid)
#     D=ones(nreg)
    

#     physics=Physics( 
#         flux=function(f,u,edge)
#         f[1]=Du[1]-u[2]
#         end,
#         storage=function(f,u,node)
#         f[1]=0
#         end
#     )
#     tfsystem=System(system.grid,physics,unknown_storage=:dense)
#     enable_species!(tfsystem,1,[i for i=1:num_cellregions(system.grid)])
#     return TestFunctionFactory{Tv}(system,tfsystem)
# end
