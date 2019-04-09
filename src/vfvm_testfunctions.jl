################################################
"""
    mutable struct TestFunctionFactory{Tv}

Data structure containing DenseSystem used to calculate
test functions for boundary flux calculations.

"""
mutable struct TestFunctionFactory{Tv}
    system::AbstractSystem{Tv}
    tfsystem::DenseSystem{Tv}
end

#
# Physics data for testfunction system
#
mutable struct TestFunctionPhysics <: Physics
    storage
    flux
    reaction
    TestFunctionPhysics()=new()
end

################################################
"""
    function TestFunctionFactory(system::AbstractSystem{Tv}) where Tv

Constructor for TestFunctionFactory,
"""
function TestFunctionFactory(system::AbstractSystem{Tv}) where Tv
    physics=TestFunctionPhysics()
    physics.flux=function(physics,edge,f,uk,ul)
        f[1]=uk[1]-ul[1]
    end
    physics.storage=function(physics,node,f,u)
        f[1]=u[1]
    end
    physics.reaction=function(physics,node,f,u)
    end
    tfsystem=DenseSystem(system.grid,physics,1)
    add_species(tfsystem,1,[i for i=1:num_cellregions(system.grid)])
    return TestFunctionFactory{Tv}(system,tfsystem)
end



############################################################################
"""
    function testfunction(factory::TestFunctionFactory{Tv}, bc0, bc1) where Tv

Create testfunction which as Dirichlet zero boundary conditions  for boundary
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

    _eval_and_assemble(factory.tfsystem,u,u,f,Inf)
    inidirichlet!(factory.tfsystem,u)
    lufact=LinearAlgebra.lu(factory.tfsystem.matrix)
    ldiv!(vec(u),lufact,vec(f))
    return vec(u)
end


############################################################################
"""
    function integrate(this::AbstractSystem{Tv},tf::Vector{Tv},U::AbstractArray{Tv}) where Tv

Calculate test function integral.
"""

function integrate(this::AbstractSystem{Tv},tf::Vector{Tv},U::AbstractMatrix{Tv}, Uold::AbstractMatrix{Tv}, tstep::Real) where Tv
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    stor=zeros(Tv,nspecies)
    storold=zeros(Tv,nspecies)
    tstepinv=1.0/tstep
    node=Node{Tv}()
    edge=Edge{Tv}()
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    edge_cutoff=1.0e-12
    
    
    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        edge.region=reg_cell(grid,icell)
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end
            fill!(edge,grid,iedge,icell)
            @views this.physics.flux(this.physics,edge,res,U[:,edge.nodeK], U[:,edge.nodeL])
            for ispec=1:nspecies
                if this.node_dof[ispec,edge.nodeK]==ispec && this.node_dof[ispec,edge.nodeL]==ispec
                    integral[ispec]+=edge_factors[iedge]*res[ispec]*(tf[edge.nodeK]-tf[edge.nodeL])
                end
            end
        end
        
        for inode=1:num_nodes_per_cell(grid)
            fill!(node,grid,inode,icell)
            @views begin
                this.physics.reaction(this.physics,node,res,U[:,node.index])
                this.physics.storage(this.physics,node,stor,U[:,node.index])
                this.physics.storage(this.physics,node,storold,Uold[:,node.index])
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

integrate(this::AbstractSystem{Tv},tf::Vector{Tv},U::AbstractMatrix{Tv}) where Tv=integrate(this,tf,U,U,Inf)



