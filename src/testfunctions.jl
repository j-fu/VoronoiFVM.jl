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
function integrate(this::AbstractSystem{Tv},tf::Vector{Tv},U::AbstractArray{Tv}) where Tv
    grid=this.grid
    nspecies=num_species(this)
    integral=zeros(Tv,nspecies)
    res=zeros(Tv,nspecies)
    node=Node{Tv}()
    edge=Edge{Tv}()
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    edge_cutoff=1.0e-12
    K1::Int32=1
    KN::Int32=nspecies
    L1::Int32=nspecies+1
    LN::Int32=2*nspecies

    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        edge.region=reg_cell(grid,icell)
        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end
            @views begin
                K=celledgenode(grid,1,iedge,icell)
                L=celledgenode(grid,2,iedge,icell)
                edge.index=iedge
                edge.nodeK=K
                edge.nodeL=L
                edge.coordL=nodecoord(grid,L)
                edge.coordK=nodecoord(grid,K)
                
                #Set up argument for fluxwrap
                this.physics.flux(this.physics,edge,res,U[:,K], U[:,L])
                for ispec=1:nspecies
                    if this.node_dof[ispec,K]==ispec && this.node_dof[ispec,L]==ispec
                        integral[ispec]+=edge_factors[iedge]*res[ispec]*(tf[K]-tf[L])
                    end
                end
                
            end
        end
        
        for inode=1:num_nodes_per_cell(grid)
            K=cellnode(grid,inode,icell)
            node.index=K
            node.coord=nodecoord(grid,K)
            this.physics.reaction(this.physics,node,res,U[:,K])
            for ispec=1:nspecies
                if this.node_dof[ispec,K]==ispec
                    integral[ispec]+=node_factors[inode]*res[ispec]*tf[K]
                end
            end
        end
    end
    return integral
end

