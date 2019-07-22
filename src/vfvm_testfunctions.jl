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
$(SIGNATURES)

Constructor for TestFunctionFactory from System
"""
function TestFunctionFactory(system::AbstractSystem{Tv}) where Tv
    physics=Physics( 
        flux=function(f,u,edge,data)
        f[1]=u[1]-u[2]
        end,
        storage=function(f,u,node,data)
        f[1]=u[1]
        end
    )
    tfsystem=DenseSystem(system.grid,physics)
    enable_species!(tfsystem,1,[i for i=1:num_cellregions(system.grid)])
    return TestFunctionFactory{Tv}(system,tfsystem)
end



############################################################################
"""
$(SIGNATURES)

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

    _eval_and_assemble(factory.tfsystem,u,u,f,Inf,
                       factory.tfsystem.physics.storage,
                       factory.tfsystem.physics.bstorage,
                       factory.tfsystem.physics.source)

    _initialize!(u,factory.tfsystem)
    lufact=LinearAlgebra.lu(factory.tfsystem.matrix)
    ldiv!(vec(u),lufact,vec(f))
    return vec(u)
end


############################################################################
"""
$(SIGNATURES)

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
    node=Node{Tv}(this)
    edge=Edge{Tv}(this)
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
            _fill!(edge,grid,iedge,icell)

            for ispec=1:nspecies
                UKL[ispec]=U[ispec,edge.nodeK]
                UKL[ispec+nspecies]=U[ispec,edge.nodeL]
            end
            @views this.physics.flux(res,UKL,edge, this.physics.data)
            for ispec=1:nspecies
                if this.node_dof[ispec,edge.nodeK]==ispec && this.node_dof[ispec,edge.nodeL]==ispec
                    integral[ispec]+=edge_factors[iedge]*res[ispec]*(tf[edge.nodeK]-tf[edge.nodeL])
                end
            end
        end
        
        for inode=1:num_nodes_per_cell(grid)
            _fill!(node,grid,inode,icell)
            @views begin
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

############################################################################
"""
$(SIGNATURES)

Calculate test function integral.
"""
integrate(this::AbstractSystem{Tv},tf::Vector{Tv},U::AbstractMatrix{Tv}) where Tv=integrate(this,tf,U,U,Inf)



