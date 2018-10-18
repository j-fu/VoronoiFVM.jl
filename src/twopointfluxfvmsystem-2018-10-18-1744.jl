# Packages for Autodiff magic. These need to be installed via Pkg
using ForwardDiff, DiffResults
using IterativeSolvers

# These are in the standard distro
using SparseArrays
using LinearAlgebra
using Printf


const Dirichlet=1.0e30

"""
    Main structure holding data for system solution
"""
mutable struct TwoPointFluxFVMSystem
    geometry::FVMGraph
    parameters::FVMParameters
    num_bnodes::Array{Int64,1}
    num_dof::Int64
    num_bulk_dof::Int64
    num_bregion_nodes::Array{Int64,1}
    num_bspecies::Array{Int64,1}
    bdof_offset::Array{Int64,1}
    bnode_index::Array{Int64,1}
    boundary_values::Array{Float64,2}
    boundary_factors::Array{Float64,2}
    matrix::SparseArrays.SparseMatrixCSC
    residual::Array{Float64,1}
    update::Array{Float64,1}
    TwoPointFluxFVMSystem(geometry::FVMGraph,parameters::FVMParameters)=TwoPointFluxFVMSystem(new(),geometry,parameters)
end


function TwoPointFluxFVMSystem(this::TwoPointFluxFVMSystem,
                               geometry::FVMGraph, # Geometry
                               parameters::FVMParameters# user parameter
                               )
    this.geometry=geometry
    this.parameters=parameters
    
    this.num_bregion_nodes=zeros(Int64,geometry.num_bregions)
    this.bnode_index=zeros(Int64,length(geometry.bnode_nodes))

    for i=1:length(geometry.bnode_nodes)
        ireg=geometry.bnode_regions[i]
        this.num_bregion_nodes[ireg]+=1
        this.bnode_index[i]=this.num_bregion_nodes[ireg]
    end
    this.num_bspecies=zeros(Int64,geometry.num_bregions)

    for ibreg=1:min(geometry.num_bregions,length(parameters.num_bspecies))
        this.num_bspecies[ibreg]=parameters.num_bspecies[ibreg]
    end
    


    # Arrays for boundary data
    this.boundary_values=zeros(parameters.num_species,geometry.num_bregions)
    this.boundary_factors=zeros(parameters.num_species,geometry.num_bregions)
    
    this.num_dof=parameters.num_species*geometry.num_nodes
    this.num_bulk_dof=this.num_dof

    this.bdof_offset=zeros(Int64,geometry.num_bregions+1)
    for ibreg=1:geometry.num_bregions
        this.bdof_offset[ibreg]=this.num_dof+1
        this.num_dof+=this.num_bspecies[ibreg]*this.num_bregion_nodes[ibreg]
    end
    this.bdof_offset[geometry.num_bregions+1]=this.num_dof+1
    print(this.bdof_offset)

    # Empty sparse matrix
    this.matrix=SparseArrays.spzeros(this.num_dof,this.num_dof) # Jacobi matrix
    
    # Iteration data
    # These are created as 1D arrays because they must fit to sparse matrix
    this.residual=Array{Float64,1}(undef,this.num_dof)
    this.update=Array{Float64,1}(undef,this.num_dof)
    return this
end


"""
    function unknowns(this::TwoPointFluxFVMSystem)

Create a vector of unknowns for a given system
"""
function unknowns(this::TwoPointFluxFVMSystem)
    return Array{Float64,1}(undef,this.num_dof)
end

function bulk_unknowns(this::TwoPointFluxFVMSystem,U::Array{Float64,1})
    V=view(U,1:this.num_bulk_dof)
    return reshape(V,this.parameters.num_species,this.geometry.num_nodes)
end

function boundary_unknowns(this::TwoPointFluxFVMSystem, U::Array{Float64,1}, ibc::Int)
    V=view(U,this.bdof_offset[ibc]:this.bdof_offset[ibc+1]-1)
    return reshape(V,this.num_bspecies[ibc],this.bdof_offset[ibc+1]-this.bdof_offset[ibc])
end


"""
    Initialize dirichlet boundary valuesfor solution
"""
function inidirichlet!(this::TwoPointFluxFVMSystem,U0)
    U=bulk_unknowns(this,U0)
    geometry=this.geometry
    nbnodes=length(geometry.bnode_nodes)
    for ibnode=1:nbnodes
        ibreg=geometry.bnode_regions[ibnode]
        inode=geometry.bnode_nodes[ibnode]
        for ispec=1:this.parameters.num_species
            if this.boundary_factors[ispec,ibreg]==Dirichlet
                U[ispec,inode]=this.boundary_values[ispec,ibreg]
            end
        end
    end
end

"""
    function integrate(this::TwoPointFluxFVMSystem,F::Function,U)

Integrate solution vector over domain
"""
function integrate(this::TwoPointFluxFVMSystem,F::Function,U0)
    U=bulk_unknowns(this,U0)
    nnodes=this.geometry.num_nodes
    nspec=this.parameters.num_species
    nodefac=this.geometry.node_factors
    integral=zeros(nspec)
    res=zeros(nspec)
    for inode=1:nnodes
        F(this.parameters,res,U[:,inode])
        for ispec=1:nspec
            integral[ispec]+=nodefac[inode]*res[ispec]
        end
    end
    return integral
end



"""
  Nonlinear operator evaluation + Jacobian assembly

  for K=1...n
     f_K = sum_(L neigbor of K) flux(u[K],U[L])*edgefac[K,L]
            + (reaction(U[K])- source(X[K]))*nodefac[K]
            + (storage(U[K])- storage(UOld[K])*nodefac[K]/tstep

    # M is correspondig Jacobi matrix of derivatives. 


"""
function eval_and_assemble(this::TwoPointFluxFVMSystem,
                           U0, # Actual solution iteration
                           UOld0, # Old timestep solution
                           tstep # time step size. Inf means stationary solution
                           )

    """
       Wrap API flux with function compatible to ForwardDiff
    """
    function fluxwrap!(y,u)
        flux!(y,u[1:num_species],u[num_species+1:2*num_species])
    end

    function breawrap!(y,u)
        iy=view(y,1:num_species)
        by=view(y,num_species+1:num_species+num_bspecies)
        breaction!(y,by,u[1:num_species],u[num_species+1:num_species+num_bspecies])
    end
    

    U=bulk_unknowns(this,U0)
    UOld=bulk_unknowns(this,UOld0)


    parameters=this.parameters
    # Create closures for physics functions
    # These allow to "glue" user parameters to function objects compatible
    # with the ForwardDiff module
    source!(y,x)=parameters.source(parameters,y,x)
    flux!(y,uk,ul)=parameters.flux(parameters,y,uk,ul)
    reaction!(y,x)=parameters.reaction(parameters,y,x)
    breaction!(y,by,u,bu)=parameters.breaction(parameters,y,by,u,bu)
    storage!(y,x)=parameters.storage(parameters,y,x)
    


    geometry=this.geometry
    nnodes=geometry.num_nodes
    num_species=this.parameters.num_species
    nedges=size(geometry.edge_nodes,2)
    M=this.matrix
    F=bulk_unknowns(this,this.residual)
    
    # Reset matrix + rhs
    M.nzval.=0.0
    F.=0.0

    # Assemble nonlinear term + source + storage using autodifferencing via ForwardDiff

    # struct holding diff results for storage, reaction
    result_r=DiffResults.DiffResult(Vector{Float64}(undef,num_species),Matrix{Float64}(undef,num_species,num_species))
    result_s=DiffResults.DiffResult(Vector{Float64}(undef,num_species),Matrix{Float64}(undef,num_species,num_species))

    

    bunknowns=[boundary_unknowns(this,U0,ibc) for ibc=1:geometry.num_bregions]
    result_b=[DiffResults.DiffResult(Vector{Float64}(undef,num_species+this.num_bspecies[ibc]),
                                     Matrix{Float64}(undef,num_species+this.num_bspecies[ibc],num_species+this.num_bspecies[ibc]))
              for ibc=1:geometry.num_bregions]
    
    BU=[Array{Float64,1}(undef,num_species+this.num_bspecies[ibc]) for ibc=1:geometry.num_bregions]
    BY=[Array{Float64,1}(undef,num_species+this.num_bspecies[ibc]) for ibc=1:geometry.num_bregions]

    # array providing space for function arguments
    Y=Array{Float64}(undef,num_species)
    
    # array holding source term
    src=Array{Float64}(undef,num_species)

    # array holding storage term for old solution
    oldstor=Array{Float64}(undef,num_species)

    iblock=0 # block offset

    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here.
    tstepinv=1.0/tstep 
    
    for inode=1:nnodes
        # Evaluate & differentiate reaction term
        result_r=ForwardDiff.jacobian!(result_r,reaction!,Y,U[:,inode])
        res_react=DiffResults.value(result_r)
        jac_react=DiffResults.jacobian(result_r)
        # Evaluate source term
        source!(src,geometry.node_coordinates[:,inode])

        # Evaluate & differentiate storage term
        result_s=ForwardDiff.jacobian!(result_s,storage!,Y,U[:,inode])
        res_stor=DiffResults.value(result_s)
        jac_stor=DiffResults.jacobian(result_s)

        # Evaluate storage term for old timestep
        storage!(oldstor,UOld[:,inode])

        # Assembly results and jacobians
        for i=1:num_species
            F[i,inode]+=geometry.node_factors[inode]*(res_react[i]-src[i] + (res_stor[i]-oldstor[i])*tstepinv)
            for j=1:num_species
                M[iblock+i,iblock+j]+=geometry.node_factors[inode]*(jac_react[i,j]+ jac_stor[i,j]*tstepinv)
            end
        end
        iblock+=num_species
    end
    
    # Create result struct for flux evaluation
    result=DiffResults.DiffResult(Vector{Float64}(undef,num_species),Matrix{Float64}(undef,num_species,2*num_species))
    Y=Array{Float64,1}(undef,num_species)
    UKL=Array{Float64,1}(undef,2*num_species)
    # Assemble main part


    for iedge=1:nedges
        K=geometry.edge_nodes[1,iedge]
        L=geometry.edge_nodes[2,iedge]
        # Set up argument for fluxwrap!
        UKL[1:num_species]=U[:,K]
        UKL[num_species+1:2*num_species]=U[:,L]
        result=ForwardDiff.jacobian!(result,fluxwrap!,Y,UKL)

        res=DiffResults.value(result)
        jac=DiffResults.jacobian(result)

        # Assemble flux data
        F[:,K]+=res*geometry.edge_factors[iedge]
        F[:,L]-=res*geometry.edge_factors[iedge]
        kblock=(K-1)*num_species
        lblock=(L-1)*num_species
        jl=num_species+1
        for jk=1:num_species
            for ik=1:num_species
                M[kblock+ik,kblock+jk]+=jac[ik,jk]*geometry.edge_factors[iedge]
                M[kblock+ik,lblock+jk]+=jac[ik,jl]*geometry.edge_factors[iedge]
                M[lblock+ik,kblock+jk]-=jac[ik,jk]*geometry.edge_factors[iedge]
                M[lblock+ik,lblock+jk]-=jac[ik,jl]*geometry.edge_factors[iedge]
            end
            jl+=1
        end
    end
    
    # Assemble boundary conditions
    # Dirichlet conditions are handeld via penalty method, 
    # for these, solution vectors need to be initialized
    # appropriately
    nbnodes=length(geometry.bnode_nodes)
    for ibnode=1:nbnodes
        inode=geometry.bnode_nodes[ibnode]
        ibreg=geometry.bnode_regions[ibnode]
        iblock=(inode-1)*num_species
        for ispec=1:num_species
            fac=this.boundary_factors[ispec,ibreg]
            if fac!=Dirichlet
                fac*=geometry.bnode_factors[ibnode]
            end
            F[ispec,inode]+=fac*(U[ispec,inode]-this.boundary_values[ispec,ibreg])
            M[iblock+ispec,iblock+ispec]+=fac
        end

        BU[ibreg][1:num_species]=U[:,inode]
        BU[ibreg][num_species+1:num_species+parameters.num_bspecies[ibreg]]=bunknowns[ibreg][1:parameters.num_bspecies[ibreg]]
        this.parameters.bregion=ibreg
        num_species=this.num_bspecies[ibreg]
        result_b[ibreg]=ForwardDiff.jacobian!(result_b[ibreg],breawrap!,BY[ibreg],BU[ibreg])
        res_breact=DiffResults.value(result_b[ibreg])
        jac_breact=DiffResults.jacobian(result_b[ibreg])
        # Assembly results and jacobians
        for i=1:num_species
            F[i,inode]+=geometry.bnode_factors[ibnode]*(res_breact[i])
            for j=1:num_species
                M[iblock+i,iblock+j]+=geometry.bnode_factors[ibnode]*(jac_breact[i,j])
            end
        end
    end
end



"""
    Actual solver function implementation
"""
function _solve(this::TwoPointFluxFVMSystem, oldsol::Array{Float64,1},control::FVMNewtonControl, tstep::Float64)
    solution=copy(oldsol)
    inidirichlet!(this,solution)

    # Newton iteration (quick and dirty...)
    oldnorm=1.0
    converged=false
    if control.verbose
        @printf("Start newton iteration: %s:%d\n", basename(@__FILE__),@__LINE__)
    end
    nlu=0
    lufact=nothing
    damp=control.damp_initial
    for ii=1:control.max_iterations
        eval_and_assemble(this,solution,oldsol,tstep)
        
        # Sparse LU factorization
        # Here, we seem miss the possibility to re-use the 
        # previous symbolic information
        # We hower reuse the factorization control.max_lureuse times.
        if nlu==0
            lufact=LinearAlgebra.lu(this.matrix)
            # LU triangular solve gives Newton update
            ldiv!(this.update,lufact,this.residual)
        else
            # When reusing lu factorization, we may try to iterate
            # Generally, this is advisable.
            if control.tol_linear <1.0
                bicgstabl!(this.update,this.matrix,this.residual,2,Pl=lufact,tol=control.tol_linear)
            else
                ldiv!(this.update,lufact,this.residual)
            end
        end
        nlu=min(nlu+1,control.max_lureuse)

        # vector expressions would allocate here...
        for i=1:this.num_dof
            solution[i]-=damp*this.update[i]
        end
        damp=min(damp*control.damp_growth,1.0)
        norm=LinearAlgebra.norm(this.update,Inf)/this.num_dof
        if control.verbose
            @printf("  it=%03d norm=%.5e cont=%.5e\n",ii,norm, norm/oldnorm)
        end
        if norm<control.tol_absolute
            converged=true
            break
        end
        
        oldnorm=norm
    end
    if !converged
        error("Error: no convergence")
    end
    return solution
end

"""
    System solver wrapper allowing to dispatch timing
"""
function solve(this::TwoPointFluxFVMSystem, # Finite volume system
               oldsol::Array{Float64,1}; # old time step solution resp. initial value
               control=FVMNewtonControl(), # Newton solver control information
               tstep::Float64=Inf) # Time step size. Inf means  stationary solution
    if control.verbose
        @time begin
            retval= _solve(this,oldsol,control,tstep)
        end
        return retval
    else
        return _solve(this,oldsol,control,tstep)
    end
end
