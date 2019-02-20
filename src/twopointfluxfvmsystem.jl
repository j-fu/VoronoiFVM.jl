# Packages for Autodiff magic. These need to be installed via Pkg
using ForwardDiff, DiffResults
using IterativeSolvers

# These are in the standard distro
using SparseArrays
using LinearAlgebra
using Printf


#####################################################
"""
Constant to be used as boundary condition factor 
to mark Dirichlet boundary conditons.    
"""
const Dirichlet=1.0e30

#####################################################
"""
````
mutable struct System
````
Main structure holding data for system solution

Public fields:
    
    boundary_values::Array{Float64,2} # Array of boundary values
    boundary_factors::Array{Float64,2} # Array of boundary factors
    geometry::Graph   # Geometry information: weighted graph created from grid
    physics::Physics  # Physical model 


Memory layout for solution arrays:

 bulk_dof | bregion_dof[1] |bregion_dof[2] |...


"""
mutable struct System
    boundary_values::Array{Float64,2} # Array of boundary values
    boundary_factors::Array{Float64,2} # Array of boundary factors

    geometry::Graph   # Geometry information: weighted graph created from grid
    physics::Physics  # Physical model 
    _num_dof::Int64       # Overall number of degrees of freedom
    _num_bulk_dof::Int64  # Number of degrees of freedom in the bulk
    _num_bspecies::Array{Int64,1} # Number of boundary species per boundary region
    _bdof_offset::Array{Int64,1}  # offset of boundary degrees of freedom per boundary region
    _matrix::SparseArrays.SparseMatrixCSC # System matrix
    _residual::Array{Float64,1} # Residual array
    _update::Array{Float64,1}   # Newton update array

    System(geometry::Graph,physics::Physics)=System(new(),geometry,physics)
end


#####################################################
"""
Constructor for System

````
function System(this::System,
                               geometry::Graph, # Geometry
                               physics::Physics # user data
                               )
````

Construct system from geometry (graph) and user data.

"""
function System(this::System,
                geometry::Graph, # Geometry
                physics::Physics # user data
                )
    this.geometry=geometry
    this.physics=physics
    

    this._num_bspecies=zeros(Int64,geometry.num_bregions)

    for ibreg=1:min(geometry.num_bregions,length(physics.num_bspecies))
        this._num_bspecies[ibreg]=physics.num_bspecies[ibreg]
    end
    


    # Arrays for boundary data
    this.boundary_values=zeros(physics.num_species,geometry.num_bregions)
    this.boundary_factors=zeros(physics.num_species,geometry.num_bregions)
    
    this._num_dof=physics.num_species*geometry.num_nodes
    this._num_bulk_dof=this._num_dof

    this._bdof_offset=zeros(Int64,geometry.num_bregions+1)
    for ibreg=1:geometry.num_bregions
        this._bdof_offset[ibreg]=this._num_dof
        this._num_dof+=this._num_bspecies[ibreg]*geometry.num_bregion_nodes[ibreg]
    end
    this._bdof_offset[geometry.num_bregions+1]=this._num_dof

    # Empty sparse matrix
    this._matrix=SparseArrays.spzeros(this._num_dof,this._num_dof) # Jacobi matrix
    
    # Iteration data
    # These are created as 1D arrays because they must fit to sparse matrix
    this._residual=Array{Float64,1}(undef,this._num_dof)
    this._update=Array{Float64,1}(undef,this._num_dof)
    return this
end


#####################################################
"""
````
function unknowns(this::System)
````

Create a vector of unknowns for a given system
"""
function unknowns(this::System)
    return Array{Float64,1}(undef,this._num_dof)
end

#####################################################
"""
````
function bulk_unknowns(this::System,U::Array{Float64,1})
````

Return a view  which allows to access the bulk unknowns stored
in a solution vector.
"""
function bulk_unknowns(this::System,U::Array{Float64,1})
    V=view(U,1:this._num_bulk_dof)
    return reshape(V,this.physics.num_species,this.geometry.num_nodes)
end


#####################################################
"""
````
function boundary_unknowns(this::System, U::Array{Float64,1}, ibc::Int)
````

Return a view  which allows to access the unknwons in a solution vector
which correspond to a certain boundary condition.
"""
function boundary_unknowns(this::System, U::Array{Float64,1}, ibc::Int)
    V=view(U,this._bdof_offset[ibc]+1:this._bdof_offset[ibc+1])
    return reshape(V,this._num_bspecies[ibc],this._bdof_offset[ibc+1]-this._bdof_offset[ibc])
end


#####################################################
"""
````
function inidirichlet!(this::System,U0)
````

Initialize dirichlet boundary values for solution
"""
function inidirichlet!(this::System,U0)
    U=bulk_unknowns(this,U0)
    geometry=this.geometry
    nbnodes=length(geometry.bnode_nodes)
    for ibnode=1:nbnodes
        ibreg=geometry.bnode_regions[ibnode]
        inode=geometry.bnode_nodes[ibnode]
        for ispec=1:this.physics.num_species
            if this.boundary_factors[ispec,ibreg]==Dirichlet
                U[ispec,inode]=this.boundary_values[ispec,ibreg]
            end
        end
    end
end

##########################################################
"""
````
function integrate(this::System,F::Function,U)
````

Integrate solution vector over domain. Returns an `Array{Int64,1}`
containing the integral for each species.
"""
function integrate(this::System,F::Function,U0)
    U=bulk_unknowns(this,U0)
    nnodes=this.geometry.num_nodes
    nspec=this.physics.num_species
    nodefac=this.geometry.node_factors
    integral=zeros(nspec)
    res=zeros(nspec)
    node=Node()
    for inode=1:nnodes
        node.index=inode
        node.coord=view(this.geometry.node_coordinates,:,inode)
        F(this.physics,node,res,U[:,inode])
        for ispec=1:nspec
            integral[ispec]+=nodefac[inode]*res[ispec]
        end
    end
    return integral
end



#####################################################
"""
Nonlinear operator evaluation + Jacobian assembly

````
  for K=1...n
     f_K = sum_(L neigbor of K) flux(u[K],U[L])*edgefac[K,L]
            + (reaction(U[K])- source(X[K]))*nodefac[K]
            + (storage(U[K])- storage(UOld[K])*nodefac[K]/tstep

    # M is correspondig Jacobi matrix of derivatives. 
````

"""
function eval_and_assemble(this::System,
                           U0, # Actual solution iteration
                           UOld0, # Old timestep solution
                           tstep # time step size. Inf means stationary solution
                           )

    """
    Wrap API flux with function compatible to ForwardDiff
    """
    function fluxwrap!(y,u)
        uk=view(u,1:num_species)
        ul=view(u,num_species+1:2*num_species)
        flux!(y,uk,ul)
    end

    """
    Wrap API boundary reaction with function compatible to ForwardDiff
    """
    function breawrap!(y,u)
        iy=view(y,1:num_species)
        by=view(y,num_species+1:num_species+num_bspecies)
        iu=view(u,1:num_species)
        bu=view(u,num_species+1:num_species+num_bspecies)
        breaction!(iy,by,iu,bu)
    end
    

    U=bulk_unknowns(this,U0)
    UOld=bulk_unknowns(this,UOld0)

    physics=this.physics
    node=Node()
    edge=Edge()
    # Create closures for physics functions
    # These allow to "glue" user physics to function objects compatible
    # with the ForwardDiff module
    # cf. http://www.juliadiff.org/ForwardDiff.jl/stable/user/limitations.html 
    source!(y)=physics.source(physics,node,y)
    flux!(y,uk,ul)=physics.flux(physics,edge,y,uk,ul)
    reaction!(y,x)=physics.reaction(physics,node,y,x)
    storage!(y,x)=physics.storage(physics,node,y,x)
    


    geometry=this.geometry
    nnodes=geometry.num_nodes
    num_species=this.physics.num_species
    num_bspecies=0
    nedges=size(geometry.edge_nodes,2)
    M=this._matrix
    F=bulk_unknowns(this,this._residual)
    
    # Reset matrix + rhs
    M.nzval.=0.0
    this._residual.=0.0

    # Assemble nonlinear term + source + storage using autodifferencing via ForwardDiff

    # struct holding diff results for storage, reaction
    result_r=DiffResults.DiffResult(Vector{Float64}(undef,num_species),Matrix{Float64}(undef,num_species,num_species))
    result_s=DiffResults.DiffResult(Vector{Float64}(undef,num_species),Matrix{Float64}(undef,num_species,num_species))

    
    if this.physics.breaction != prototype_breaction! ||  this.physics.bstorage != prototype_bstorage! 
        breaction!(y,by,u,bu)=physics.breaction(physics,node,y,by,u,bu)
        bstorage!(by,bu)=physics.bstorage(physics,node,by,bu)

        result_br=[DiffResults.DiffResult(Vector{Float64}(undef,num_species+this._num_bspecies[ibc]),
                                          Matrix{Float64}(undef,num_species+this._num_bspecies[ibc],num_species+this._num_bspecies[ibc]))
                   for ibc=1:geometry.num_bregions]
        
        result_bs=[DiffResults.DiffResult(Vector{Float64}(undef,this._num_bspecies[ibc]),
                                          Matrix{Float64}(undef,this._num_bspecies[ibc],this._num_bspecies[ibc]))
                   for ibc=1:geometry.num_bregions]
        
        BU=[Array{Float64,1}(undef,num_species+this._num_bspecies[ibc]) for ibc=1:geometry.num_bregions]
        BY=[Array{Float64,1}(undef,num_species+this._num_bspecies[ibc]) for ibc=1:geometry.num_bregions]
        
        BUS=[Array{Float64,1}(undef,this._num_bspecies[ibc]) for ibc=1:geometry.num_bregions]
        BYS=[Array{Float64,1}(undef,this._num_bspecies[ibc]) for ibc=1:geometry.num_bregions]
    end

    # array providing space for function arguments
    Y=Array{Float64}(undef,num_species)
    
    # array holding source term
    src=zeros(Float64,num_species)

    # array holding storage term for old solution
    oldstor=zeros(Float64,num_species)

    iblock=0 # block offset

    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here.
    tstepinv=1.0/tstep 

    res_react=zeros(Float64,num_species)
    jac_react=zeros(Float64,num_species,num_species)

    res_stor=zeros(Float64,num_species)
    jac_stor=zeros(Float64,num_species,num_species)

    for inode=1:nnodes
        node.index=inode
        node.coord=view(geometry.node_coordinates,:,inode)
        # Evaluate & differentiate reaction term
        if physics.reaction!=prototype_reaction!
            result_r=ForwardDiff.jacobian!(result_r,reaction!,Y,U[:,inode])
            res_react=DiffResults.value(result_r)
            jac_react=DiffResults.jacobian(result_r)
        end

        # Evaluate source term
        if physics.source!=prototype_source!
            source!(src)
        end
        
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
        edge.index=iedge
        edge.nodeK=K
        edge.nodeL=L
        edge.coordL=view(geometry.node_coordinates,:,L)
        edge.coordK=view(geometry.node_coordinates,:,K)

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
        node.index=ibnode
        node.coord=view(geometry.node_coordinates,:,ibnode)

        # Standard boundary conditions
        #
        inode=geometry.bnode_nodes[ibnode]
        ibreg=geometry.bnode_regions[ibnode]
        iblock=(inode-1)*num_species
        for ispec=1:num_species
            fac=this.boundary_factors[ispec,ibreg]
            if fac>0.0
                if fac!=Dirichlet
                    fac*=geometry.bnode_factors[ibnode]
                end
                F[ispec,inode]+=fac*(U[ispec,inode]-this.boundary_values[ispec,ibreg])
                M[iblock+ispec,iblock+ispec]+=fac
            end
        end
        
        if this.physics.breaction != prototype_breaction! 
            ibxnode=0
            ibxblock=0
            num_bspecies=this._num_bspecies[ibreg]

            BU[ibreg][1:num_species]=U[:,inode]
            if num_bspecies>0
                ibxnode=geometry.bnode_index[ibnode]
                ibxblock=this._bdof_offset[ibreg]+(ibxnode-1)*num_bspecies
                BU[ibreg][num_species+1:num_species+num_bspecies]=U0[ibxblock+1:ibxblock+num_bspecies]
            end
            
            this.physics.bregion=ibreg
            result_br[ibreg]=ForwardDiff.jacobian!(result_br[ibreg],breawrap!,BY[ibreg],BU[ibreg])
            res_breact=DiffResults.value(result_br[ibreg])
            jac_breact=DiffResults.jacobian(result_br[ibreg])
            
            # Assembly results and jacobians
            for i=1:num_species
                F[i,inode]+=geometry.bnode_factors[ibnode]*(res_breact[i])
                for j=1:num_species
                    M[iblock+i,iblock+j]+=geometry.bnode_factors[ibnode]*(jac_breact[i,j])
                end
            end
            
            if num_bspecies>0
                for i=1:num_bspecies
                    this._residual[ibxblock+i]+=res_breact[num_species+i]
                    for j=1:num_bspecies
                        M[ibxblock+i,ibxblock+j]+=jac_breact[num_species+i,num_species+j]
                    end
                    for j=1:num_species
                        M[ibxblock+i,iblock+j]+=jac_breact[num_species+i,j]
                        M[iblock+j,ibxblock+i]+=geometry.bnode_factors[ibnode]*jac_breact[j,num_species+i]
                    end
                end
            end
        end
        
        if this.physics.bstorage != prototype_bstorage! 
            ibxnode=0
            ibxblock=0
            num_bspecies=this._num_bspecies[ibreg]
            if num_bspecies>0
                ibxnode=geometry.bnode_index[ibnode]
                ibxblock=this._bdof_offset[ibreg]+(ibxnode-1)*num_bspecies
                bu=view(U0,ibxblock+1:ibxblock+num_bspecies)
                

                this.physics.bregion=ibreg
                result_bs[ibreg]=ForwardDiff.jacobian!(result_bs[ibreg],bstorage!,BYS[ibreg],bu)
                res_bstor=DiffResults.value(result_bs[ibreg])
                jac_bstor=DiffResults.jacobian(result_bs[ibreg])
            
                buold=view(UOld0,ibxblock+1:ibxblock+num_bspecies)
                # Evaluate storage term for old timestep
                bstorage!(BYS[ibreg],buold)

                for i=1:num_bspecies
                    this._residual[ibxblock+i]+=(res_bstor[i]-BYS[ibreg][i])*tstepinv
                    for j=1:num_bspecies
                        M[ibxblock+i,ibxblock+j]+=jac_bstor[i,j]*tstepinv
                    end
                end
            end
        end
    end
end


################################################################
"""
Actual solver function implementation
"""
function _solve(this::System, oldsol::Array{Float64,1},control::NewtonControl, tstep::Float64)
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
    tolx=0.0
    for ii=1:control.max_iterations
        eval_and_assemble(this,solution,oldsol,tstep)
        
        # Sparse LU factorization
        # Here, we seem miss the possibility to re-use the 
        # previous symbolic information
        # We however reuse the factorization control.max_lureuse times.
        if nlu==0
            lufact=LinearAlgebra.lu(this._matrix)
            # LU triangular solve gives Newton update
            ldiv!(this._update,lufact,this._residual)
        else
            # When reusing lu factorization, we may try to iterate
            # Generally, this is advisable.
            if control.tol_linear <1.0
                bicgstabl!(this._update,this._matrix,this._residual,2,Pl=lufact,tol=control.tol_linear)
            else
                ldiv!(this._update,lufact,this._residual)
            end
        end
        nlu=min(nlu+1,control.max_lureuse)
        # vector expressions would allocate here...
        for i in eachindex(solution)
            solution[i]-=damp*this._update[i]
        end
        damp=min(damp*control.damp_growth,1.0)
        norm=LinearAlgebra.norm(this._update,Inf)/this._num_dof
        if tolx==0.0
            tolx=norm*control.tol_relative
        end
        if control.verbose
            @printf("  it=%03d norm=%.5e cont=%.5e\n",ii,norm, norm/oldnorm)
        end
        if norm<control.tol_absolute || norm <tolx
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

################################################################
"""
Solution method for instance of System

````
function solve(
    this::System, # Finite volume system
    oldsol::Array{Float64,1};    # old time step solution resp. initial value
    control=NewtonControl(),  # Solver control information
    tstep::Float64=Inf           # Time step size. Inf means  stationary solution
    )
````
Perform solution of stationary system (if `tstep==Inf`) or implicit Euler time
step system. 

"""
function solve(
    this::System, # Finite volume system
    oldsol::Array{Float64,1}; # old time step solution resp. initial value
    control=NewtonControl(), # Newton solver control information
    tstep::Float64=Inf          # Time step size. Inf means  stationary solution
    )

    if control.verbose
        @time begin
            retval= _solve(this,oldsol,control,tstep)
        end
        return retval
    else
        return _solve(this,oldsol,control,tstep)
    end
end


