# Packages for Autodiff magic. These need to be installed via Pkg
using ForwardDiff, DiffResults

# These are in the standard distro
using SparseArrays
using LinearAlgebra
using Printf


"""
   struct containing control data for the nonlinar solver
"""
mutable struct FVMNewtonControl
    tolerance::Float64 # Tolerance (in terms of norm of Newton update)
    damp::Float64      # Initial damping parameter
    maxiter::Int32     # Maximum number of iterations
    verbose::Bool      # Verbosity
    function FVMNewtonControl()
        new(1.0e-10,1.0,100,true)
    end
end


"""
   Abstract type for user problem data.
 
   Derived types must contain field number_of_species::Int64
"""
abstract type FVMParameters end


"""
   Problem data type for default parameter function
"""
mutable struct DefaultParameters <: FVMParameters
    number_of_species::Int64
end

"""
    Default source term
"""
function default_source!(this::FVMParameters, f,x)
    for i=1:this.number_of_species
        f[i]=0
    end
end

"""
    Default reaction term
"""
function default_reaction!(this::FVMParameters, f,u)
    for i=1:this.number_of_species
        f[i]=0
    end
end

"""
    Default flux term
"""
function default_flux!(this::FVMParameters, f,uk,ul)
    for i=1:this.number_of_species
        f[i]=uk[i]-ul[i]
    end
end



"""
    Default storage term
"""
function default_storage!(this::FVMParameters, f,u)
    for i=1:this.number_of_species
        f[i]=u[i]
    end
end


const Dirichlet=1.0e30

"""
    Main structure holding data for system solution
"""
struct TwoPointFluxFVMSystem
    geometry::FVMGraph
    number_of_species::Int64
    source!::Function
    reaction!::Function
    storage!::Function
    flux!::Function
    boundary_values::Array{Float64,2}
    boundary_factors::Array{Float64,2}
    matrix::SparseArrays.SparseMatrixCSC
    residual::Array{Float64,1}
    update::Array{Float64,1}
    function TwoPointFluxFVMSystem(geometry::FVMGraph; # Geometry
                                   parameters::FVMParameters=DefaultParameters(1), # user parameters
                                   source::Function=default_source!,
                                   reaction::Function=default_reaction!,
                                   storage::Function=default_storage!,
                                   flux::Function=default_flux!)

        number_of_species=parameters.number_of_species

        # Create closures for physics functions
        # These allow to "glue" user parameters to function objects compatible
        # with the ForwardDiff module
        _source!(y,x)=source(parameters,y,x)
        _flux!(y,uk,ul)=flux(parameters,y,uk,ul)
        _reaction!(y,x)=reaction(parameters,y,x)
        _storage!(y,x)=storage(parameters,y,x)

        # Arrays for boundary data
        boundary_values=zeros(number_of_species,geometry.NumberOfBoundaryRegions)
        boundary_factors=zeros(number_of_species,geometry.NumberOfBoundaryRegions)



        # Empty sparse matrix
        matrix=SparseArrays.spzeros(geometry.NumberOfNodes*number_of_species,geometry.NumberOfNodes*number_of_species) # Jacobi matrix

        # Iteration data
        # These are created as 1D arrays because they must fit to sparse matrix
        residual=Array{Float64,1}(undef,geometry.NumberOfNodes*number_of_species)
        update=Array{Float64,1}(undef,geometry.NumberOfNodes*number_of_species)


        new(geometry,
            number_of_species,
            _source!,
            _reaction!,
            _storage!,
            _flux!,
            boundary_values,
            boundary_factors,
            matrix,
            residual,
            update)
    end
end

"""
    Create a vector of unknowns for a given system
"""
function unknowns(fvsystem::TwoPointFluxFVMSystem)
    return Array{Float64,2}(undef,fvsystem.number_of_species,fvsystem.geometry.NumberOfNodes)
end


"""
    Initialize dirichlet boundary valuesfor solution
"""
function inidirichlet!(fvsystem::TwoPointFluxFVMSystem,U)
    geom=fvsystem.geometry
    nbnodes=length(geom.BoundaryNodes)
    for ibnode=1:nbnodes
        ibreg=geom.BoundaryRegions[ibnode]
        for ispec=1:fvsystem.number_of_species
            if fvsystem.boundary_factors[ispec,ibreg]==Dirichlet
                U[ispec,ibnode]=fvsystem.boundary_values[ispec,ibreg]
            end
        end
    end
end


"""
  Nonlinear operator evaluation + Jacobian assembly

  for K=1...n
     f_K = sum_(L neigbor of K) flux(u[K],U[L])*edgefac[K,L]
            + (reaction(U[K])- source(X[K]))*nodefac[K]
            + (storage(U[K])- storage(UOld[K])*nodefac[K]/tstep

    # M is correspondig Jacobi matrix of derivatives. 


"""
function eval_and_assemble(fvsystem::TwoPointFluxFVMSystem,
                           U, # Actual solution iteration
                           UOld, # Old timestep solution
                           tstep # time step size. Inf means stationary solution
                           )

    """
       Wrap API flux with function compatible to ForwardDiff
    """
    function fluxwrap!(y,u)
        fvsystem.flux!(y,u[1:number_of_species],u[number_of_species+1:2*number_of_species])
    end
    
    geom=fvsystem.geometry
    nnodes=geom.NumberOfNodes
    number_of_species=fvsystem.number_of_species
    nedges=size(geom.Edges,2)
    M=fvsystem.matrix
    F=reshape(fvsystem.residual,number_of_species,nnodes)
    
    # Reset matrix + rhs
    M.nzval.=0.0
    F.=0.0

    # Assemble nonlinear term + source + storage using autodifferencing via ForwardDiff

    # struct holding diff results for storage, reaction
    result_r=DiffResults.DiffResult(Vector{Float64}(undef,number_of_species),Matrix{Float64}(undef,number_of_species,number_of_species))
    result_s=DiffResults.DiffResult(Vector{Float64}(undef,number_of_species),Matrix{Float64}(undef,number_of_species,number_of_species))

    # array providing space for function arguments
    Y=Array{Float64}(undef,number_of_species)
    
    # array holding source term
    src=Array{Float64}(undef,number_of_species)

    # array holding storage term for old solution
    oldstor=Array{Float64}(undef,number_of_species)

    iblock=0 # block offset

    # Inverse of timestep
    # According to Julia documentation, 1/Inf=0 which
    # comes handy to write compact code here.
    tstepinv=1.0/tstep 
    
    for inode=1:nnodes
        # Evaluate & differentiate reaction term
        result_r=ForwardDiff.jacobian!(result_r,fvsystem.reaction!,Y,U[:,inode])
        res_react=DiffResults.value(result_r)
        jac_react=DiffResults.jacobian(result_r)

        # Evaluate source term
        fvsystem.source!(src,geom.Nodes[:,inode])

        # Evaluate & differentiate storage term
        result_s=ForwardDiff.jacobian!(result_s,fvsystem.storage!,Y,U[:,inode])
        res_stor=DiffResults.value(result_s)
        jac_stor=DiffResults.jacobian(result_s)
       
        # Evaluate storage term for old timestep
        fvsystem.storage!(oldstor,UOld[:,inode])
        
        # Assembly results and jacobians
        for i=1:number_of_species
            F[i,inode]+=geom.NodeFactors[inode]*(res_react[i]-src[i] + (res_stor[i]-oldstor[i])*tstepinv)
            for j=1:number_of_species
                M[iblock+i,iblock+j]+=geom.NodeFactors[inode]*(jac_react[i,j]+ jac_stor[i,j]*tstepinv)
            end
        end
        iblock+=number_of_species
    end
    
    # Create result struct for flux evaluation
    result=DiffResults.DiffResult(Vector{Float64}(undef,number_of_species),Matrix{Float64}(undef,number_of_species,2*number_of_species))
    Y=Array{Float64,1}(undef,number_of_species)
    UKL=Array{Float64,1}(undef,2*number_of_species)
    # Assemble main part
    for iedge=1:nedges
        K=geom.Edges[1,iedge]
        L=geom.Edges[2,iedge]
        # Set up argument for fluxwrap!
        UKL[1:number_of_species]=U[:,K]
        UKL[number_of_species+1:2*number_of_species]=U[:,L]
        result=ForwardDiff.jacobian!(result,fluxwrap!,Y,UKL)

        res=DiffResults.value(result)
        jac=DiffResults.jacobian(result)

        # Assemble flux data
        F[:,K]+=res*geom.EdgeFactors[iedge]
        F[:,L]-=res*geom.EdgeFactors[iedge]
        kblock=(K-1)*number_of_species
        lblock=(L-1)*number_of_species
        jl=number_of_species+1
        for jk=1:number_of_species
            for ik=1:number_of_species
                M[kblock+ik,kblock+jk]+=jac[ik,jk]*geom.EdgeFactors[iedge]
                M[kblock+ik,lblock+jk]+=jac[ik,jl]*geom.EdgeFactors[iedge]
                M[lblock+ik,kblock+jk]-=jac[ik,jk]*geom.EdgeFactors[iedge]
                M[lblock+ik,lblock+jk]-=jac[ik,jl]*geom.EdgeFactors[iedge]
            end
            jl+=1
        end
    end
    
    # Assemble boundary conditions
    # Dirichlet conditions are handeld via penalty method, 
    # for these, solution vectors need to be initialized
    # appropriately
    nbnodes=length(geom.BoundaryNodes)
    for ibnode=1:nbnodes
        inode=geom.BoundaryNodes[ibnode]
        ibreg=geom.BoundaryRegions[ibnode]
        iblock=(inode-1)*number_of_species
        for ispec=1:number_of_species
            fac=fvsystem.boundary_factors[ispec,ibreg]
            if fac!=Dirichlet
                fac*=geom.BoundaryNodeFactors[ibnode]
            end
            F[ispec,inode]+=fac*(U[ispec,inode]-fvsystem.boundary_values[ispec,ibreg])
            M[iblock+ispec,iblock+ispec]+=fac
        end
    end
end



"""
    Actual solver function implementation
"""
function _solve(fvsystem::TwoPointFluxFVMSystem, oldsol::Array{Float64,2},control::FVMNewtonControl, tstep::Float64)
    
    nunknowns=fvsystem.geometry.NumberOfNodes*fvsystem.number_of_species

    solution=copy(oldsol)
    inidirichlet!(fvsystem,solution)
    solution_r=reshape(solution,nunknowns)
    residual=fvsystem.residual
    update=fvsystem.update

    # Newton iteration (quick and dirty...)
    oldnorm=1.0
    converged=false
    if control.verbose
        @printf("Start newton iteration: %s:%d\n", basename(@__FILE__),@__LINE__)
    end
    for ii=1:control.maxiter
        eval_and_assemble(fvsystem,solution,oldsol,tstep)
        
        # Sparse LU factorization
        # Here, we seem miss the possibility to re-use the 
        # previous symbolic information
        # !!! may be there is such a call
        lufact=LinearAlgebra.lu(fvsystem.matrix)
        
        # LU triangular solve gives Newton update
        # !!! is there a version wich does not allocate ?
        update=lufact\residual # DU is the Newton update

        # vector expressions would allocate here...
        for i=1:nunknowns
            solution_r[i]-=control.damp*update[i]
        end

        norm=LinearAlgebra.norm(update)/nunknowns
        if control.verbose
            @printf("  it=%03d norm=%.5e cont=%.5e\n",ii,norm, norm/oldnorm)
        end
        if norm<control.tolerance
            converged=true
            break
        end
        
        oldnorm=norm
    end
    if !converged
        println("error: no convergence")
        exit(1)
    end
    return solution
end

"""
    System solver wrapper allowing to dispatch timing
"""
function solve(fvsystem::TwoPointFluxFVMSystem, # Finite volume system
               oldsol::Array{Float64,2}; # old time step solution resp. initial value
               control=FVMNewtonControl(), # Newton solver control information
               tstep::Float64=Inf) # Time step size. Inf means  stationary solution
    if control.verbose
        @time begin
            retval= _solve(fvsystem,oldsol,control,tstep)
        end
        return retval
    else
        return _solve(fvsystem,oldsol,control,tstep)
    end
end
