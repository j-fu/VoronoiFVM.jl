# Autodiff magic, get it via Pkg.add
using ForwardDiff, DiffResults, Printf


# These are in the standard distro
using SparseArrays
using LinearAlgebra

mutable struct FVMNewtonControl
     tolerance::Float64
     damp::Float64
     maxiter::Int32
     function FVMNewtonControl()
         new(1.0e-10,1.0,100)
     end
end



abstract type FVMParameters end

mutable struct DefaultParameters <: FVMParameters
    number_of_species::Int64
end

function default_source!(this::FVMParameters, f,x)
    for i=1:this.number_of_species
        f[i]=0
    end
end

function default_reaction!(this::FVMParameters, f,u)
    for i=1:this.number_of_species
        f[i]=0
    end
end

function default_flux!(this::FVMParameters, f,uk,ul)
    for i=1:this.number_of_species
        f[i]=uk[i]-ul[i]
    end
end



struct TwoPointFluxFVMSystem
    geometry::FVMGraph
    number_of_species::Int64
    source!::Function
    reaction!::Function
    flux!::Function
    dirichlet_values::Array{Float64,2}
    matrix::SparseArrays.SparseMatrixCSC
    residual::Array{Float64,1}
    update::Array{Float64,1}
    function TwoPointFluxFVMSystem(geometry::FVMGraph; 
                                   parameters::FVMParameters=DefaultParameters(1),
                                   source::Function=default_source!,
                                   reaction::Function=default_reaction!,
                                   flux::Function=default_flux!)
        number_of_species=parameters.number_of_species
        _source!(y,x)=source(parameters,y,x)
        _flux!(y,uk,ul)=flux(parameters,y,uk,ul)
        _reaction!(y,x)=reaction(parameters,y,x)
        # Set up solution data
        matrix=SparseArrays.spzeros(geometry.NumberOfNodes*number_of_species,geometry.NumberOfNodes*number_of_species) # Jacobi matrix
        residual=Array{Float64,1}(undef,geometry.NumberOfNodes*number_of_species)
        update=Array{Float64,1}(undef,geometry.NumberOfNodes*number_of_species)
        dirichlet_values=Array{Float64,2}(undef,number_of_species,length(geometry.BPoints))
        new(geometry,number_of_species,_source!,_reaction!,_flux!,dirichlet_values,matrix,residual,update)
    end
end

function unknowns(fvsystem::TwoPointFluxFVMSystem)
    return Array{Float64,2}(undef,fvsystem.number_of_species,fvsystem.geometry.NumberOfNodes)
end





#
# Nonlinear operator evaluation + Jacobian assembly
#
function eval_and_assemble(fvsystem::TwoPointFluxFVMSystem,U)
    dirichlet_penalty=1.0e30
    
    function fluxwrap!(y,u)
        fvsystem.flux!(y,u[1:number_of_species],u[number_of_species+1:2*number_of_species])
    end
    
    geom=fvsystem.geometry
    nnodes=geom.NumberOfNodes
    number_of_species=fvsystem.number_of_species
    nedges=size(geom.Edges,2)
    M=fvsystem.matrix
    F=reshape(fvsystem.residual,number_of_species,nnodes)
    #  for K=1...n
    #  f_K = sum_(L neigbor of K) eps (U[K]-U[L])*edgefac[K,L]
    #        + (reaction(U[K])- source(X[K]))*nodefac[K]
    # M is correspondig Jacobi matrix of derivatives. 
    
    # Reset matrix
    M.nzval.=0.0
    F.=0.0
    # Assemble nonlinear term + source using autodifferencing via ForwardDiff
    result=DiffResults.DiffResult(Vector{Float64}(undef,number_of_species),Matrix{Float64}(undef,number_of_species,number_of_species))
    Y=Array{Float64}(undef,number_of_species)
    iblock=0
    for inode=1:nnodes
        result=ForwardDiff.jacobian!(result,fvsystem.reaction!,Y,U[:,inode])
        res=DiffResults.value(result)
        fvsystem.source!(Y,geom.Points[:,inode])
        F[:,inode]=geom.NodeFactors[inode]*(res-Y)
        jac=DiffResults.jacobian(result)
        for i=1:number_of_species
            for j=1:number_of_species
                M[iblock+i,iblock+j]+=geom.NodeFactors[inode]*jac[i,j]
            end
        end
        iblock+=number_of_species
    end
    
    result=DiffResults.DiffResult(Vector{Float64}(undef,number_of_species),Matrix{Float64}(undef,number_of_species,2*number_of_species))
    Y=Array{Float64,1}(undef,number_of_species)
    UKL=Array{Float64,1}(undef,2*number_of_species)
    # Assemble main part
    for iedge=1:nedges
        K=geom.Edges[1,iedge]
        L=geom.Edges[2,iedge]
        UKL[1:number_of_species]=U[:,K]
        UKL[number_of_species+1:2*number_of_species]=U[:,L]
        result=ForwardDiff.jacobian!(result,fluxwrap!,Y,UKL)
        res=DiffResults.value(result)
        jac=DiffResults.jacobian(result)
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
    
    # Assemble boundary value part using Dirichlet penalty method
    nbc=length(geom.BPoints)
    for i=1:nbc
        ibc=geom.BPoints[i]
        F[:,ibc]+=dirichlet_penalty*(U[:,ibc]-fvsystem.dirichlet_values[:,i])
        iblock=(ibc-1)*number_of_species
        for ib=1:number_of_species
            M[iblock+ib,iblock+ib]+=dirichlet_penalty
        end
    end
end




# Input: geometry, physics
# Output: Array of solution values

function solve(fvsystem::TwoPointFluxFVMSystem, inival::Array{Float64,2};control=FVMNewtonControl())
    
    nunknowns=fvsystem.geometry.NumberOfNodes*fvsystem.number_of_species
    solution=copy(inival)
    solution_r=reshape(solution,nunknowns)
    residual=fvsystem.residual
    update=fvsystem.update

    # Newton iteration (quick and dirty...)
    oldnorm=1.0
    converged=false
    for ii=1:control.maxiter
        eval_and_assemble(fvsystem,solution)
        
        # Sparse LU factorization
        # Here, we miss the possibility to re-use the 
        # previous symbolic information
        # !!! may be there is such a call
        lufact=LinearAlgebra.lu(fvsystem.matrix)
        
        # LU triangular solve gives Newton update
        # !!! is there a version wich does not allocate ?
        update=lufact\residual # DU is the Newton update

        # vector expressions would allocate, we might
        # miss 
        for i=1:nunknowns
            solution_r[i]-=control.damp*update[i]
        end

        norm=LinearAlgebra.norm(update)
        @printf("it=%03d norm=%.5e cont=%.5e\n",ii,norm, norm/oldnorm)
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

