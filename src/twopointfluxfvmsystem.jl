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



abstract type FVMPhysics end

struct TwoPointFluxFVMSystem
    geometry::FVMGraph
    nspec::Int64
    source!::Function
    reaction!::Function
    flux!::Function
    dirichlet_values::Array{Float64,2}
    matrix::SparseArrays.SparseMatrixCSC
    residual::Array{Float64,1}
    update::Array{Float64,1}
    function TwoPointFluxFVMSystem(geometry::FVMGraph,problem::FVMPhysics)
        nspec=problem.nspec(problem)
        source!(y,x)=problem.source(problem,y,x)
        flux!(y,uk,ul)=problem.flux(problem,y,uk,ul)
        reaction!(y,x)=problem.reaction(problem,y,x)
        # Set up solution data
        matrix=SparseArrays.spzeros(geometry.NumberOfNodes*nspec,geometry.NumberOfNodes*nspec) # Jacobi matrix
        residual=Array{Float64,1}(undef,geometry.NumberOfNodes*nspec)
        update=Array{Float64,1}(undef,geometry.NumberOfNodes*nspec)
        dirichlet_values=Array{Float64,2}(undef,nspec,length(geometry.BPoints))
        new(geometry,nspec,source!,reaction!,flux!,dirichlet_values,matrix,residual,update)
    end
end

function unknowns(fvsystem::TwoPointFluxFVMSystem)
    return Array{Float64,2}(undef,fvsystem.nspec,fvsystem.geometry.NumberOfNodes)
end





#
# Nonlinear operator evaluation + Jacobian assembly
#
function eval_and_assemble(fvsystem::TwoPointFluxFVMSystem,U)
    dirichlet_penalty=1.0e30
    
    function fluxwrap!(y,u)
        fvsystem.flux!(y,u[1:nspec],u[nspec+1:2*nspec])
    end
    
    geom=fvsystem.geometry
    nnodes=geom.NumberOfNodes
    nspec=fvsystem.nspec
    nedges=size(geom.Edges,2)
    M=fvsystem.matrix
    F=reshape(fvsystem.residual,nspec,nnodes)
    #  for K=1...n
    #  f_K = sum_(L neigbor of K) eps (U[K]-U[L])*edgefac[K,L]
    #        + (reaction(U[K])- source(X[K]))*nodefac[K]
    # M is correspondig Jacobi matrix of derivatives. 
    
    # Reset matrix
    M.nzval.=0.0
    F.=0.0
    # Assemble nonlinear term + source using autodifferencing via ForwardDiff
    result=DiffResults.DiffResult(Vector{Float64}(undef,nspec),Matrix{Float64}(undef,nspec,nspec))
    Y=Array{Float64}(undef,nspec)
    iblock=0
    for inode=1:nnodes
        result=ForwardDiff.jacobian!(result,fvsystem.reaction!,Y,U[:,inode])
        res=DiffResults.value(result)
        fvsystem.source!(Y,geom.Points[:,inode])
        F[:,inode]=geom.NodeFactors[inode]*(res-Y)
        jac=DiffResults.jacobian(result)
        for i=1:nspec
            for j=1:nspec
                M[iblock+i,iblock+j]+=geom.NodeFactors[inode]*jac[i,j]
            end
        end
        iblock+=nspec
    end
    
    result=DiffResults.DiffResult(Vector{Float64}(undef,nspec),Matrix{Float64}(undef,nspec,2*nspec))
    Y=Array{Float64,1}(undef,nspec)
    UKL=Array{Float64,1}(undef,2*nspec)
    # Assemble main part
    for iedge=1:nedges
        K=geom.Edges[1,iedge]
        L=geom.Edges[2,iedge]
        UKL[1:nspec]=U[:,K]
        UKL[nspec+1:2*nspec]=U[:,L]
        result=ForwardDiff.jacobian!(result,fluxwrap!,Y,UKL)
        res=DiffResults.value(result)
        jac=DiffResults.jacobian(result)
        F[:,K]+=res*geom.EdgeFactors[iedge]
        F[:,L]-=res*geom.EdgeFactors[iedge]

        kblock=(K-1)*nspec
        lblock=(L-1)*nspec
        for ik=1:nspec
            jk=1
            jl=nspec+1
            for j=1:nspec
                M[kblock+ik,kblock+jk]+=jac[ik,jk]*geom.EdgeFactors[iedge]
                M[kblock+ik,lblock+jk]+=jac[ik,jl]*geom.EdgeFactors[iedge]
                M[lblock+ik,kblock+jk]-=jac[ik,jk]*geom.EdgeFactors[iedge]
                M[lblock+ik,lblock+jk]-=jac[ik,jl]*geom.EdgeFactors[iedge]
                jk+=1
                jl+=1
            end
        end
    end
    
    # Assemble boundary value part using Dirichlet penalty method
    nbc=length(geom.BPoints)
    for i=1:nbc
        ibc=geom.BPoints[i]
        F[:,ibc]+=dirichlet_penalty*(U[:,ibc]-fvsystem.dirichlet_values[:,i])
        iblock=(ibc-1)*nspec
        for ib=1:nspec
            M[iblock+ib,iblock+ib]+=dirichlet_penalty
        end
    end
end




# Input: geometry, physics
# Output: Array of solution values

function solve(fvsystem::TwoPointFluxFVMSystem, inival::Array{Float64,2};control=FVMNewtonControl())
    
    nunknowns=fvsystem.geometry.NumberOfNodes*fvsystem.nspec
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

