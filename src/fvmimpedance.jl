abstract type AbstractImpedanceSystem{Tv <: Number} end

mutable struct ImpedanceSystem{Tv} <: AbstractImpedanceSystem{Tv}
    sys::AbstractSystem{Tv}
    storderiv::AbstractMatrix{Tv}
    matrix::AbstractMatrix{Complex{Tv}}
    ibc::Integer
    ispec::Integer
    F::AbstractMatrix{Complex{Tv}}
    U0::AbstractMatrix{Tv}
    ImpedanceSystem{Tv}() where Tv =new()
end




function ImpedanceSystem(sys::AbstractSystem{Tv}, U0::AbstractMatrix, xispec,xibc) where Tv
    this=ImpedanceSystem{Tv}()

    this.sys=sys
    this.storderiv=spzeros(Tv,num_dof(sys), num_dof(sys))
    this.matrix=SparseMatrixCSC(sys.matrix.m,
                                sys.matrix.n,
                                sys.matrix.colptr,
                                sys.matrix.rowval,
                                complex(sys.matrix.nzval)
                                )
    this.ispec=xispec
    this.ibc=xibc
    this.U0=U0
    
    F=complex(unknowns(sys))
    this.F=F
    
    grid=sys.grid
    
    physics::Physics=sys.physics
    node::Node=Node{Tv}()
    bnode::BNode=BNode{Tv}()
    nspecies::Int32=num_species(sys)
    
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    bnode_factors=zeros(Tv,num_nodes_per_bface(grid))
    
    @inline function storagewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.storage(physics,node,y,u)
    end
    
    @inline function bstoragewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.bstorage(physics,bnode,y,u)
    end
    
    F.=0.0
    UK=Array{Tv,1}(undef,nspecies)
    Y=Array{Tv,1}(undef,nspecies)
    
    # structs holding diff results for storage, reaction,  flux ...
    result_s=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    
    
    # Main cell loop
    for icell=1:num_cells(grid)
        # set up form factors
        cellfactors!(grid,icell,node_factors,edge_factors)
        
        # set up data for callbacks
        node.region=reg_cell(grid,icell)
        
        for inode=1:num_nodes_per_cell(grid)
            fill!(node,grid, inode,icell)
            @views begin
                UK[1:nspecies]=U0[:,node.index]
            end
            
            # Evaluate & differentiate storage term
            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK)
            jac_stor=DiffResults.jacobian(result_s)


            K=node.index
            for idof=firstnodedof(F,K):lastnodedof(F,K)
                ispec=spec(F,idof,K)
                for jdof=firstnodedof(F,K):lastnodedof(F,K)
                    jspec=spec(F,jdof,K)
                    addnz(this.storderiv,idof,jdof,jac_stor[ispec,jspec],node_factors[inode])
                end
            end

        end
        
    end
    
    for ibface=1:num_bfaces(grid)
        bfacefactors!(grid,ibface,bnode_factors)
        ibreg=grid.bfaceregions[ibface]
        bnode.region=ibreg
        for ibnode=1:num_nodes_per_bface(grid)
            @views begin
                fill!(bnode,grid,ibnode,ibface)
                
                UK[1:nspecies]=U0[:,bnode.index]
            end
            if ibreg==xibc
                for ispec=1:nspecies # should involve only rspecies
                    if ispec!=xispec
                        continue
                    end
                    idof=dof(F,ispec,bnode.index)
                    if idof>0 
                        fac=sys.boundary_factors[ispec,ibreg]
                        val=sys.boundary_values[ispec,ibreg]
                        if fac==Dirichlet
                            F[ispec,bnode.index]+=fac
                        else
                            F[ispec,bnode.index]+=fac*bnode_factors[ibnode]
                        end
                    end
                end
            end
            
            if isdefined(physics, :bstorage) # should involve only bspecies
                # Evaluate & differentiate storage term
                ForwardDiff.jacobian!(result_s,bstoragewrap,Y,UK)
                res_bstor=DiffResults.value(result_s)
                jac_bstor=DiffResults.jacobian(result_s)


                K=bnode.index
                for idof=firstnodedof(F,K):lastnodedof(F,K)
                    ispec=spec(F,idof,K)
                    for jdof=firstnodedof(F,K):lastnodedof(F,K)
                        jspec=spec(F,jdof,K)
                        addnz(this.storderiv,idof,jdof,jac_bstor[ispec,jspec],bnode_factors[ibnode])
                    end
                end
            end
            
        end
    end

    
    return this
end

# correct!
function solve(this::ImpedanceSystem{Tv}, ω) where Tv
    iω=ω*1im
    nspecies::Int32=num_species(this.sys)
    for i=1:length(this.matrix.nzval)
        this.matrix.nzval[i]=complex(this.sys.matrix.nzval[i])
    end
    U0=this.U0
    for inode=1:num_nodes(this.sys.grid)
        for idof=firstnodedof(U0,inode):lastnodedof(U0,inode)
            for jdof=firstnodedof(U0,inode):lastnodedof(U0,inode)
                addnz(this.matrix,idof,jdof,this.storderiv[idof,jdof],iω)
            end
        end
    end
    U=copy(this.F)
    lufact=LinearAlgebra.lu(this.matrix)
    ldiv!(values(U),lufact,values(this.F))
    return U
end


function integrate(this::AbstractImpedanceSystem{Tv},tf::Vector{Tv},ω::Tv,UZ::AbstractMatrix{Complex{Tv}}) where Tv
    iω=ω*1im
    sys=this.sys
    grid=this.sys.grid
    nspecies=num_species(sys)
    integral=zeros(Complex{Tv},nspecies)
    res=zeros(Tv,nspecies)
    stor=zeros(Tv,nspecies)
    node=Node{Tv}()
    edge=Edge{Tv}()
    node_factors=zeros(Tv,num_nodes_per_cell(grid))
    edge_factors=zeros(Tv,num_edges_per_cell(grid))
    edge_cutoff=1.0e-12
    U0=this.U0
    physics=sys.physics


    @inline function fluxwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        @views physics.flux(physics,edge,y,u[K1:KN],u[L1:LN])
    end

    @inline function reactionwrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.reaction(physics,node,y,u)
    end

    @inline function storagewrap(y::AbstractVector, u::AbstractVector)
        y.=0
        physics.storage(physics,node,y,u)
    end

    K1::Int32=1
    KN::Int32=nspecies
    L1::Int32=nspecies+1
    LN::Int32=2*nspecies

    # istoragew+= S'*uz*nodefac*tfc
    # ireactionw+= R'*uz*nodefac*tfc
    # iflux+=(FK'*uzk+FL'*uzl)*(tfcK-tfcL)*edgefac

    # integral=istoragew+ireactionw+iflux
    
    # structs holding diff results for storage, reaction,  flux ...
    result_r=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_s=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    result_flx=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,2*nspecies))

    jac_react=zeros(Tv,nspecies,nspecies)

    # Arrays for gathering solution data
    UK=Array{Tv,1}(undef,nspecies)
    UKL=Array{Tv,1}(undef,2*nspecies)
    Y=Array{Tv,1}(undef,nspecies)

    for icell=1:num_cells(grid)
        cellfactors!(grid,icell,node_factors,edge_factors)
        for inode=1:num_nodes_per_cell(grid)
            fill!(node,grid,inode,icell)
            @views  UK[1:nspecies]=U0[:,node.index]

            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK)
            jac_stor=DiffResults.jacobian(result_s)
            if isdefined(physics, :reaction)
                ForwardDiff.jacobian!(result_r,reactionwrap,Y,UK)
                jac_react=DiffResults.jacobian(result_r)
            end
            K=node.index
            for idof=firstnodedof(U0,K):lastnodedof(U0,K)
                ispec=spec(U0,idof,K)
                for jdof=firstnodedof(U0,K):lastnodedof(U0,K)
                    jspec=spec(U0,jdof,K)
                    
                    integral[ispec]+=(node_factors[inode]*tf[K]
                                      *(jac_stor[ispec,jspec]*iω+jac_react[ispec,jspec])*UZ[jspec,K])
                end
            end
        end


        for iedge=1:num_edges_per_cell(grid)
            if edge_factors[iedge]<edge_cutoff
                continue
            end
            fill!(edge,grid,iedge,icell)
            @views UKL[K1:KN]=U0[:,edge.nodeK]
            @views UKL[L1:LN]=U0[:,edge.nodeL]

            ForwardDiff.jacobian!(result_flx,fluxwrap,Y,UKL)
            jac=DiffResults.jacobian(result_flx)
            
            K=edge.nodeK
            L=edge.nodeL
            fac=edge_factors[iedge]*(tf[K]-tf[L])
            for idofK=firstnodedof(U0,K):lastnodedof(U0,K)
                ispec=spec(U0,idofK,K)
                idofL=dof(U0,ispec,L)
                if idofL==0
                    continue
                end
                for jdofK=firstnodedof(U0,K):lastnodedof(U0,K)
                    jspec=spec(U0,jdofK,K)
                    jdofL=dof(U0,jspec,L)
                    if jdofL==0
                        continue
                    end
                    integral[ispec]+=fac*(jac[ispec,jspec]*UZ[jspec,K]+jac[ispec,jspec+nspecies]*UZ[jspec,L])

                end
            end
        end
        
    end
    
    return integral
end
