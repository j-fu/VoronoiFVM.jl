#
# Interface to VoronoiFVMDiffEq.jl
#
"""
    $(TYPEDEF)

History information for DiffEqHistory
$(TYPEDFIELDS)
"""
Base.@kwdef mutable struct DiffEqHistory 
    """ number of combined jacobi/rhs evaluations"""
    nd=0
    """ number of combined jacobi evaluations"""
    njac=0
    """ number of rhs evaluations"""
    nf=0
end

details(h::DiffEqHistory)=(nd=h.nd,njac=h.njac,nf=h.nf)
Base.summary(h::DiffEqHistory)=details(h)


#Evaluate functiaon and Jacobian at u if they have not been evaluated before at u.
#See https://github.com/SciML/DifferentialEquations.jl/issues/521 for discussion of another way to do this.
function _eval_res_jac!(sys,u,t)
    uhash=hash(u)
    if uhash!=sys.uhash
        ur=reshape(u,sys)
        eval_and_assemble(sys,ur,ur,sys.residual,value(t),Inf,0.0,zeros(0))
        sys.uhash=uhash
        sys.history.nd+=1
    end
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
rhs function for DifferentialEquations.jl.
"""
function eval_rhs!(du, u, sys,t)
    _eval_res_jac!(sys,u,t)
    du.=-vec(sys.residual)
    sys.history.nf+=1
    nothing
end

"""
$(SIGNATURES)

Interpret the  discrete problem as an ODE/DAE problem. Provide the 
jacobi matrix calculation function for DifferentialEquations.jl.
"""
function eval_jacobian!(J, u, sys,t)
    _eval_res_jac!(sys,u,t)
    # Need to implement broadcast for ExtendableSparse.
    J.=-sys.matrix.cscmatrix
    sys.history.njac+=1
    nothing
end

"""
$(SIGNATURES)

Calculate the mass matrix for use with DifferentialEquations.jl.
"""
function mass_matrix(system::AbstractSystem{Tv, Ti}) where {Tv, Ti}
    physics=system.physics
    grid=system.grid
    data=physics.data
    node=Node(system)
    bnode=BNode(system)
    nspecies=num_species(system)
    cellnodefactors=system.cellnodefactors
    bfacenodefactors=system.bfacenodefactors
    bgeom=grid[BFaceGeometries][1]
    nbfaces=num_bfaces(grid)

    if isdata(data)
        isbstorage=(physics.bstorage!=nofunc)

        storagewrap= function(y, u)
            y.=0
            physics.storage(y,u,node,data)
            nothing
        end
        
        bstoragewrap=function(y, u)
            y.=0
            physics.bstorage(y,u,bnode,data)
            nothing
        end
        
    else
        isbstorage=(physics.bstorage!=nofunc2)

        storagewrap= function(y, u)
            y.=0
            physics.storage(y,u,node)
            nothing
        end
        bstoragewrap=function(y, u)
            y.=0
            physics.bstorage(y,u,bnode)
            nothing
        end
    end        

    result_s=DiffResults.DiffResult(Vector{Tv}(undef,nspecies),Matrix{Tv}(undef,nspecies,nspecies))
    Y=Array{Tv,1}(undef,nspecies)
    UK=Array{Tv,1}(undef,nspecies)
    cfg_s=ForwardDiff.JacobianConfig(storagewrap, Y, UK)
    cfg_bs=ForwardDiff.JacobianConfig(bstoragewrap, Y, UK)

    U=unknowns(system,inival=1)
    M=unknowns(system,inival=0)
    
    geom=grid[CellGeometries][1]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]
    bfacenodes=grid[BFaceNodes]

    for icell=1:num_cells(grid)
        for inode=1:num_nodes(geom)
            _fill!(node,inode,icell)
            for ispec=1:nspecies
                UK[ispec]=U[ispec,node.index]
            end
            ForwardDiff.jacobian!(result_s,storagewrap,Y,UK,cfg_s)
            res_stor=DiffResults.value(result_s)
            jac_stor=DiffResults.jacobian(result_s)
            K=node.index
            for idof=_firstnodedof(M,K):_lastnodedof(M,K)
                ispec=_spec(M,idof,K)
                M[ispec,K]+=cellnodefactors[inode,icell]*jac_stor[ispec,ispec]
            end
        end
    end
    
    if isbstorage
        for ibface=1:nbfaces
            for ibnode=1:num_nodes(bgeom)
                _fill!(bnode,ibnode,ibface)
                for ispec=1:nspecies
                    UK[ispec]=U[ispec,bnode.index]
                end
                ForwardDiff.jacobian!(result_s,bstoragewrap,Y,UK,cfg_bs)
                jac_bstor=DiffResults.jacobian(result_s)
                for idof=_firstnodedof(M,K):_lastnodedof(M,K)
                    ispec=_spec(M,idof,K)
                    M[ispec,K]+=bfacenodefactors[ibnode,ibface]*jac_bstor[ispec,ispec]
                end
                
            end
        end
    end
    
    return Diagonal(vec(M))
end

"""
$(SIGNATURES)
Prepare system for use with VoronoiFVMDiffEq.

- `jacval`: value at which to evaluate jacobian to obtatin prototype
- `tjac`: time moment for jacobian
 
Returns a prototype for the jacobian.
"""
function prepare_diffeq!(sys,jacval,tjac)
    sys.history=DiffEqHistory()
    _complete!(sys, create_newtonvectors=true)
    _eval_res_jac!(sys,jacval,tjac)
    flush!(sys.matrix)
    sys.matrix.cscmatrix
end
