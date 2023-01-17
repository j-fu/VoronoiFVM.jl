"""
$(TYPEDEF)

Abstract type for impedance system.
"""
abstract type AbstractImpedanceSystem{Tv <: Number} end


"""
$(TYPEDEF)

Concrete type for impedance system.

$(TYPEDFIELDS)
"""
mutable struct ImpedanceSystem{Tv} <: AbstractImpedanceSystem{Tv}
    
    """
    Nonzero pattern of time domain system matrix
    """
    sysnzval::AbstractVector{Complex{Tv}}

    """
    Derivative of storage term
    """
    storderiv::AbstractMatrix{Tv}

    """
    Complex matrix of impedance system
    """
    matrix::AbstractMatrix{Complex{Tv}}

    """
    Right hand side of impedance system
    """
    F::AbstractMatrix{Complex{Tv}}

    """
    Stationary state
    """
    U0::AbstractMatrix{Tv}

    ImpedanceSystem{Tv}() where Tv =new()
end


"""
$(SIGNATURES)

Construct impedance system from time domain system `sys` and steady state solution `U0`
"""
function ImpedanceSystem(system::AbstractSystem{Tv,Tc,Ti}, U0::AbstractMatrix; λ0=0.0) where {Tv,Tc,Ti}
    residual=copy(U0)

    if system.num_parameters>0
        params=[λ0]
    else
        params=zeros(0)
    end
    
    # Ensure that system.matrix contains the jacobian at U0
    # Moreover, here we use the fact that for large time step sizes,
    # the contribution from the time derivative (which should not belong to the
    # main part of the impedance matrix) is small. We also pass the steady state
    # value as the "old  timestep" value.
    # An advantage of this approach is the fact that this way, we get the
    # nonzero pattern for the iω term right (as opposite to passing Inf as time step size)
    eval_and_assemble(system,U0,U0,residual,0.0,1.0e30,0.0,params)

    impedance_system=ImpedanceSystem{Tv}()

    # catch the nonzero values of the system matrix
    impedance_system.sysnzval=complex(nonzeros(system.matrix))
    m,n=size(system.matrix)

    # initialize storage derivative with zero
    impedance_system.storderiv=spzeros(Tv,m,n)

    # create sparse matrix with the same nonzero pattern of original matrix
    impedance_system.matrix=SparseMatrixCSC(m,n,
                                SparseArrays.getcolptr(system.matrix),
                                SparseArrays.getrowval(system.matrix),
                                copy(impedance_system.sysnzval)
                                )
    impedance_system.U0=U0

    # Initialize right hand side of impedance system
    impedance_system.F=unknowns(Complex{Float64},system)
 
    matrix=impedance_system.matrix
    storderiv=impedance_system.storderiv
    F=impedance_system.F
    grid=system.grid
    
    physics=system.physics

    # Prepare calculation of derivative of the storage part
    data=physics.data
    node=Node(system)
    bnode=BNode(system)
    nspecies=num_species(system)
        
    UK=Array{Tv,1}(undef,nspecies)
    stor_eval = ResJacEvaluator(physics,:storage,UK,node,nspecies)
    bstor_eval = ResJacEvaluator(physics,:bstorage,UK,node,nspecies)

    F.=0.0
    if system.num_parameters>0
        F.-=system.dudp[1]
    end
    

    geom=grid[CellGeometries][1]

    bgeom=grid[BFaceGeometries][1]

    asm_res(idof,ispec)=nothing
    asm_param(idof,ispec,iparam)= nothing

    # Interior cell loop for building up storage derivative
    for icell=1:num_cells(grid)
        for inode=1:num_nodes(geom)
            _fill!(node,inode,icell)
            @views UK[1:nspecies]=U0[:,node.index]
            
            # Evaluate & differentiate storage term at U0
            evaluate!(stor_eval,UK)
            jac_stor=jac(stor_eval)

            # Sort it into storderiv matrix.
            K=node.index
            asm_jac(idof,jdof,ispec,jspec)=updateindex!(storderiv,+,jac_stor[ispec,jspec]*system.cellnodefactors[inode,icell],idof,jdof)
            assemble_res_jac(node,system,asm_res,asm_jac,asm_param)
        end
    end
    
    # Boundary face loop for building up storage derivative
    # and right hand side contribution from boundary condition
    if isdefined(physics, :bstorage) # should involve only bspecies
        for ibface=1:num_bfaces(grid)
            for ibnode=1:num_nodes(bgeom)
                _fill!(bnode,ibnode,ibface)
                @views UK[1:nspecies]=U0[:,bnode.index]
                # Evaluate & differentiate storage term
                evaluate!(bstor_eval,UK)
                jac_bstor=jac(bstor_eval)
                K=bnode.index
                asm_jac(idof,jdof,ispec,jspec)=updateindex!(storderiv,+,jac_bstor[ispec,jspec]*system.bfacenodefactors[ibnode,ibface],idof,jdof)
                assemble_res_jac(bnode,system,asm_res,asm_jac,asm_param)
            end
        end
    end

    return impedance_system
end    

"""
$(SIGNATURES)

Construct impedance system from time domain system `sys` and steady state solution `U0`
under the assumption of a periodic perturbation of species `excited_spec` at  boundary `excited_bc`.
"""
function ImpedanceSystem(system::AbstractSystem{Tv,Tc,Ti}, U0::AbstractMatrix, excited_spec,excited_bc) where {Tv,Tc,Ti}
    impedance_system=ImpedanceSystem(system,U0)
    grid=system.grid
    bgeom=grid[BFaceGeometries][1]
    bnode=BNode(system)
    bfaceregions::Vector{Ti}=grid[BFaceRegions]
    nspecies=num_species(system)
    F=impedance_system.F
    
    
    for ibface=1:num_bfaces(grid)
        ibreg=bfaceregions[ibface]
        for ibnode=1:num_nodes(bgeom)
            _fill!(bnode,ibnode,ibface)

            # Set right hand side for excited boundary conditions
            # We don't need to put the penalty term to storderiv
            # as it is already in the main part of the matrix
            if ibreg==excited_bc
                for ispec=1:nspecies 
                    if ispec!=excited_spec
                        continue
                    end
                    idof=dof(F,ispec,bnode.index)
                    if idof>0 
                        fac=system.boundary_factors[ispec,ibreg]
                        if fac==Dirichlet
                            F[ispec,bnode.index]+=fac
                        else
                            F[ispec,bnode.index]+=fac*system.bfacenodefactors[ibnode,ibface]
                        end
                    end
                end
            end
        end
    end
    
    return impedance_system
end

"""
$(SIGNATURES)

Create a vector of unknowns of the impedance system
"""
unknowns(impedance_system::ImpedanceSystem{Tv}) where Tv=copy(impedance_system.F)


"""
$(SIGNATURES)

Solve the impedance system for given frequency `ω`.
"""
function solve!(UZ::AbstractMatrix{Complex{Tv}},impedance_system::ImpedanceSystem{Tv}, ω) where Tv
    iω=ω*1im
    matrix=impedance_system.matrix
    storderiv=impedance_system.storderiv
    # Reset the matrix nonzero values
    nzval=nonzeros(matrix)
    nzval.=impedance_system.sysnzval
    # Add the ω dependent term to main diagonal.
    # This makes the implicit assumption that storderiv does not
    # introduce new matrix entries.

    rowval=storderiv.rowval
    colptr=storderiv.colptr
    inzval=storderiv.nzval

    for i=1:length(colptr)-1
        for k=colptr[i]:colptr[i+1]-1
            j=rowval[k]
            updateindex!(matrix,+,inzval[k]*iω,i,j)
        end
    end
    
    lufact=LinearAlgebra.lu(matrix)
    ldiv!(values(UZ),lufact,values(impedance_system.F))
end


"""
$(SIGNATURES)

Calculate the derivative of the scalar measurement functional at steady state U0

Usually, this functional is  a test function integral.  Initially,
we assume that its value depends on all unknowns of the system.
"""
function measurement_derivative(system::AbstractSystem,measurement_functional,U0)

    # Create a sparse 1×ndof matrix assuming that the functional
    # depends on all unknowns in the system
    ndof=num_dof(system)
    colptr=[i for i in 1:ndof+1]
    rowval=[1 for i in 1:ndof]
    nzval=[1.0 for in in 1:ndof]
    jac=SparseMatrixCSC(1,ndof,colptr,rowval,nzval)

    # See https://github.com/JuliaDiff/SparseDiffTools.jl

    # Color the matrix for automtic differentiation
    colors = matrix_colors(jac)

    # Use Julia automatic differentiation for the calculation of the Jacobian
    forwarddiff_color_jacobian!(jac, measurement_functional, values(U0), colorvec = colors)

    # Drop any zero entries 
    dropzeros!(jac)
    
    return jac
end


"""
````
impedance(impedance_system,ω, U0 ,
          excited_spec, excited_bc, excited_bcval,
           dmeas_stdy,
           dmeas_tran 
           )
    
````

Calculate impedance.

-    ω:  frequency 
-    U0: steady state slution
-    dmeas_stdy: Derivative of steady state part of measurement functional
-    dmeas_tran  Derivative of transient part of the measurement functional

"""
function impedance(impedance_system::ImpedanceSystem, # frequency domain system
                   ω,    # frequency 
                   U0 ,  # steady state slution
                   dmeas_stdy, # Derivative of steady state part of measurement functional
                   dmeas_tran  # Derivative of transient part of the measurement functional
                   )
    
    iω=1im*ω
    # solve impedance system
    UZ=unknowns(impedance_system)
    solve!(UZ,impedance_system,ω)
 
    # obtain measurement in frequency  domain
    m_stdy=dmeas_stdy*values(UZ)
    m_tran=dmeas_tran*values(UZ)

    # Calculate complex measurement 
    z=m_stdy[1]+iω*m_tran[1]

    return 1/z
end

"""
$(SIGNATURES)

Calculate reciprocal value of impedance.
-  excited_spec,excited_bc,excited_bcval are ignored.

!!! warning
   This is deprecated: use [`impedance`](@ref).
"""

function freqdomain_impedance(impedance_system::ImpedanceSystem, # frequency domain system
                   ω,    # frequency 
                   U0 ,  # steady state slution
                   excited_spec,excited_bc,excited_bcval,
                   dmeas_stdy, # Derivative of steady state part of measurement functional
                   dmeas_tran  # Derivative of transient part of the measurement functional
                   )
    z=impedance(impedance_system,ω,U0,dmeas_stdy,dmeas_tran)
    1/z
end
