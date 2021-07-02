##########################################################
"""
$(TYPEDEF)
    
Abstract type for finite volume system structure.
"""
abstract type AbstractSystem{Tv<:Number, Ti <:Integer, Tm <:Integer} end

##################################################################
"""
$(TYPEDEF)
    
Structure holding data for finite volume system solution.

Information on species distribution is kept in sparse or dense matrices
matrices and, correspondingly,  the solution array is of type SparseSolutionArray
or matrix, respectively.

In the case of sparse unknown storage, the system matrix handles exactly those
degrees of freedom which correspond to unknowns. However, handling
of the sparse matrix structures for the bookkeeping of the unknowns
creates overhead.

$(TYPEDFIELDS)
"""
mutable struct System{Tv,Ti, Tm, TSpecMat<:AbstractMatrix, TSolArray<:AbstractMatrix} <: AbstractSystem{Tv,Ti, Tm}

    """
    Grid
    """
    grid::ExtendableGrid{Tv,Ti}

    """
    Physics data
    """
    physics::Physics

    """
    Array of boundary values 
    """
    boundary_values::Array{Tv,2} 

    """
    Array of boundary factors 
    """
    boundary_factors::Array{Tv,2}

    """
    Matrix containing species numbers for inner regions
    """
    region_species::TSpecMat

    """
    Matrix containing species numbers for boundary regions
    """
    bregion_species::TSpecMat


    """
    Matrix containing degree of freedom numbers for each node
    """
    node_dof::TSpecMat


    discontspec::TSpecMat

    """
    Jacobi matrix for nonlinear problem
    """
    matrix::ExtendableSparseMatrix{Tv,Tm}

    """
    Matrix factorization
    """
    factorization::Union{Nothing,ExtendableSparse.AbstractFactorization{Tv,Tm}}
    
    """
    Flag which says if the number of unknowns per node is constant
    """
    species_homogeneous::Bool

    """
    Solution vector holding Newton update
    """
    update::TSolArray

    """
    Solution vector holding Newton residual
    """
    residual::TSolArray

    """
    Precomputed geometry factors for cell nodes
    """
    cellnodefactors::Array{Tv,2}
    
    """
    Precomputed geometry factors for cell edges
    """
    celledgefactors::Array{Tv,2}

    """
    Precomputed geometry factors for boundary nodes
    """
    bfacenodefactors::Array{Tv,2}

    """
    Precomputed geometry factors for boundary edges
    """
    bfaceedgefactors::Array{Tv,2}

    """
    Sparse matrix for generic operator handling
    """
    generic_matrix::SparseMatrixCSC

    """
    Sparse matrix colors for generic operator handling
    """
    generic_matrix_colors::Vector

    """
    Hash value of latest unknowns vector the assembly was called with
    """
    uhash::UInt64
    
    """
    Data for allocation check
    """
    allocs::Int
    
    System{Tv,Ti,Tm, TSpecMat, TSolArray}() where {Tv,Ti,Tm, TSpecMat, TSolArray} = new()
end

"""
    const DenseSystem

Type alias for system with dense matrix based species management
"""
const DenseSystem = System{Tv,Ti,Tm,Matrix{Ti},Matrix{Tv}} where {Tv, Ti, Tm}


"""
    const SparseSystem

Type alias for system with sparse matrix based species management
"""
const SparseSystem = System{Tv,Ti,Tm,SparseMatrixCSC{Ti,Ti},SparseSolutionArray{Tv,Ti} } where {Tv, Ti, Tm}


default_check_allocs()= haskey(ENV,"VORONOIFVM_CHECK_ALLOCS") ? parse(Bool,ENV["VORONOIFVM_CHECK_ALLOCS"]) :  false

##################################################################
"""
````
function System(grid,physics;
                unknown_storage=:dense;
                matrixindextype=Int32)
````

Create Finite Volume System. 

- `grid`: 1D/2D/3D discretization grid
- `physics`: Physics struct containing node and edge callbacks
- `unknown_storage`: string or symbol:
     - `:dense` :  solution vector is an  `nspecies` x `nnodes`  dense matrix
     - `:sparse` :  solution vector is an `nspecies` x `nnodes`  sparse matrix
- `matrixindextype` : Index type for sparse matrices created in the system
"""
function System(grid,physics::Physics; unknown_storage=:dense, matrixindextype=Int64, check_allocs=default_check_allocs())

    Tv=coord_type(grid)
    Ti=index_type(grid)
    Tm=matrixindextype

    if Symbol(unknown_storage)==:dense
        system=System{Tv,Ti,Tm,Matrix{Ti}, Matrix{Tv}}()
    elseif Symbol(unknown_storage)==:sparse
        system=System{Tv,Ti,Tm,SparseMatrixCSC{Ti,Ti},SparseSolutionArray{Tv,Ti}}()
    else
        throw("specify either unknown_storage=:dense  or unknown_storage=:sparse")
    end

    maxspec=0
    system.grid=grid
    system.physics=physics
    system.region_species=spzeros(Ti,Int16,maxspec,num_cellregions(grid))
    system.bregion_species=spzeros(Ti,Int16,maxspec,num_bfaceregions(grid))
    system.node_dof=spzeros(Ti,Tm,maxspec,num_nodes(grid))
    system.discontspec=spzeros(Ti,0,num_cellregions(grid))
    system.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    system.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    system.species_homogeneous=false
    system.uhash=0x0
    system.allocs=-1000
    system.factorization=nothing
    check_allocs!(system,check_allocs)
    return system
end

"""
$(SIGNATURES)

Constructor for DenseSystem.
"""
DenseSystem(grid,physics::Physics; matrixindextype=Int64)=System(grid,physics,matrixindextype=matrixindextype,unknown_storage=:dense)


"""
$(SIGNATURES)

Constructor for SparseSystem.
"""
SparseSystem(grid,physics::Physics; matrixindextype=Int64)=System(grid,physics,matrixindextype=matrixindextype,unknown_storage=:sparse)


###################################################################################################
# Test if allocated is zero and if allocation check is enabled
# However we need to be aware that @allocated reports allocations
# during the compilation phase. So we need to wait for at least
# one run of the system in order enact the checking.
function _check_allocs(system, allocated)
       if system.allocs>=0 # we had a couple of runs before to bridge the compilation phase
           system.allocs=allocated
           return system.allocs==0 
       elseif system.allocs > -100 # probably still in compiling phase
           system.allocs=system.allocs+1
           return true
       else
           # otherwise, checking has been switched off.
           return true
       end
end


"""
```
check_allocs!(system,true_or_false)
```

Enable/disable checking for time-consuming allocations in the assembly loop. 
By default,  this check  is switched off.  By setting  the environment
variable `ENV["VORONOIFVM_CHECK_ALLOCS"]="true"`, this  default can be
changed.

Unless  the   matrix  pattern  changes,  there   shouldn't  occur  any
allocations in this loop. The check  method is aware of matrix pattern
changes. As a consequence, allocations in the assembly loop are mostly
due to type instabilities in physics callbacks, see the the discussion
[here](../runexamples/#Performance-with-closures).  Type instabilities
can be debugged via the `@time`  macro applied to expressions in a
physics callback.

The following  cases provide some ideas  where to look for  reasons of
the problem and possible remedies:

Case 1: a parameter changes its value, and Julia is not sure about the type.
```julia
eps=1.0

flux(f,_u,edge)
    u=unkowns(edge,_u)
    f[1]=eps*(u[1,1]-[1,2])
end
... solve etc ...
eps=2.0
```
Remedy: use a type annotation `eps::Float64=...` to signalize your intent to Julia.
This behaviour is explained in the [Julia documentation](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-captured).



Case 2: variables in the closure have the same name as a variable
introduced in a callback.
```julia
flux(f,_u,edge)
    u=unkowns(edge,_u)
    f[1]=(u[1,1]-[1,2])
end

... create etc ...

u=solve(...)
```
Remedy: rename e.g. `u=solve()` to `sol=solve()`

"""
function check_allocs!(system::AbstractSystem,chk::Bool)
    if chk
        system.allocs=-2
    else
        system.allocs=-1000
    end
    system
end

##################################################################

# Constant to be used as boundary condition factor 
# to mark Dirichlet boundary conditons.    
const Dirichlet=1.0e30


#################################################################

"""
    addzrows(matrix,maxrow)

Return matrix with number of rows increased to maxrow, and set
the new elements to zero.
"""
function addzrows(matrix::Matrix,maxrow)
    nrow,ncol=size(matrix)
    if maxrow<=nrow
        return matrix
    end
    newmatrix=zeros(eltype(matrix),maxrow,ncol)
    for icol=1:ncol
        for irow=1:nrow
            newmatrix[irow,icol]=matrix[irow,icol]
        end
    end
    newmatrix
end

function addzrows(matrix::SparseMatrixCSC,maxrow)
    nrow,ncol=size(matrix)
    if maxrow<=nrow
        return matrix
    end
    SparseMatrixCSC(maxrow, matrix.n,matrix.colptr,matrix.rowval,matrix.nzval)
end



"""
     increase_num_species!(system,maxspec)

Increase number of species in system to maxspec by adding new rows to all  relevant
matrices.
"""
function increase_num_species!(system,maxspec)
    system.region_species=addzrows(system.region_species,maxspec)
    system.bregion_species=addzrows(system.bregion_species,maxspec)
    system.node_dof=addzrows(system.node_dof,maxspec)
    system.boundary_values=addzrows(system.boundary_values,maxspec)
    system.boundary_factors=addzrows(system.boundary_factors,maxspec)
end



struct Species
    i::Int
end


function enable_discontinuous_species!(sys,spec::Species,regions)
    sys.discontspec=addzrows(sys.discontspec,spec.i)
    nspec=num_species(sys)
    for ireg ∈ regions
        nspec=nspec+1
        enable_species!(sys,nspec,[ireg])
        sys.discontspec[spec.i,ireg]=nspec
    end
end

function enable_species!(sys,spec::Species,regions)
    sys.discontspec=addzrows(sys.discontspec,spec.i)
    nspec=num_species(sys)
    nspec=nspec+1
    enable_species!(sys,nspec,regions)
    for ireg ∈ regions
        sys.discontspec[spec.i,ireg]=nspec
    end
end


function subgrids(spec::Species, sys)
    grid=sys.grid
    subgrids=Vector{ExtendableGrid}(undef,0)
    for ireg=1:num_cellregions(grid)
        ispec=sys.discontspec[spec.i,ireg]
        if ispec>0
            push!(subgrids,subgrid(grid,[ireg]))
        end
    end
    subgrids
end

function views(U, spec::Species, subgrids,sys)
    grid=sys.grid
    projections=Vector[]
    j=1
    for ireg=1:num_cellregions(grid)
        ispec=sys.discontspec[spec.i,ireg]
        if ispec>0
            push!(projections,view(U[ispec,:],subgrids[j]))
            j=j+1
        end
    end
    projections
end

# # """

# # For a discontinuous quantity, we need to have a different
# # species number for each region.
# # """
# function enable_discontinuous_species(sys,species::Species, regions)
#     spec.regspec=
#     nspec=num_species(sys)
#     for ireg ∈ regions
#         nspec=nspec+1
#         enable_species!(sys,nspec,[ireg])
#         dspec.regspec[ireg]=nspec
#     end
# end

# function continuous_species(sys,species::Species, regions)
#     spec=Species(sys)
#     nspec=num_species(sys)
#     nspec=nspec+1
#     enable_species!(sys,nspec,regions)
#     for ireg ∈ regions
#         dspec.regspec[ireg]=nspec
#     end
# end

# function boundary_species(sys,bregions)
#     spec=Species(sys)
#     nspec=num_species(sys)
#     nspec=nspec+1
#     enable_boundary_species!(sys,nspec,bregions)
#     for ireg ∈ bregions
#         dspec.regspec[ireg]=nspec
#     end
# end


##################################################################
"""
````
is_boundary_species(AbstractSystem, ispec) -> Bool
````

Check if species number corresponds to a boundary species.
"""
function is_boundary_species(system::AbstractSystem, ispec)
    isbspec=false
    for ibreg=1:num_bfaceregions(system.grid)
        if system.bregion_species[ispec,ibreg]>0
            isbspec=true
        end
    end
    return isbspec
end

##################################################################
"""
````
is_bulk_species(AbstractSystem, ispec) -> Bool
````

Check if species number corresponds to a bulk species.
"""
function is_bulk_species(system::AbstractSystem, ispec)
    isrspec=false
    for ixreg=1:num_cellregions(system.grid)
        if system.region_species[ispec,ixreg]>0
            isrspec=true
        end
    end
    return isrspec
end

##################################################################
"""
````
enable_species!(system,ispec,regions)
````

Add species `ispec` to a list of bulk regions. Species numbers for
bulk and boundary species have to be distinct.
"""
function enable_species!(system::AbstractSystem,ispec::Integer, regions::AbstractVector)
    if ispec>num_species(system)
        increase_num_species!(system,ispec)
    end
    
    if is_boundary_species(system,ispec)
        throw(DomainError(ispec,"Species is already boundary species"))
    end
    _cellregions=cellregions(system.grid)
    _cellnodes=cellnodes(system.grid)
    for i in eachindex(regions)
        ireg=regions[i]
        system.region_species[ispec,ireg]=ispec
        for icell=1:num_cells(system.grid)
            if _cellregions[icell]==ireg
                for iloc=1:num_targets(_cellnodes,icell)
                    iglob=_cellnodes[iloc,icell]
                    system.node_dof[ispec,iglob]=ispec
                end
            end
        end
    end
end


##################################################################
"""
````
enable_boundary_species!(system,ispec,regions)
````

Add species `ispec` to a list of boundary regions. Species numbers for
bulk and boundary species have to be distinct.

"""
function enable_boundary_species!(system::AbstractSystem, ispec::Integer, bregions::AbstractVector)
    if ispec>num_species(system)
        increase_num_species!(system,ispec)
    end
    if is_bulk_species(system,ispec)
        throw(DomainError(ispec,"Species is already bulk species"))
    end
    _bfaceregions=bfaceregions(system.grid)
    _bfacenodes=bfacenodes(system.grid)
    
    for i in eachindex(bregions)
        ireg=bregions[i]
        system.bregion_species[ispec,ireg]=ispec
        for ibface=1:num_bfaces(system.grid)
            if _bfaceregions[ibface]==ireg
                for iloc=1:size(_bfacenodes,1)
                    iglob=_bfacenodes[iloc,ibface]
                    system.node_dof[ispec,iglob]=ispec
                end
            end
        end
    end
end


# Create matrix in system and figure out if species
# distribution is homgeneous
function _complete!(system::AbstractSystem{Tv,Ti, Tm};create_newtonvectors=false) where {Tv,Ti, Tm}

    if isdefined(system,:matrix)
        return
    end
    system.matrix=ExtendableSparseMatrix{Tv,Tm}(num_dof(system), num_dof(system))
    system.species_homogeneous=true
    species_added=false
    for inode=1:size(system.node_dof,2)
        for ispec=1:size(system.node_dof,1)
            if system.node_dof[ispec,inode]==ispec
                species_added=true
            else
                system.species_homogeneous=false
            end
        end
    end
    
    if (!species_added)
        error("No species enabled.\n Call enable_species(system,species_number, list_of_regions) at least once.")
    end
    if create_newtonvectors
        system.residual=unknowns(system)
        system.update=unknowns(system)
    end
    update_grid!(system)
    if has_generic_operator(system)
        if has_generic_operator_sparsity(system) 
            system.generic_matrix=system.physics.generic_operator_sparsity(system)
        else
            generic_operator(f,u)=system.physics.generic_operator(f,u,system)
            input=rand(num_dof(system))
            output=similar(input)
            tdetect=@elapsed begin
                sparsity_pattern = jacobian_sparsity(generic_operator,output,input)
                system.generic_matrix = Float64.(sparse(sparsity_pattern))
            end
            println("sparsity detection for generic operator: $(tdetect) s")
            if nnz(system.generic_matrix)==0
                error("Sparsity detection failed: no pattern found")
            end
        end
        tdetect=@elapsed begin
            system.generic_matrix_colors = matrix_colors(system.generic_matrix)
        end
        println("matrix coloring for generic operator: $(tdetect) s")
    end
end

"""
$(SIGNATURES)
Set generic operator sparsity, in the case where a generic operator has been
defined in physics.
"""
function generic_operator_sparsity!(system::AbstractSystem, sparsematrix::SparseMatrixCSC)
    system.generic_matrix=sparsematrix
end


"""
````
update_grid!(system; grid=system.grid)
````

Update grid (e.g. after rescaling of coordinates).
"""
function update_grid!(system::AbstractSystem{Tv,Ti,Tm}; grid=system.grid) where {Tv,Ti,Tm}
    
    geom        = grid[CellGeometries][1]
    csys        = grid[CoordinateSystem]
    coord       = grid[Coordinates]
    cellnodes   = grid[CellNodes]
    cellregions = grid[CellRegions]
    bgeom       = grid[BFaceGeometries][1]
    bfacenodes  = grid[BFaceNodes]
    nbfaces     = num_bfaces(grid)
    ncells      = num_cells(grid)

    system.cellnodefactors  = zeros(Tv, num_nodes(geom), ncells)
    system.celledgefactors  = zeros(Tv, num_edges(geom), ncells)
    system.bfacenodefactors = zeros(Tv, num_nodes(bgeom), nbfaces)
    system.bfaceedgefactors = zeros(Tv, num_edges(bgeom), nbfaces)
    
    for icell=1:ncells
        @views cellfactors!(geom, csys, coord, cellnodes, icell,
                            system.cellnodefactors[:,icell],system.celledgefactors[:,icell])
    end
    
    for ibface=1:nbfaces
        @views bfacefactors!(bgeom,csys,coord,bfacenodes,ibface,
                             system.bfacenodefactors[:,ibface], system.bfaceedgefactors[:, ibface])
    end
end

##################################################################
"""
$(SIGNATURES)
    
Check if degree of freedom is defined.
"""
isdof(system::AbstractSystem,ispec,inode)= system.node_dof[ispec,inode]==ispec


##################################################################
"""
$(SIGNATURES)

Number of species in system
"""
num_species(system::AbstractSystem) = size(system.node_dof,1)



##################################################################
"""
$(SIGNATURES)

Retrieve user data record.
"""
data(system::AbstractSystem) = system.physics.data


##################################################################
"""
$(SIGNATURES)

Set Dirichlet boundary condition for species ispec at boundary ibc:

``u_{ispec}=v`` on ``\\Gamma_{ibc}``
"""
function boundary_dirichlet!(system::AbstractSystem, ispec, ibc, v)
    system.boundary_factors[ispec,ibc]=Dirichlet
    system.boundary_values[ispec,ibc]=v
end


# just return the first which comes into mind.
# we need to ensure homgeneity of bregion-region structure.
# if that is true, this works.
function get_cellregion(sys,ibc)
    grid=sys.grid
    bfregions=grid[BFaceRegions]
    cregions=grid[CellRegions]
    for ibface=1:num_bfaces(sys.grid)
        if bfregions[ibface]==ibc
            bfcells=grid[BFaceCells]
            return cregions[bfcells[ibface,1]]
        end
    end
    return 0
end

boundary_dirichlet!(sys::AbstractSystem, spec::Species, ibc, v)=boundary_dirichlet!(sys,
                                                                                   sys.discontspec[spec.i,get_cellregion(sys,ibc)],
                                                                                   ibc,
                                                                                   v)

##################################################################
"""
$(SIGNATURES)

Set Neumann boundary condition for species ispec at boundary ibc:

``\\mathrm{flux}_{ispec}\\cdot \\vec n=v`` on ``\\Gamma_{ibc}``
"""
function boundary_neumann!(system::AbstractSystem, ispec, ibc, v)
    system.boundary_factors[ispec,ibc]=0.0
    system.boundary_values[ispec,ibc]=v
end


##################################################################
"""
$(SIGNATURES)

Set Robin boundary condition for species ispec at boundary ibc:

``\\mathrm{flux}_{ispec}\\cdot \\vec n + \\alpha u_{ispec}=v`` on ``\\Gamma_{ibc}``

"""
function boundary_robin!(system::AbstractSystem, ispec, ibc, α, v)
    system.boundary_factors[ispec,ibc]=α
    system.boundary_values[ispec,ibc]=v
end

##################################################################
"""
$(SIGNATURES)

Number of species (size of first dimension) of solution array.
"""
num_species(a::AbstractArray)=size(a,1)







#
# Initialize Dirichlet BC
#
function _initialize_dirichlet!(U::AbstractMatrix,system::AbstractSystem)
    _bfaceregions=bfaceregions(system.grid)
    _bfacenodes=bfacenodes(system.grid)
    for ibface=1:num_bfaces(system.grid)
        ibreg=_bfaceregions[ibface]
        for ispec=1:num_species(system)
            if system.boundary_factors[ispec,ibreg]≈ Dirichlet
                for inode=1:dim_grid(system.grid)
                    U[ispec,_bfacenodes[inode,ibface]]=system.boundary_values[ispec,ibreg]
                end
            end
        end
    end
end




function _initialize!(U::AbstractMatrix,system::AbstractSystem)
    _initialize_dirichlet!(U,system)
    _initialize_inactive_dof!(U,system)
end



function _eval_and_assemble_inactive_species(system::AbstractSystem,U,Uold,F) end

function _eval_and_assemble_inactive_species(system::DenseSystem,U,Uold,F)
    if system.species_homogeneous
        return
    end
    for inode=1:size(system.node_dof,2)
        for ispec=1:size(system.node_dof,1)
            if !isdof(system,ispec,inode)
                F[ispec,inode]+= U[ispec,inode]-Uold[ispec,inode];
                idof=dof(F,ispec,inode)
                system.matrix[idof,idof]+=1.0
            end
        end
    end
end

function _initialize_inactive_dof!(U::AbstractMatrix,system::AbstractSystem) end

function _initialize_inactive_dof!(U::DenseSolutionArray,system::DenseSystem)
    if system.species_homogeneous
        return
    end
    for inode=1:size(system.node_dof,2)
        for ispec=1:size(system.node_dof,1)
            if !isdof(system,ispec,inode)
                U[ispec,inode]=0
            end
        end
    end
end


function Base.show(io::IO,sys::AbstractSystem)
    str=@sprintf("%s(num_species=%d)",typeof(sys),num_species(sys))
    println(io,str)
end

#####################################################
has_generic_operator(sys::AbstractSystem) = sys.physics.generic_operator!=nofunc_generic
has_generic_operator_sparsity(sys::AbstractSystem) =  sys.physics.generic_operator_sparsity!=nofunc_generic_sparsity



##################################################################
"""
$(SIGNATURES)

Number of degrees of freedom for system.
"""
function num_dof(system::AbstractSystem) end

num_dof(system::SparseSystem)= nnz(system.node_dof)

num_dof(system::DenseSystem)= length(system.node_dof)


"""
$(SIGNATURES)

Create a solution vector for system.
If inival is not specified, the entries of the returned vector are undefined.
"""
function unknowns(system::AbstractSystem;inival=undef) end

unknowns(sys::SparseSystem{Tv,Ti,Tm};inival=undef) where {Tv,Ti, Tm}=unknowns(Tv,sys,inival=inival)

unknowns(system::DenseSystem{Tv,Ti,Tm};inival=undef) where {Tv,Ti, Tm} = unknowns(Tv,system,inival=inival)


"""
$(SIGNATURES)

Create a solution vector for system with elements of type `Tu`.
If inival is not specified, the entries of the returned vector are undefined.
"""
function unknowns(Tu::Type, sys::AbstractSystem; inival=undef) end

function unknowns(Tu::Type, system::SparseSystem;inival=undef)
    a0=Array{Tu}(undef,num_dof(system))
    if inival!=undef
        fill!(a0,inival)
    end
    return SparseSolutionArray(SparseMatrixCSC(system.node_dof.m,
                                               system.node_dof.n,
                                               system.node_dof.colptr,
                                               system.node_dof.rowval,
                                               a0
                                               )
    )
end

function unknowns(Tu::Type, sys::DenseSystem; inival=undef)
    a=Array{Tu}(undef,size(sys.node_dof)...)
    if inival!=undef
        fill!(a,inival)
    end
    return a
end





"""
$(SIGNATURES)

Reshape vector to fit as solution to system.
"""
function Base.reshape(v,system::AbstractSystem) end

Base.reshape(v::DenseSolutionArray,system::DenseSystem)=v

Base.reshape(v::SparseSolutionArray,sys::SparseSystem)=v

function Base.reshape(v::AbstractVector, sys::DenseSystem)
    @assert  length(v)==num_dof(sys)
    nspec=num_species(sys)
    reshape(v,Int64(nspec),Int64(length(v)/nspec))
end

function Base.reshape(v::AbstractVector,system::SparseSystem)
    @assert  length(v)==num_dof(system)
    SparseSolutionArray(SparseMatrixCSC(system.node_dof.m,
                                        system.node_dof.n,
                                        system.node_dof.colptr,
                                        system.node_dof.rowval,
                                        Vector(v)
                                        )
                        )
end



"""
$(SIGNATURES)

Calculate norm, paying attention to species distribution over regions
"""
function LinearAlgebra.norm(system::AbstractSystem,u,p) end

function LinearAlgebra.norm(system::DenseSystem,u,p)
    _initialize_inactive_dof!(u,system)
    norm(u,p)
end

LinearAlgebra.norm(system::SparseSystem,u,p)=norm(u.node_dof.nzval,p)



