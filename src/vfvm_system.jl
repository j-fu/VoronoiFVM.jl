"""
$(TYPEDEF)

Structure holding data for finite volume system.
"""
mutable struct System{Tv, Tc, Ti, Tm, TSpecMat<:AbstractMatrix, TSolArray<:AbstractMatrix} <: AbstractSystem{Tv, Tc, Ti, Tm}

    """
    Grid
    """
    grid::ExtendableGrid{Tc,Ti}

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


    
    """
    - :multidiagonal  (currently disabled)
    - :sparse
    - :banded
    - :tridiagonal
    """
    matrixtype::Symbol

    
    """
    Jacobi matrix for nonlinear problem
    """
    matrix::Union{ExtendableSparseMatrix{Tv,Tm},
                  Tridiagonal{Tv, Vector{Tv}},
#                  MultidiagonalMatrix,
                  BandedMatrix{Tv}}

    """
    Matrix factorization
    """
    factorization::Union{Nothing,ExtendableSparse.AbstractFactorization{Tv,Tm}}
    
    """
    Flag which says if the number of unknowns per node is constant
    """
    species_homogeneous::Bool

    """
    Number of quantities defined on system
    """
    num_quantities::Ti


    """
    Number of parameter the system depnds on.
    """
    num_parameters::Ti

    """
    Parameter derivative (vector of solution arrays)
    """
    dudp::Vector{TSolArray}
    
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


    """
    History record for last solution process
    """
    history
    
    System{Tv,Tc,Ti,Tm, TSpecMat, TSolArray}() where {Tv, Tc, Ti,Tm, TSpecMat, TSolArray} = new()
end

"""
    const DenseSystem

Type alias for system with dense matrix based species management

"""
const DenseSystem = System{Tv,Tc, Ti,Tm,Matrix{Ti},Matrix{Tv}} where {Tv, Tc, Ti, Tm}


isdensesystem(s::System{Tv, Tc, Ti, Tm, TSpecMat, TSolArray}) where {Tv, Tc, Ti, Tm, TSpecMat, TSolArray} = TSolArray<:Matrix


"""
    const SparseSystem

Type alias for system with sparse matrix based species management

"""
const SparseSystem = System{Tv,Tc,Ti,Tm,SparseMatrixCSC{Ti,Ti},SparseSolutionArray{Tv,Ti} } where {Tv,Tc, Ti, Tm}

default_check_allocs()= haskey(ENV,"VORONOIFVM_CHECK_ALLOCS") ? parse(Bool,ENV["VORONOIFVM_CHECK_ALLOCS"]) :  false

##################################################################
"""
````
System(grid; kwargs...)
````
    
Create structure of type [`VoronoiFVM.System{Tv,Ti, Tm, TSpecMat<:AbstractMatrix, TSolArray<:AbstractMatrix}`](@ref)  holding data for finite volume system solution. 

Parameters: 

- `grid::ExtendableGrid`: 1, 2 or 3D computational grid

Keyword arguments:
- `species`: vector of integer species indices. Added to all grid regions, avoiding the need to call [`enable_species!`](@ref) for this default case.
              If it is kept empty, species have be added to the system after creation via  [`enable_species!`](@ref).
- `unknown_storage`: string or symbol.  
    Information  on  species  distribution  is kept  in  sparse  or  dense
    matrices matrices and, correspondingly, the  solution array is of type
    SparseSolutionArray  or matrix,  respectively. In  the case  of sparse
    unknown storage,  the system matrix  handles exactly those  degrees of
    freedom which correspond to unknowns.  However, handling of the sparse
    matrix  structures  for  the   bookkeeping  of  the  unknowns  creates
    overhead.
     - `:dense` :  solution vector is an  `nspecies` x `nnodes`  dense matrix
     - `:sparse` :  solution vector is an `nspecies` x `nnodes`  sparse matrix
- `matrixindextype`: Integer type. Index type for sparse matrices created in the system.

Physics keyword arguments:
- `flux`: Function.     Flux between neigboring control volumes: `flux(f,u,edge)` or `flux(f,u,edge,data)`
    should return in `f[i]` the flux of species i along the edge joining circumcenters
    of neigboring control volumes.  For species i,`u[i,1]` and `u[i,2]` contain the unknown values at the corresponding ends of the edge.
- `storage`: Function.  Storage term (term under time derivative): `storage(f,u,node)` or `storage(f,u,node,data)` 
    It should return in `f[i]` the storage term for the i-th equation. `u[i]` contains the value of
    the i-th unknown.
- `reaction`:  Function. Reaction term:  `reaction(f,u,node)` or `reaction(f,u,node,data)` 
    It should return in `f[i]` the reaction term for the i-th equation. `u[i]` contains the value of
    the i-th unknown.
- `source`:  Function. Source term: `source(f,node)` or `source(f,node,data)`.
    It should return the in `f[i]` the value of the source term for the i-th equation.
- `bflux`:  Function. Flux between neighboring control volumes on the boundary
- `breaction` Function.  Boundary reaction term:  `breaction(f,u,node)` or `breaction(f,u,node,data)` 
    Similar to reaction, but restricted to the inner or outer boundaries.
- `bcondition` Function. Alias for `breaction`.
- `bsource`: Function. Boundary source term: `bsource(f,node)` or `bsource(f,node,data)`.
    It should return in `f[i]` the value of the source term for the i-th equation.
- `bstorage`: Function.  Boundary storage term: `bstorage(f,u,node)` or `bstorage(f,u,node,data)` 
    Similar to storage, but restricted to the inner or outer boundaries.
- `generic_operator`: Funtion.  Generic operator  `generic_operator(f,u,sys)`. 
    This operator acts on the full solution `u` of a system. Sparsity
    is detected automatically  unless `generic_operator_sparsity` is given.
-  `generic_operator_sparsity`:  Function defining the sparsity structure of the generic operator.
    This should return the sparsity pattern of the `generic_operator`.
-  `nparams`: number of parameters the system is depending on, and with respect to which the derivatives
    need to be obtained
-  `data`:  User data (parameters).
    This allows to pass various parameters to the callback functions. If `data` is given, all callback functions
    should accept a last `data` argument. Otherwise, no data are passed explicitely, and constitutive callbacks can
    take parameters from the closure where the function is defined.

-  `matrixtype`: :default, :sparse, :tridiagonal, :banded

"""
function System(grid::ExtendableGrid;
                valuetype=coord_type(grid),
                indextype=index_type(grid),
                species=Int[],
                unknown_storage=:dense,
                matrixindextype=Int64,
                matrixtype=:sparse,
                check_allocs=default_check_allocs(),
                nparams=0,
                kwargs...)

    Tv=valuetype
    Tc=coord_type(grid)
    Ti=indextype
    Tm=matrixindextype
    
    
    
    if Symbol(unknown_storage)==:dense
        system=System{Tv,Tc,Ti,Tm,Matrix{Ti}, Matrix{Tv}}()
    elseif Symbol(unknown_storage)==:sparse
        system=System{Tv,Tc,Ti,Tm,SparseMatrixCSC{Ti,Ti},SparseSolutionArray{Tv,Ti}}()
    else
        throw("specify either unknown_storage=:dense  or unknown_storage=:sparse")
    end

    
    maxspec=0
    system.grid=grid
    system.region_species=spzeros(Ti,Int16,maxspec,num_cellregions(grid))
    system.bregion_species=spzeros(Ti,Int16,maxspec,num_bfaceregions(grid))
    system.node_dof=spzeros(Ti,Tm,maxspec,num_nodes(grid))
    system.boundary_values=zeros(Tv,maxspec,num_bfaceregions(grid))
    system.boundary_factors=zeros(Tv,maxspec,num_bfaceregions(grid))
    system.species_homogeneous=false
    system.num_quantities=0
    system.uhash=0x0
    system.matrixtype=matrixtype
    system.allocs=-1000
    system.factorization=nothing
    system.history=nothing
    system.num_parameters=nparams
    
    check_allocs!(system,check_allocs)

    physics!(system; kwargs...)
    enable_species!(system;species)
    return system
end


"""
````
System(X; kwargs...)
````
Create an [1D grid from vector X](https://j-fu.github.io/ExtendableGrids.jl/stable/gridconstructors/) and call  [`VoronoiFVM.System(grid::ExtendableGrid; kwargs...)`](@ref).
"""
System(X::AbstractVector; kwargs...)= System(simplexgrid(X); kwargs...)


"""
````
System(X,Y; kwargs...)
````
Create a [2D grid from vectors X,Y ](https://j-fu.github.io/ExtendableGrids.jl/stable/gridconstructors/) and call  [`VoronoiFVM.System(grid::ExtendableGrid; kwargs...)`](@ref).
"""
System(X::AbstractVector, Y::AbstractVector; kwargs...)= System(simplexgrid(X,Y); kwargs...)


"""
````
System(X,Y, Z; kwargs...)
````
Create a [3D grid from vectors X,Y,Z ](https://j-fu.github.io/ExtendableGrids.jl/stable/gridconstructors/) and call  [`VoronoiFVM.System(grid::ExtendableGrid; kwargs...)`](@ref).
"""
System(X::AbstractVector, Y::AbstractVector, Z::AbstractVector; kwargs...)= System(simplexgrid(X,Y,Z); kwargs...)


"""
    physics!(system,physics)

Replace System's physics data
"""
function physics!(system,physics)
    system.physics=physics
    system
end

"""
   physics!(system; kwargs...)

Replace System's physics data.
"""
function physics!(system; kwargs...)
    kwdict=Dict(kwargs)

    if haskey(kwdict,:bcondition)
        if haskey(kwdict,:breaction)
            error("specify either bcondition or breaction")
        end
        kwdict[:breaction]=kwdict[:bcondition]
        delete!(kwdict,:bcondition)
    end

    if haskey(kwdict,:source) && isa(kwdict[:source], AbstractArray)
        src=kwdict[:source]
        if isa(src,AbstractVector)
            kwdict[:source]= (y,node,args...) -> @views y[1]=src[node.index]
        else
            kwdict[:source]= (y,node,args...) -> @views y.=src[:,node.index]
        end
    end
    
    physics!(system,Physics(; kwdict...))
end


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
    if maxspec<=num_species(system)
        return
    end
    
    if isdefined(system,:matrix)
        error("Unable to increase number of species to $(maxspec).\nPlease add species before first solver run.")
    end
    
    system.region_species=addzrows(system.region_species,maxspec)
    system.bregion_species=addzrows(system.bregion_species,maxspec)
    system.node_dof=addzrows(system.node_dof,maxspec)
    system.boundary_values=addzrows(system.boundary_values,maxspec)
    system.boundary_factors=addzrows(system.boundary_factors,maxspec)
end





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
Once a species has been added, it cannot be removed.
"""
function enable_species!(system::AbstractSystem,ispec::Integer, regions::AbstractVector)
    increase_num_species!(system,ispec)
    
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


"""
````
enable_species!(system; kwargs...)
````

Keyword arguments:
- `species`: Integer or vector of integers. Species to be added to the system.
- `regions`: Vector of integers. Regions, where these species shall be added.If `nothing`, they are added to all species.
Once a species has been added, it cannot be removed.
"""
function enable_species!(sys::AbstractSystem; species=nothing, regions=nothing)
    if regions==nothing
        regions=collect(1:num_cellregions(sys.grid))
    end
    
    if isa(species,Number)
        species=[species]
    end
    
    for ispec ∈ species
        enable_species!(sys,ispec,regions)
    end
    sys
end


##################################################################
"""
````
enable_boundary_species!(system,ispec,regions)
````

Add species `ispec` to a list of boundary regions. Species numbers for
bulk and boundary species have to be distinct.
Once a species has been added, it cannot be removed.
"""
function enable_boundary_species!(system::AbstractSystem, ispec::Integer, bregions::AbstractVector)
    increase_num_species!(system,ispec)

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
function _complete!(system::AbstractSystem{Tv,Tc,Ti, Tm};create_newtonvectors=false) where {Tv,Tc,Ti, Tm}

    if isdefined(system,:matrix)
        return
    end


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
    

    
    nspec=size(system.node_dof,1)
    n=num_dof(system)
    
    matrixtype=system.matrixtype
    #    matrixtype=:sparse
    # Sparse even in 1D is not bad, 
    
    if matrixtype==:default
        if !isdensesystem(system)
            matrixtype=:sparse
        else
            if nspec==1
                matrixtype=:tridiagonal
            else
                matrixtype=:banded
            end                
        end
    end

    if matrixtype==:tridiagonal
        system.matrix=Tridiagonal(zeros(Tv,n-1),zeros(Tv,n),zeros(Tv,n-1))
    elseif matrixtype==:banded
        system.matrix=BandedMatrix{Tv}(Zeros(n,n), (2*nspec-1,2*nspec-1))
    # elseif matrixtype==:multidiagonal
    #     system.matrix=mdzeros(Tv,n,n,[-1,0,1]; blocksize=nspec)
    else # :sparse
        system.matrix=ExtendableSparseMatrix{Tv,Tm}(n,n)
    end

    
    if create_newtonvectors
        system.residual=unknowns(system)
        system.update=unknowns(system)
        if system.num_parameters>0
            system.dudp=[unknowns(system) for i=1:system.num_parameters]
        end
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
                sparsity_pattern = Symbolics.jacobian_sparsity(generic_operator,output,input)
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
function update_grid!(system::AbstractSystem{Tv,Tc,Ti,Tm}; grid=system.grid) where {Tv,Tc,Ti,Tm}
    
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

    function barrier(csys::Type{T}) where T
        nalloc=@allocated for icell=1:ncells
            @views cellfactors!(geom, csys, coord, cellnodes, icell,
                                system.cellnodefactors[:,icell],system.celledgefactors[:,icell])
        end
        nalloc>0 &&  @warn "$nalloc allocations in cell factor calculation"
            
        nalloc=@allocated  for ibface=1:nbfaces
            @views bfacefactors!(bgeom,csys,coord,bfacenodes,ibface,
                                 system.bfacenodefactors[:,ibface], system.bfaceedgefactors[:, ibface])
        end
        nalloc>0 &&  @warn "$nalloc allocations in bface factor calculation"

    end
    barrier(csys)
end

##################################################################
"""
$(SIGNATURES)
    
Check if degree of freedom is defined.
"""
isdof(system::AbstractSystem,ispec,inode)= system.node_dof[ispec,inode]==ispec



#
# DOF ASSEMBLY API
#
_firstnodedof(sys::DenseSystem{Tv},K::Integer) where Tv = (K-1)*num_species(sys)+1
_lastnodedof(sys::DenseSystem{Tv},K::Integer) where Tv = K*num_species(sys)
_species_of_dof(sys::DenseSystem{Tv},idof,K) where Tv =   idof-(K-1)*num_species(sys)

_firstnodedof(sys::SparseSystem, K::Integer) =sys.node_dof.colptr[K]
_lastnodedof(sys::SparseSystem, K::Integer) =sys.node_dof.colptr[K+1]-1
_species_of_dof(sys::SparseSystem,idof,K) =sys.node_dof.rowval[idof]


"""
$(SIGNATURES)

Assemble residual and jacobian for node functions. Parameters:

- `system`: System to be worked with
- `node`: node
- `asm_res(idof,ispec)`: e.g. assemble local ispec to global degree of freedom in unknowns
- `asm_jac(idof,jdof,ispec,jspec)`: e.g.  assemble entry `ispec,jspec` of local jacobian into entry `idof,jdof` of global matrix
- `asm_param(idof,ispec,iparam)` shall assemble parameter derivatives
"""
function assemble_res_jac(node::Node,system::AbstractSystem,asm_res::R,asm_jac::J, asm_param::P) where {R,J,P}
    K=node.index
    ireg=node.region
    for idof=_firstnodedof(system,K):_lastnodedof(system,K)
        ispec=_species_of_dof(system,idof,K)
        if system.region_species[ispec,ireg]>0 # it is not enough to know if the species are defined...
            asm_res(idof,ispec)
            for jdof=_firstnodedof(system,K):_lastnodedof(system,K)
                jspec=_species_of_dof(system,jdof,K)
                if system.region_species[jspec,ireg]>0
                    asm_jac(idof,jdof,ispec,jspec)
                end
            end
        end
        for iparam=1:system.num_parameters
            asm_param(idof,ispec,iparam) 
        end
    end
end

"""
$(SIGNATURES)

Assemble residual and jacobian for boundary node functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
function assemble_res_jac(bnode::BNode,system::AbstractSystem, asm_res::R,asm_jac::J, asm_param::P) where {R,J,P}
    F=system.residual
    K=bnode.index
    ireg=bnode.region
    for idof=_firstnodedof(system,K):_lastnodedof(system,K)
        ispec=_species_of_dof(system,idof,K)
        if isdof(system,ispec,K)
            asm_res(idof,ispec)
            for jdof=_firstnodedof(system,K):_lastnodedof(system,K)
                jspec=_species_of_dof(system,jdof,K)
                if isdof(system,jspec,K)
                    asm_jac(idof,jdof,ispec,jspec)
                end
            end
        end
        for iparam=1:system.num_parameters
            asm_param(idof,ispec,iparam) 
        end
    end
end


"""
$(SIGNATURES)

Assemble residual for node functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
function assemble_res(node::Node, system::AbstractSystem, asm_res::R) where {R}
    F=system.residual
    K=node.index
    ireg=node.region
    for idof=_firstnodedof(system,K):_lastnodedof(system,K)
        ispec=_species_of_dof(system,idof,K)
        if system.region_species[ispec,ireg]>0
            asm_res(idof,ispec)
        end
    end
end

"""
$(SIGNATURES)

Assemble residual for boundary node functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
function assemble_res(bnode::BNode, system::AbstractSystem, asm_res::R) where {R}
    F=system.residual
    K=bnode.index
    ireg=node.region
    for idof=_firstnodedof(system,K):_lastnodedof(system,K)
        ispec=_species_of_dof(system,idof,K)
        if isdof(system,ispec,K)
            asm_res(idof,ispec)
        end
    end
end


"""
$(SIGNATURES)

Assemble residual and jacobian for edge (flux) functions. Parameters:

- `system`: System to be worked with
- `node`: node
- `asm_res(idofK,idofL,ispec)`: e.g. assemble local ispec to global degrees of freedom in unknowns
- `asm_jac(idofK,jdofK,idofL,jdofL,ispec,jspec)`: e.g.  assemble entry `ispec,jspec` of local jacobian into entry four entries defined by `idofK` and `idofL` of global matrix
- `asm_param(idofK,idofL,ispec,iparam)` shall assemble parameter derivatives
"""
function assemble_res_jac(edge::Edge,system::AbstractSystem, asm_res::R,asm_jac::J, asm_param::P ) where {R,J,P}
    F=system.residual
    K=edge.node[1]
    L=edge.node[2]
    ireg=edge.region
    
    for idofK=_firstnodedof(system,K):_lastnodedof(system,K)
        ispec=_species_of_dof(system,idofK,K)
        idofL=dof(F,ispec,L)
        if idofL==0
            continue
        end
        if system.region_species[ispec,ireg]<=0
            continue
        end

        asm_res(idofK,idofL,ispec)
                
        for jdofK=_firstnodedof(system,K):_lastnodedof(system,K)
            jspec=_species_of_dof(system,jdofK,K)
            if system.region_species[jspec,ireg]<=0
                continue
            end
            jdofL=dof(F,jspec,L)
            if jdofL==0
                continue
            end

            asm_jac(idofK,jdofK,idofL,jdofL,ispec,jspec)
        end
        
        for iparam=1:system.num_parameters
            asm_param(idofK,idofL,ispec,iparam)            
        end
    end

end


"""
$(SIGNATURES)

Assemble residual for edge (flux) functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
function assemble_res(edge::Edge,system::AbstractSystem, asm_res::R) where {R}
    F=system.residual
    K=edge.node[1]
    L=edge.node[2]
    ireg=edge.region
    
    for idofK=_firstnodedof(system,K):_lastnodedof(system,K)
        ispec=_species_of_dof(system,idofK,K)
        idofL=dof(F,ispec,L)
        if idofL==0
            continue
        end
        if system.region_species[ispec,ireg]<=0
            continue
        end
        asm_res(idofK,idofL,ispec)
    end
end



"""
$(SIGNATURES)

Assemble residual and jacobian for boundary edge (flux) functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
function assemble_res_jac(bedge::BEdge,system::AbstractSystem,asm_res::R,asm_jac::J, asm_param::P ) where {R,J,P}
    F=system.residual
    K   = bedge.node[1]
    L   = bedge.node[2]
    
    
    for idofK = _firstnodedof(system, K):_lastnodedof(system, K)
        ispec =_species_of_dof(system, idofK, K)
        if !isdof(system,ispec,K)
            continue
        end
        
        idofL = dof(F, ispec, L)
        if idofL == 0
            continue
        end
        
        asm_res(idofK,idofL,ispec)
        
        for jdofK = _firstnodedof(system,K):_lastnodedof(system,K)
            jspec = _species_of_dof(system,jdofK,K)
            if !isdof(system,jspec,K)
                continue
            end
            
            jdofL = dof(F,jspec,L)
            if jdofL == 0
                continue
            end
            asm_jac(idofK,jdofK,idofL,jdofL,ispec,jspec)
        end
    end
    
    for iparam=1:system.num_parameters
        asm_param(idofK,idofL,ispec,iparam)            
    end
end
        



    
"""
$(SIGNATURES)

Assemble residual for boundary edge (flux) functions.
See [`assemble_res_jac`](@ref) for more explanations.
"""
function assemble_res(bedge::BEdge,system::AbstractSystem,asm_res::R) where {R}
    F=system.residual
    K   = bedge.node[1]
    L   = bedge.node[2]
    for idofK = _firstnodedof(system, K):_lastnodedof(system, K)
        ispec =_species_of_dof(system, idofK, K)
        if !isdof(system,ispec,K)
            continue
        end
        
        idofL = dof(F, ispec, L)
        if idofL == 0
            continue
        end

        asm_res(idofK,idofL,ispec)
    end
end
    


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

!!! info  
    Starting with version 0.14, it is preferable to define boundary condtitions within the `bcondition` physics callback
"""
function boundary_dirichlet!(system::AbstractSystem, ispec, ibc, v)
    increase_num_species!(system,ispec)
    system.boundary_factors[ispec,ibc]=Dirichlet
    system.boundary_values[ispec,ibc]=v
end

"""
      boundary_dirichlet!(system; kwargs...)
Keyword argument version:
- `species`: species number
- `region`: region number
- `value`: value

!!! info  
    Starting with version 0.14, it is preferable to define boundary condtions within the `bcondition` physics callback
"""
boundary_dirichlet!(system::AbstractSystem; species=1, region=1,value=0)=boundary_dirichlet!(system,species,region,value)


"""
     boundary_dirichlet!(y,u,bnode,ispec,ireg,val)

Set Dirichlet boundary condition for species ispec at boundary ibc.
"""
function boundary_dirichlet!(y,u,bnode,ispec,ireg,val; penalty=bnode.Dirichlet)
    if  bnode.region == ireg
        y[ispec] += penalty*(u[ispec]-val)
        # just for call during initialization, so we can convert from dual number
        bnode.dirichlet_value[ispec]=value(val)
    end
    nothing
end

"""
     boundary_dirichlet!(y,u,bnode, args...; kwargs...)
Keyword argument version:
- `species`: species number. Default: 1
- `region`: boundary region number. By default, all boundary regions.
- `value`: value
"""
boundary_dirichlet!(y,u,bnode,args...; species=1, region=bnode.region,value=0,penalty=bnode.Dirichlet) = boundary_dirichlet!(y,u,bnode,species,region,value;penalty)


"""
       ramp(t; kwargs...)
Ramp function for specifying time dependent boundary conditions

Keyword arguments:
- `dt`: Tuple: start and end time of ramp. Default: `(0,0.1)`
- `du`: Tuple: values at start and end time. Default: `(0,0)`
"""
function ramp(t;dt=(0,0.1),du=(0,0))
    (t,ubegin,uend,tbegin,tend)=promote(Float64(t),du[1],du[2],dt[1],dt[2])
    if t<tbegin
        return ubegin
    elseif t<tend
        return ubegin+(uend-ubegin)*(t-tbegin)/(tend-tbegin)
    else
        return uend
    end
end


##################################################################
"""
$(SIGNATURES)

Set Neumann boundary condition for species ispec at boundary ibc:

``\\mathrm{flux}_{ispec}\\cdot \\vec n=v`` on ``\\Gamma_{ibc}``
!!! info  
    Starting with version 0.14, it is preferable to define boundary condtitions within the `bcondition` physics callback
"""
function boundary_neumann!(system::AbstractSystem, ispec, ibc, v)
    increase_num_species!(system,ispec)
    system.boundary_factors[ispec,ibc]=0.0
    system.boundary_values[ispec,ibc]=v
end

"""
      boundary_neumann!(system; kwargs...)
Keyword argument version:
- `species`: species number
- `region`: region number
- `value`: value
!!! info  
    Starting with version 0.14, it is preferable to define boundary condtitions within the `bcondition` physics callback
"""
boundary_neumann!(system::AbstractSystem; species=0, region=0,value=0)=boundary_neumann!(system,species,region,value)


"""
     boundary_neumann!(y,u,bnode,ispec,ireg,val)

Set Neumann boundary condition for species ispec at boundary ibc.
"""
boundary_neumann!(y,u,bnode,ispec,ireg,val) =    bnode.region == ireg ? y[ispec] -= val : nothing

"""
     boundary_neumann!(y,u,bnode, args...; kwargs...)
Keyword argument version:
- `species`: species number. Default: 1
- `region`: boundary region number. By default, all boundary regions.
- `value`: value
"""
boundary_neumann!(y,u,bnode,args...; species=1, region=bnode.region,value=0) = boundary_neumann!(y,u,bnode,species,region,value)


##################################################################
"""
$(SIGNATURES)

Set Robin boundary condition for species ispec at boundary ibc:

``\\mathrm{flux}_{ispec}\\cdot \\vec n + \\alpha u_{ispec}=v`` on ``\\Gamma_{ibc}``

!!! info  
    Starting with version 0.14, it is preferable to define boundary condtitions within the `bcondition` physics callback
"""
function boundary_robin!(system::AbstractSystem, ispec, ibc, α, v)
    increase_num_species!(system,ispec)
    system.boundary_factors[ispec,ibc]=α
    system.boundary_values[ispec,ibc]=v
end

"""
      boundary_robin!(system; kwargs...)
Keyword argument version:
- `species`: species number
- `region`: region number
- `factor`: factor
- `value`: value
!!! info  
    Starting with version 0.14, it is preferable to define boundary condtitions within the `bcondition` physics callback
"""
boundary_robin!(system::AbstractSystem; species=0, region=0,factor=0,value=0)=boundary_robin!(system,species,region,factor,value)


"""
     boundary_robin!(y,u,bnode,ispec,ireg,fac,val)

Set Robin boundary condition for species ispec at boundary ibc.
"""
boundary_robin!(y,u,bnode,ispec,ireg,fac,val) =      bnode.region == ireg ? y[ispec] += fac*u[ispec]-val : nothing


"""
     boundary_robin!(y,u,bnode, args...; kwargs...)
Keyword argument version:
- `species`: species number. Default: 1
- `region`: boundary region number. By default, all boundary regions.
- `factor`: factor
- `value`: value
"""
boundary_robin!(y,u,bnode,args...; species=1, region=bnode.region,factor=0,value=0) = boundary_robin!(y,u,bnode,species,region,factor,value)

##################################################################
"""
$(SIGNATURES)

Number of species (size of first dimension) of solution array.
"""
num_species(a::AbstractArray)=size(a,1)







#
# Initialize Dirichlet BC
#
function _initialize_dirichlet!(U::AbstractMatrix,system::AbstractSystem{Tv, Tc, Ti, Tm}; time=0.0, λ=0.0, params::Vector{Tp}=Float64[]) where{Tv,Tp,Tc,Ti,Tm}
    _bfaceregions=bfaceregions(system.grid)
    _bfacenodes=bfacenodes(system.grid)
    nspecies=num_species(system)

    # set up bnode
    bnode = BNode(system,time,λ,params)
    bnodeparams=(bnode,)
    data=system.physics.data
    if isdata(data)
        bnodeparams=(bnode,data,)
    end

    # setup unknowns to be passed
    UK=zeros(Tv,num_species(system)+length(params))
    for iparm=1:length(params)
        UK[nspecies+iparm]=params[iparm]
    end
    u=unknowns(bnode,UK)
    
    # right hand side to be passed
    y=rhs(bnode,zeros(Tv,num_species(system)))

    # loop over all boundary faces
    for ibface=1:num_bfaces(system.grid)
        ibreg=_bfaceregions[ibface]

        # loop over all nodes of boundary face        
        for inode=1:dim_grid(system.grid)
            # Set Diichlet values to uninitialized
            _fill!(bnode,inode,ibface)
            bnode.dirichlet_value.=Inf
            jnode=_bfacenodes[inode,ibface]
            # set up solution vector, call boundary reaction
            @views UK[1:nspecies].=U[:,jnode]
            system.physics.breaction(y,u,bnodeparams...)

            # Check for Dirichlet bc
            for ispec=1:nspecies
                # Dirichlet bc given in breaction
                if !isinf(bnode.dirichlet_value[ispec])
                    U[ispec,jnode]=bnode.dirichlet_value[ispec]
                end
                
                # Dirichlet bc given after system creation (old API)
                if system.boundary_factors[ispec,ibreg]≈ Dirichlet
                    U[ispec,jnode]=system.boundary_values[ispec,ibreg]
                end
            end
        end
    end
end



function _initialize!(U::AbstractMatrix,system::AbstractSystem; time=0.0, λ=0.0, params=Number[] )
    _initialize_dirichlet!(U,system; time, λ,params)
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

Detect if array fits to the system.
"""
isunknownsof(u::Any, sys::AbstractSystem)=false
isunknownsof(u::DenseSolutionArray, sys::DenseSystem) = size(u) == size(sys.node_dof)
isunknownsof(u::SparseSolutionArray, sys::SparseSystem) = size(u) == size(sys.node_dof)


"""
$(SIGNATURES)

Create a solution vector for system.
If inival is not specified, the entries of the returned vector are undefined.
"""
unknowns(sys::AbstractSystem{Tv};inival=undef) where {Tv}=unknowns(Tv,sys,inival=inival)


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

LinearAlgebra.norm(system::SparseSystem,u,p)=LinearAlgebra.norm(u.node_dof.nzval,p)

######################################
# History
"""
    history(sys)

Return solver history from last `solve` call, if `log` was set to true.
See  see [`NewtonSolverHistory`](@ref), [`TransientSolverHistory`](@ref).
"""
history(sys::AbstractSystem)=sys.history


"""
    history_details(sys)

Return details of solver history from last `solve` call, if `log` was set to true.
See [`details`](@ref).
"""
history_details(sys::AbstractSystem)=details(sys.history)


"""
    history_summary(sys)

Return summary of solver history from last `solve` call, if `log` was set to true.
"""
history_summary(sys::AbstractSystem)=summary(sys.history)

####################################################################
# LEGACY

"""
````
System(grid,physics; kwargs...)
````
Create system with physics record.
!!! info  
    Starting with version 0.14, all physics data can be passed directly to the system constructor
"""
function System(grid::ExtendableGrid,physics::Physics;
                valuetype=coord_type(grid),
                indextype=index_type(grid),
                unknown_storage=:dense,
                matrixindextype=Int64,
                check_allocs=default_check_allocs(), kwargs...)

    system=System(grid; valuetype,indextype, unknown_storage, matrixindextype, check_allocs, kwargs...)
    physics!(system,physics)
end


"""
$(SIGNATURES)

Constructor for DenseSystem.
!!! compat  
    Will be removed in future versions
"""
DenseSystem(grid,physics::Physics; matrixindextype=Int64)=System(grid,physics,matrixindextype=matrixindextype,unknown_storage=:dense)


"""
$(SIGNATURES)

Constructor for SparseSystem.
!!! compat  
    Will be removed in future versions
"""
SparseSystem(grid,physics::Physics; matrixindextype=Int64)=System(grid,physics,matrixindextype=matrixindextype,unknown_storage=:sparse)




