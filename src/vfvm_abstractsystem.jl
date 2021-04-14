##########################################################
"""
$(TYPEDEF)
    
Abstract type for finite volume system structure.
"""
abstract type AbstractSystem{Tv<:Number, Ti <:Integer, Tm <:Integer} end



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
function System(grid,physics::Physics; unknown_storage=:dense, matrixindextype=Int32, check_allocs=false)
    if Symbol(unknown_storage)==:dense
        sys=DenseSystem(grid,physics, matrixindextype=matrixindextype)
    elseif Symbol(unknown_storage)==:sparse
        sys=SparseSystem(grid,physics, matrixindextype=matrixindextype)
    else
        throw("specify either unknown_storage=:dense  or unknown_storage=:sparse")
    end
    check_allocs!(sys,check_allocs)
end


"""
```
check_allocs!(system,true_or_false)
```

Enable/disable checking for time-consuming allocations in the assembly loop. Unless the matrix
pattern changes, there shouldn't occur any allocations in this loop. The check
method is aware of matrix pattern changes. As a consequence, allocations  in the assembly loop are
mostly due to type instabilities in physics callbacks.

By default, this check is switched on. 

Type instabilities can be introduced by variables in the closure of some physics callback.
They can be debugged via the `@time` macro applied to the expressions in a physics callback.

The following cases provide some ideas where to look for reasons of the problem and possible remedies:

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


Checking can be switched off via `check_allocs!(system,false)`.
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
        throw(DomainError(ispec,"Number of species exceeded"))
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
        throw(DomainError(ispec,"Number of species exceeded"))
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
function update_grid!(system::AbstractSystem{Tv,Ti}; grid=system.grid) where {Tv,Ti}
    
    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]
    bgeom=grid[BFaceGeometries][1]
    bfacenodes=grid[BFaceNodes]
    nbfaces=num_bfaces(grid)
    ncells=num_cells(grid)

    system.cellnodefactors=zeros(Tv,num_nodes(geom),ncells)
    system.celledgefactors=zeros(Tv,num_edges(geom),ncells)
    system.bfacenodefactors=zeros(Tv,num_nodes(bgeom),nbfaces)

    for icell=1:ncells
        @views cellfactors!(geom,csys,coord,cellnodes,icell,system.cellnodefactors[:,icell],system.celledgefactors[:,icell])
    end
    
    for ibface=1:nbfaces
        @views bfacefactors!(bgeom,csys,coord,bfacenodes,ibface,system.bfacenodefactors[:,ibface])
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
num_species(system::AbstractSystem) = system.physics.num_species



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




function Base.show(io::IO,sys::AbstractSystem)
    str=@sprintf("%s(num_species=%d)",typeof(sys),sys.physics.num_species)
    println(io,str)
end

#####################################################
has_generic_operator(sys::AbstractSystem) = sys.physics.generic_operator!=nofunc_generic
has_generic_operator_sparsity(sys::AbstractSystem) =  sys.physics.generic_operator_sparsity!=nofunc_generic_sparsity
