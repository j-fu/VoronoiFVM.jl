##########################################################
"""
$(TYPEDEF)
    
Abstract type for finite volume system structure
"""
abstract type AbstractSystem{Tv<:Number, Ti <:Integer} end



##################################################################
"""
$(TYPEDSIGNATURES)

Create Finite Volume System. 

- `grid`: 1D/2D/3D discretization grid
- `physics`: Physics struct containing node and edge callbacks
- `unknown_storage`: string or symbol:
     - `:dense` :  solution vector is an  `nspecies` x `nnodes`  dense matrix
     - `:sparse` :  solution vector is an `nspecies` x `nnodes`  sparse matrix

"""
function System(grid,physics::Physics; unknown_storage=:sparse)
    if Symbol(unknown_storage)==:dense
        return DenseSystem(grid,physics)
    elseif Symbol(unknown_storage)==:sparse
        return SparseSystem(grid,physics)
    else
        throw("specify either unknown_storage=:dense  or unknown_storage=:sparse")
    end
end

##################################################################
"""
Constant to be used as boundary condition factor 
to mark Dirichlet boundary conditons.    
"""
const Dirichlet=1.0e30


##################################################################
"""
$(TYPEDSIGNATURES)

Check if species number corresponds to boundary species.
"""
function is_boundary_species(this::AbstractSystem, ispec::Integer)
    isbspec=false
    for ibreg=1:num_bfaceregions(this.grid)
        if this.bregion_species[ispec,ibreg]>0
            isbspec=true
        end
    end
    return isbspec
end

##################################################################
"""
$(TYPEDSIGNATURES)

Check if species number corresponds bulk species.
"""
function is_bulk_species(this::AbstractSystem, ispec::Integer)
    isrspec=false
    for ixreg=1:num_cellregions(this.grid)
        if this.region_species[ispec,ixreg]>0
            isrspec=true
        end
    end
    return isrspec
end

##################################################################
"""
$(TYPEDSIGNATURES)

Add species to a list of bulk regions. Species numbers for
bulk and boundary species have to be distinct.
"""
function enable_species!(this::AbstractSystem,ispec::Integer, regions::AbstractVector)
    if ispec>num_species(this)
        throw(DomainError(ispec,"Number of species exceeded"))
    end
    if is_boundary_species(this,ispec)
        throw(DomainError(ispec,"Species is already boundary species"))
    end
    _cellregions=cellregions(this.grid)
    _cellnodes=cellnodes(this.grid)
    for i in eachindex(regions)
        ireg=regions[i]
        this.region_species[ispec,ireg]=ispec
        for icell=1:num_cells(this.grid)
            if _cellregions[icell]==ireg
                for iloc=1:size(_cellnodes,1)
                    iglob=_cellnodes[iloc,icell]
                    this.node_dof[ispec,iglob]=ispec
                end
            end
        end
    end
end


##################################################################
"""
$(TYPEDSIGNATURES)

Add species to a list of boundary regions. Species numbers for
bulk and boundary species have to be distinct.

"""
function enable_boundary_species!(this::AbstractSystem, ispec::Integer, bregions::AbstractVector)
    if ispec>num_species(this)
        throw(DomainError(ispec,"Number of species exceeded"))
    end
    if is_bulk_species(this,ispec)
        throw(DomainError(ispec,"Species is already bulk species"))
    end
    _bfaceregions=bfaceregions(this.grid)
    _bfacenodes=bfacenodes(this.grid)
    
    for i in eachindex(bregions)
        ireg=bregions[i]
        this.bregion_species[ispec,ireg]=ispec
        for ibface=1:num_bfaces(this.grid)
            if _bfaceregions[ibface]==ireg
                for iloc=1:size(_bfacenodes,1)
                    iglob=_bfacenodes[iloc,ibface]
                    this.node_dof[ispec,iglob]=ispec
                end
            end
        end
    end
end


# Create matrix in system and figure out if species
# distribution is homgeneous
function _complete!(this::AbstractSystem{Tv,Ti};create_newtonvectors=false) where {Tv,Ti}

    if isdefined(this,:matrix)
        return
    end
    this.matrix=ExtendableSparseMatrix{Tv,Ti}(num_dof(this), num_dof(this))
    this.species_homogeneous=true
    species_added=false
    for inode=1:size(this.node_dof,2)
        for ispec=1:size(this.node_dof,1)
            if this.node_dof[ispec,inode]==ispec
                species_added=true
            else
                this.species_homogeneous=false
            end
        end
    end
    
    if (!species_added)
        error("No species enabled.\n Call enable_species(system,species_number, list_of_regions) at least once.")
    end
    if create_newtonvectors
        this.residual=unknowns(this)
        this.update=unknowns(this)
    end
    update_grid!(this)
end


"""
$(SIGNATURES)

Update grid (e.g. after rescaling of coordinates).
"""
function update_grid!(this::AbstractSystem{Tv,Ti};grid=this.grid) where{Tv, Ti}

    geom=grid[CellGeometries][1]
    csys=grid[CoordinateSystem]
    coord=grid[Coordinates]
    cellnodes=grid[CellNodes]
    cellregions=grid[CellRegions]
    bgeom=grid[BFaceGeometries][1]
    bfacenodes=grid[BFaceNodes]
    nbfaces=num_bfaces(grid)
    ncells=num_cells(grid)

    this.cellnodefactors=zeros(Tv,num_nodes(geom),ncells)
    this.celledgefactors=zeros(Tv,num_edges(geom),ncells)
    this.bfacenodefactors=zeros(Tv,num_nodes(bgeom),nbfaces)

    for icell=1:ncells
        @views cellfactors!(geom,csys,coord,cellnodes,icell,this.cellnodefactors[:,icell],this.celledgefactors[:,icell])
    end
    
    for ibface=1:nbfaces
        @views bfacefactors!(bgeom,csys,coord,bfacenodes,ibface,this.bfacenodefactors[:,ibface])
    end
end

##################################################################
"""
$(TYPEDSIGNATURES)
    
Check if degree of freedom is defined.
"""
isdof(this::AbstractSystem,ispec,inode)= this.node_dof[ispec,inode]==ispec ? true : false


##################################################################
"""
$(TYPEDSIGNATURES)

Number of species in system
"""
num_species(this::AbstractSystem) = this.physics.num_species



##################################################################
"""
$(TYPEDSIGNATURES)

Retrieve user data record.
"""
data(this::AbstractSystem) = this.physics.data


##################################################################
"""
$(TYPEDSIGNATURES)

Set Dirichlet boundary conditon for species ispec at boundary ibc.
"""
function boundary_dirichlet!(this::AbstractSystem, ispec::Integer, ibc::Integer, val)
    this.boundary_factors[ispec,ibc]=Dirichlet
    this.boundary_values[ispec,ibc]=val
end


##################################################################
"""
$(TYPEDSIGNATURES)

Set Neumann boundary conditon for species ispec at boundary ibc.
"""
function boundary_neumann!(this::AbstractSystem, ispec::Integer, ibc::Integer, val)
    this.boundary_factors[ispec,ibc]=0.0
    this.boundary_values[ispec,ibc]=val
end


##################################################################
"""
$(TYPEDSIGNATURES)

Set Robin boundary conditon for species ispec at boundary ibc.
"""
function boundary_robin!(this::AbstractSystem, ispec::Integer, ibc::Integer,fac, val)
    this.boundary_factors[ispec,ibc]=fac
    this.boundary_values[ispec,ibc]=val
end



##################################################################
"""
$(TYPEDSIGNATURES)

Number of species (size of first dimension) of solution array.
"""
num_species(a::AbstractArray)=size(a,1)





#
# Initialize Dirichlet BC
#
function _initialize_dirichlet!(U::AbstractMatrix,this::AbstractSystem)
    _bfaceregions=bfaceregions(this.grid)
    _bfacenodes=bfacenodes(this.grid)
    for ibface=1:num_bfaces(this.grid)
        ibreg=_bfaceregions[ibface]
        for ispec=1:num_species(this)
            if this.boundary_factors[ispec,ibreg]≈ Dirichlet
                for inode=1:dim_grid(this.grid)
                    U[ispec,_bfacenodes[inode,ibface]]=this.boundary_values[ispec,ibreg]
                end
            end
        end
    end
end




function _initialize!(U::AbstractMatrix,this::AbstractSystem)
    _initialize_dirichlet!(U,this)
    _initialize_inactive_dof!(U,this)
end




function Base.show(io::IO,sys::AbstractSystem)
    str=@sprintf("%s(num_species=%d)",typeof(sys),sys.physics.num_species)
    println(io,str)
end

