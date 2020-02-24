##########################################################
"""
$(TYPEDEF)
    
Abstract type for physics.
"""
abstract type AbstractPhysics end



##########################################################
"""
$(TYPEDEF)
    
Abstract type for user data.

"""
abstract type AbstractData end

#
# Dummy data 
#
mutable struct NoData <: AbstractData
    NoData()=new() 
end
const nodata=NoData()

#
# Dummy callbacks
#
function nofunc(f,u,node,data)
end

function default_storage(f,u,node,data)
    f.=u
end

function nofunc2(f,u,node)
end

function default_storage2(f,u,node)
    f.=u
end


##########################################################
"""
$(TYPEDEF)
    
Physics data record.

$(TYPEDFIELDS) 
"""
struct Physics{Flux<:Function,
               Reaction<:Function,
               Storage<:Function,
               Source<:Function,
               BReaction<:Function,
               BStorage<:Function,
               D<:AbstractData} <: AbstractPhysics
    """
    Flux between neigboring control volumes
    """
    flux::Flux

    """
    Storage term (term under time derivative)
    """
    storage::Storage

    """
    Reaction term
    """
    reaction::Reaction

    """"
    Source term
    """
    source::Source

    """
    Boundary reaction term
    """
    breaction::BReaction

    """
    Boundary storage term
    """
    bstorage::BStorage

    """
    User data (parameters)
    """
    data::D

    """ 
    Number of species
    """
    num_species::Int8
end

##########################################################
"""
$(TYPEDSIGNATURES)

Constructor for physics data with default values.
"""
function Physics(;num_species=1,
                 data=nodata,
                 flux::Function=nofunc,
                 reaction::Function=nofunc,
                 storage::Function=default_storage,
                 source::Function=nofunc,
                 breaction::Function=nofunc,
                 bstorage::Function=nofunc
                 )
    
    return Physics(flux,
                   storage,
                   reaction,
                   source,
                   breaction,
                   bstorage,
                   data,
                   Int8(num_species)
                   )
end


##########################################################
"""
$(TYPEDSIGNATURES)

Constructor for physics data with default values.
"""
function FVMPhysics(;num_species=1,
                    data=nodata,
                    flux::Function=nofunc2,
                    reaction::Function=nofunc2,
                    storage::Function=default_storage2,
                    source::Function=nofunc2,
                    breaction::Function=nofunc2,
                    bstorage::Function=nofunc2
                    )
    
    return Physics(flux,
                   storage,
                   reaction,
                   source,
                   breaction,
                   bstorage,
                   data,
                   Int8(num_species)
                   )
end





function Base.show(io::IO,physics::AbstractPhysics)
    str=@sprintf("VoronoiFVM.Physics(num_species=%d",physics.num_species)
    if physics.data!=nodata
        str=str*", data=$(typeof(physics.data))"
    end
    function addfunc(func,name)
        if func!=nofunc
            str=str*", $(name)=$(nameof(func))"
        end
    end

    for name in fieldnames(typeof(physics))
        if (name!=:num_species)  && (name!=:data) && getfield(physics,name)!=nofunc
             str=str*", $(name)=$(nameof(getfield(physics,name)))"
        end
    end
    str=str*")"
    println(io,str)
end
