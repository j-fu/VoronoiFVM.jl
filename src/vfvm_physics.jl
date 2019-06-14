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
# Dummy callback
#
function nofunc(f,u,node, data)
end


##########################################################
"""
$(TYPEDEF)
    
Physics data record.

$(FIELDS) 
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
                 storage::Function=nofunc,
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


