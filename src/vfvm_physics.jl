##########################################################
"""
   abstract type AbstractPhysics
    
Abstract type for physics.
"""
abstract type AbstractPhysics end



##########################################################
"""
    abstract type Data end
    
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
    struct Physics{Flux<:Function,
                   Reaction<:Function,
                   Storage<:Function,
                   Source<:Function,
                   BReaction<:Function,
                   BStorage<:Function,
                   D<:Data} <: AbstractPhysics
    flux::Flux
    storage::Storage
    reaction::Reaction
    source::Source
    breaction::BReaction
    bstorage::BStorage
    data::D
    num_species::Int8
end
    
Physics data record 
"""
struct Physics{Flux<:Function,
               Reaction<:Function,
               Storage<:Function,
               Source<:Function,
               BReaction<:Function,
               BStorage<:Function,
               D<:AbstractData} <: AbstractPhysics
    flux::Flux
    storage::Storage
    reaction::Reaction
    source::Source
    breaction::BReaction
    bstorage::BStorage
    data::D
    num_species::Int8
end

##########################################################
"""
    function Physics(;num_species=1,
                      data=nodata,
                      flux::Function=nofunc,
                       reaction::Function=nofunc,
                       storage::Function=nofunc,
                       source::Function=nofunc,
                       breaction::Function=nofunc,
                       storage::Function=nofunc
                     )

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


