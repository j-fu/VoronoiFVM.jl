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

!!! Deprecate
There is no benefit in inheriting from `AbstractData` and this type
will be removed in a future release.
"""
abstract type AbstractData end

struct NoData <: AbstractData
end

struct DummyData <: AbstractData
end

isdata(::NoData)=false
isdata(::Any)=true

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

function nofunc_generic(f,u,sys)
end

function nofunc_generic_sparsity(sys)
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
               GenericOperator<:Function,
               GenericOperatorSparsity<:Function,
               Data} <: AbstractPhysics
    """
    Flux between neigboring control volumes: `flux(f,_u,edge)` or `flux(f,_u,edge,data)`
    should return in `f[i]` the flux of species i along the edge joining circumcenters
    of neigboring control volumes. `u=unknowns(_u)` returns a 2D array such that for species i,
    `u[i,1]` and `u[i,2]` contain the unknown values at the corresponding ends of the edge.
    """
    flux::Flux

    """
    Storage term (term under time derivative): `storage(f,u,node)` or `storage(f,u,node,data)` 

    It should return in `f[i]` the storage term for the i-th equation. `u[i]` contains the value of
    the i-th unknown.
    """
    storage::Storage

    """
    Reaction term:  `reaction(f,u,node)` or `reaction(f,u,node,data)` 

    It should return in `f[i]` the reaction term for the i-th equation. `u[i]` contains the value of
    the i-th unknown.
    """
    reaction::Reaction

    """
    Source term: `source(f,node)` or `source(f,node,data)`.

    It should return the in `f[i]` the value of the source term for the i-th equation.
    """
    source::Source

    """
    Boundary reaction term:  `breaction(f,u,node)` or `breaction(f,u,node,data)` 
    Similar to reaction, but restricted to the inner or outer boundaries.
    """
    breaction::BReaction

    """
    Boundary storage term: `bstorage(f,u,node)` or `bstorage(f,u,node,data)` 
    Similar to storage, but restricted to the inner or outer boundaries.
    """
    bstorage::BStorage

    """
    Generic operator  `generic_operator(f,u,sys)`. 
    This operator acts on the full solution `u` of a system. Sparsity
    is detected automatically  unless `generic_operator_sparsity` is given.
    """
    generic_operator::GenericOperator

    
    """
    Function defining the sparsity structure of the generic operator.
    This should return the sparsity pattern of the `generic_operator`.
    """
    generic_operator_sparsity::GenericOperatorSparsity

    
    """
    User data (parameters).
    This allows to pass various parameters to the callback functions.
    """
    data::Data

    """
    Number of species including boundary species.
    """
    num_species::Int8

end

##########################################################
"""
$(SIGNATURES)

Constructor for physics data with default values for the
constitutive callback functions.

There are two variants of this constructor. It `data` is given, all callback functions
should accept a last `data` argument. Otherwise, no data are passed explicitely, and it
is assumed that constitutive callbacks take parameters from the closure where the function
is defined.
"""
function Physics(;num_species=1,
                 data=NoData(),
                 flux::Function=nofunc,
                 reaction::Function=nofunc,
                 storage::Function=default_storage,
                 source::Function=nofunc,
                 breaction::Function=nofunc,
                 bstorage::Function=nofunc,
                 generic::Function=nofunc_generic,
                 generic_sparsity::Function=nofunc_generic_sparsity
                 )
    if !isdata(data)
        flux==nofunc ? flux=nofunc2 : true
        reaction==nofunc ? reaction=nofunc2 : true
        storage==default_storage ? storage=default_storage2 : true
        source==nofunc ? source=nofunc2 : true
        breaction==nofunc ? breaction=nofunc2 : true
        bstorage==nofunc ? bstorage=nofunc2 : true
    end
    
    return Physics(flux,
                   storage,
                   reaction,
                   source,
                   breaction,
                   bstorage,
                   generic,
                   generic_sparsity,
                   data,
                   Int8(num_species)
                   )
end

"""
$(SIGNATURES)

Check if physics object has data
"""
hasdata(physics::Physics)=isdata(physics.data)


"""
$(SIGNATURES)

Show physics object
"""
function Base.show(io::IO,physics::AbstractPhysics)
    str=@sprintf("VoronoiFVM.Physics(num_species=%d",physics.num_species)
    if isdata(physics.data)
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
