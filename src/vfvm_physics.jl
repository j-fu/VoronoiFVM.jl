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
!!! compat  
    Will be removed in future versions
"""
abstract type AbstractData end

isdata(::Nothing)=false
isdata(::Any)=true

#
# Dummy callbacks
#
function nofunc(f,u,node,data)
end

function nosrc(f,node,data)
end

function default_storage(f,u,node,data)
    f.=u
end

function nofunc2(f,u,node)
end

function nosrc2(f,node)
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
````
struct Physics
````

Physics data record with the following fields:

$(TYPEDFIELDS)
"""
struct Physics{Flux<:Function,
               Reaction<:Function,
               Storage<:Function,
               Source<:Function,
               BFlux<:Function,
               BReaction<:Function,
               BSource<:Function,
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
    Flux between neighboring control volumes on the boundary
    """
    bflux::BFlux

    """
    Boundary reaction term:  `breaction(f,u,node)` or `breaction(f,u,node,data)` 
    Similar to reaction, but restricted to the inner or outer boundaries.
    """
    breaction::BReaction


    """
    Boundary source term: `bsource(f,node)` or `bsource(f,node,data)`.

    It should return in `f[i]` the value of the source term for the i-th equation.
    """
    bsource::BSource
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
````
Physics(;num_species=1,
         data=nothing,
         flux,
         reaction,
         storage,
         source,
         breaction,
         bstorage,
         generic,
         generic_sparsity
)
````

Constructor for physics data. For the meaning of the optional function
callbacks, see [`Physics`](@ref).

There are two variants of this constructor. It `data` is given, all callback functions
should accept a last `data` argument. Otherwise, no data are passed explicitely, and it
is assumed that constitutive callbacks take parameters from the closure where the function
is defined.


"""
function Physics(;num_species=1,
                 data=nothing,
                 flux::Function=nofunc,
                 reaction::Function=nofunc,
                 storage::Function=default_storage,
                 source::Function=nosrc,
                 bflux::Function=nfunc,
                 breaction::Function=nofunc,
                 bsource::Function=nosrc,
                 bstorage::Function=nofunc,
                 generic::Function=nofunc_generic,
                 generic_sparsity::Function=nofunc_generic_sparsity
                 )
    if !isdata(data)
        flux==nofunc ? flux=nofunc2 : true
        reaction==nofunc ? reaction=nofunc2 : true
        storage==default_storage ? storage=default_storage2 : true
        source==nosrc ? source=nosrc2 : true
        bflux==nofunc ? bflux=nfunc2 : true
        breaction==nofunc ? breaction=nofunc2 : true
        bsource==nosrc ? bsource=nosrc2 : true
        bstorage==nofunc ? bstorage=nofunc2 : true
    end
    
    return Physics(flux,
                   storage,
                   reaction,
                   source,
                   bflux,
                   breaction,
                   bsource,
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




"""
```
@create_physics_wrappers(physics,node,bnode,edge)
```

Create wrapper functions around physics callbacks which
fit the API of `ForwardDiff.jacobian!` and pass the `data` parameter
if necessary. These are meant to  be defined
before performing assembly loops. The macro creates the follwing variables:

- wrapper functions: `fluxwrap`,`storagewrap`,`reactionwrap`,`bstoragewrap`, `breactionwrap`
- flag variables: `issource`, `isreaction`,`isbreaction`,`isbstorage`

"""
macro create_physics_wrappers(physics,node,bnode,edge)
    return quote
        data=$(esc(physics)).data
        if isdata(data)
            global issource=($(esc(physics)).source!=nofunc)
            global isreaction=($(esc(physics)).reaction!=nofunc)
            global isbreaction=($(esc(physics)).breaction!=nofunc)
            global isbsource=($(esc(physics)).source!=nofunc)
            global isbstorage=($(esc(physics)).bstorage!=nofunc)
            
            global fluxwrap=function(y, u)
                y.=0
                $(esc(physics)).flux(y,u,$(esc(edge)),data)
                nothing
            end
            
            global reactionwrap=function (y, u)
                y.=0
                ## for ii in ..  uu[node.speclist[ii]]=u[ii]
                $(esc(physics)).reaction(y,u,$(esc(node)),data)
                ## for ii in .. y[ii]=y[node.speclist[ii]]
                nothing
            end
            
            global storagewrap= function(y, u)
                y.=0
                $(esc(physics)).storage(y,u,$(esc(node)),data)
                nothing
            end
        
            global sourcewrap=function(y)
                y.=0
                $(esc(physics)).source(y,$(esc(node)),data)
                nothing
            end

            global bfluxwrap=function(y, u)
                y.=0
                $(esc(physics)).bflux(y,u,$(esc(edge)),data)
                nothing
            end
            
            global breactionwrap=function(y, u)
                y.=0
                $(esc(physics)).breaction(y,u,$(esc(bnode)),data)
                nothing
            end

            global bsourcewrap=function(y)
                y.=0
                $(esc(physics)).bsource(y,$(esc(bnode)),data)
                nothing
            end
            
            global bstoragewrap=function(y, u)
                y.=0
                $(esc(physics)).bstorage(y,u,$(esc(bnode)),data)
                nothing
            end
        
        else
            global issource    = ($(esc(physics)).source    != nofunc2)
            global isreaction  = ($(esc(physics)).reaction  != nofunc2)
            global isbreaction = ($(esc(physics)).breaction != nofunc2)
            global isbsource   = ($(esc(physics)).source    != nofunc2)
            global isbflux     = ($(esc(physics)).bflux     != nofunc2)
            global isbstorage  = ($(esc(physics)).bstorage  != nofunc2)
            
            global fluxwrap=function(y, u)
                y.=0
                $(esc(physics)).flux(y,u,$(esc(edge)))
                nothing
            end
            
            global reactionwrap=function(y, u)
                y.=0
                ## for ii in ..  uu[node.speclist[ii]]=u[ii]
                $(esc(physics)).reaction(y,u,$(esc(node)))
                ## for ii in .. y[ii]=y[node.speclist[ii]]
                nothing
            end
            
            global storagewrap= function(y, u)
                y.=0
                $(esc(physics)).storage(y,u,$(esc(node)))
                nothing
            end
            
            global sourcewrap=function(y)
                y.=0
                $(esc(physics)).source(y,$(esc(node)))
                nothing
            end

            global bfluxwrap=function(y, u)
                y.=0
                $(esc(physics)).bflux(y,u,$(esc(edge)))
                nothing
            end
            
            global breactionwrap=function(y, u)
                y.=0
                $(esc(physics)).breaction(y,u,$(esc(bnode)))
                nothing
            end

            global bsourcewrap=function(y)
                y.=0
                $(esc(physics)).bsource(y,$(esc(bnode)))
                nothing
            end
            
            global bstoragewrap=function(y, u)
                y.=0
                $(esc(physics)).bstorage(y,u,$(esc(bnode)))
                nothing
            end
            
        end
    end
end
