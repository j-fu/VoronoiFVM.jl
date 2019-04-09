###########################################################
"""
Support macro to define a number of fields in a struct.

See 
https://discourse.julialang.org/t/inheritance-in-julia/5416/3
http://www.stochasticlifestyle.com/type-dispatch-design-post-object-oriented-programming-julia/
    
In Julia we have no inheritance, but duck typing, traits and composition, without
concepts so far. This macro allows to to define a group of fields which
are assumed to be those of a base type.
"""
macro define_field_group(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end




###########################################################
"""
Abstract type for user problem data.
 
Any derived type must contain the fields which are 
defined via the macro [`@AddPhysicsBaseClassFields`](@ref).
These are

    flux::Function # flux function
    reaction::Function # reaction term in interior
    breaction::Function # reaction term at  boundary
    source::Function # source term
    storage::Function # storage term
    bstorage::Function # storage term at boundary

"""
abstract type Physics end


###########################################################
"""
Define "Base class" fields to be pasted into all structs
of type Physics.
"""
@define_field_group AddPhysicsBaseClassFields begin
    flux::Function # flux function
    reaction::Function # reaction term in interior
    breaction::Function # reaction term at  boundary
    source::Function # source term
    storage::Function # storage term
    bstorage::Function # storage term at boundary
end




###########################################################
"""
    Prototype flux term for [`Physics`](@ref).

````
function prototype_flux!(
    physics::Physics,  # physics data
    edge::Edge,        # edge data
    f::AbstractArray,  # result vector
    uk::AbstractArray, # unknowns on "left" end of edge
    ul::AbstractArray) # unknowns on "right" end of edge
````

Evaluate edge flux function 
\$f(u_k,u_l)=(f_1(u_k,u_l)\\dots f_n(u_k,u_l))\$,
where \$u_k=(u_{k,1}\\dots u_{k,n})\$ is the value 
of unknowns at "left" end of an edge and
 \$u_l=(u_{l,1}\\dots u_{l,n})\$
is the value of unknowns at "right" end of an edge.

"""
function prototype_flux!(
    physics::Physics,  # physics data
    edge::Edge,        # edge data
    f::AbstractArray,  # result vector
    uk::AbstractArray, # unknowns on "left" end of edge
    ul::AbstractArray) # unknowns on "right" end of edge

    for i=1:edge.num_species
        f[i]=uk[i]-ul[i]
    end
end



###########################################################
"""
Prototype storage term for [`Physics`](@ref)

````
function prototype_storage!(
    physics::Physics,  # physics data
    node::Node,        # node data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )
````

Evaluate storage term vector \$f(u)=(f_1(u)\\dots f_n(u))\$
for given vector of nodal unknowns  \$u=(u_1\\dots u_n)\$.


"""
function prototype_storage!(
    physics::Physics,  # physics data
    node::Node,        # node data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )

    for i=1:node.nspecies
        f[i]=u[i]
    end
end

###########################################################
"""
Prototype source term for [`Physics`](@ref).

````
function prototype_source!(
    physics::Physics,   # physics data
    node::Node,        # node data
    f::AbstractArray, # result vector 
    )
````

Evaluate source term vector \$f(x)=(f_1(x)\\dots f_n(x))\$ in point \$x=(x_1\\dots x_d)\$.


"""
function prototype_source!(
    physics::Physics,   # physics data
    node::Node,        # node data
    f::AbstractArray, # result vector 
    )

    for i=1:node.nspecies
        f[i]=0
    end
end

###########################################################
"""
Prototype reaction term for [`Physics`](@ref).

````
function prototype_reaction!(
    physics::Physics,  # physics data
    node::Node,        # node data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )
````

Evaluate reaction term vector \$f(u)=(f_1(u)\\dots f_n(u))\$
for given vector of nodal unknowns  \$u=(u_1\\dots u_n)\$.

"""
function prototype_reaction!(
    physics::Physics,  # physics data
    node::Node,        # node data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )

    for i=1:node.nspecies
        f[i]=0
    end
end


###########################################################
"""
    Prototype boundary reaction term for [`Physics`](@ref).

````
function prototype_breaction!(
    physics::Physics,    # physics data
    node::Node,        # node data
    f::AbstractArray,  # result vector for interior species
    bf::AbstractArray, # result vector for boundary species
    u::AbstractArray,  # interior unknowns
    bu::AbstractArray  # boundary unknowns
    )
````

Evaluate reaction rates at boundary: \$f(u,bu)=(f_1(u,bu)\\dots f_n(u,bu))\$,
\$bf(u,bu)=(bf_1(u,bu)\\dots bf_{n_b}(u,bu))\$ where
\$u=(u_1\\dots u_n)\$ is the value (trace) 
of interior unknowns at the boundary and 
\$bu=(bu_1\\dots bu_{n_b})\$ is the value of the boundary unknowns.   

When called, `physics.bregion` contains the boundary region number,
and  `physics.num_bspecies[physics.bregion]` contains the number \$n_b\$ of
boundary species present in the boundary node.

"""
function prototype_breaction!(
    physics::Physics,    # physics data
    node::Node,        # node data
    f::AbstractArray,  # result vector for interior species
    bf::AbstractArray, # result vector for boundary species
    u::AbstractArray,  # interior unknowns
    bu::AbstractArray  # boundary unknowns
    )
    for i=1:node.nspecies
        f[i]=0
    end
    for i=1:physics.num_bspecies[physics.bregion]
        bf[i]=0
    end
end

###########################################################
"""
    Prototype boundary storage term for [`Physics`](@ref).

````
function prototype_bstorage!(
    physics::Physics, # physics data
    node::Node,        # node data
    bf::AbstractArray, # result vector for boundary species
    bu::AbstractArray  # boundary unknowns 
    )
````

Evaluate storage rates at boundary: 
\$bf(u,bu)=(bf_1(u,bu)\\dots bf_{n_b}(u,bu))\$ where
\$u=(u_1\\dots u_n)\$ is the value (trace) 
of interior unknowns at the boundary and 
\$bu=(bu_1\\dots bu_{n_b})\$ is the value of the boundary unknowns.   

When called, `physics.bregion` contains the boundary region number,
and  `physics.num_bspecies[physics.bregion]` contains the number \$n_b\$ of
boundary species present in the boundary node.


"""
function prototype_bstorage!(
    physics::Physics, # physics data
    node::Node,        # node data
    bf::AbstractArray, # result vector for boundary species
    bu::AbstractArray  # boundary unknowns 
    )

    for i=1:physics.num_bspecies[physics.bregion]
        bf[i]=bu[i]
    end
end




###########################################################
"""
Base class for physics data
"""
mutable struct PhysicsBase <: Physics
    @AddPhysicsBaseClassFields
    PhysicsBase(nspec::Int) = PhysicsBase(new(), nspec)
end

###########################################################
"""
Base class intialization for physics dats

````
function PhysicsBase(physics::Physics,nspec::Int)
````

"""
function PhysicsBase(this::Physics,nspec::Int)
    this.flux=prototype_flux!
    this.reaction=prototype_reaction!
    this.breaction=prototype_breaction!
    this.bstorage=prototype_bstorage!
    this.source=prototype_source!
    this.storage=prototype_storage!
    return this
end



