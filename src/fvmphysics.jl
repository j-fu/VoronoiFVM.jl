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
defined via the macro [`@AddFVMPhysicsBaseClassFields`](@ref).
These are

    num_species::Int64 # number of species in the interior
    num_bspecies::Array{Int64,1} # number of species in boundary regions
    bregion::Int64 # number of current boundary region
    flux::Function # flux function
    reaction::Function # reaction term in interior
    breaction::Function # reaction term at  boundary
    source::Function # source term
    storage::Function # storage term
    bstorage::Function # storage term at boundary

"""
abstract type FVMPhysics end


###########################################################
"""
Define "Base class" fields to be pasted into all structs
of type FVMPhysics.
"""
@define_field_group AddFVMPhysicsBaseClassFields begin
    num_species::Int64 # number of species in the interior
    num_bspecies::Array{Int64,1} # number of species in boundary regions
    bregion::Int64 # number of current boundary region
    flux::Function # flux function
    reaction::Function # reaction term in interior
    breaction::Function # reaction term at  boundary
    source::Function # source term
    storage::Function # storage term
    bstorage::Function # storage term at boundary
end




###########################################################
"""
    Prototype flux term for [`FVMPhysics`](@ref).

````
function prototype_flux!(
    this::FVMPhysics,  # physics data
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
    this::FVMPhysics,  # physics data
    f::AbstractArray,  # result vector
    uk::AbstractArray, # unknowns on "left" end of edge
    ul::AbstractArray) # unknowns on "right" end of edge

    for i=1:this.num_species
        f[i]=uk[i]-ul[i]
    end
end



###########################################################
"""
Prototype storage term for [`FVMPhysics`](@ref)

````
function prototype_storage!(
    this::FVMPhysics,  # physics data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )
````

Evaluate storage term vector \$f(u)=(f_1(u)\\dots f_n(u))\$
for given vector of nodal unknowns  \$u=(u_1\\dots u_n)\$.


"""
function prototype_storage!(
    this::FVMPhysics,  # physics data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )

    for i=1:this.num_species
        f[i]=u[i]
    end
end

###########################################################
"""
Prototype source term for [`FVMPhysics`](@ref).

````
function prototype_source!(
    this::FVMPhysics,   # physics data
    f::AbstractArray, # result vector 
    x::AbstractArray  # coordinate vector
    )
````

Evaluate source term vector \$f(x)=(f_1(x)\\dots f_n(x))\$ in point \$x=(x_1\\dots x_d)\$.


"""
function prototype_source!(
    this::FVMPhysics,   # physics data
    f::AbstractArray, # result vector 
    x::AbstractArray  # coordinate vector
    )

    for i=1:this.num_species
        f[i]=0
    end
end

###########################################################
"""
Prototype reaction term for [`FVMPhysics`](@ref).

````
function prototype_reaction!(
    this::FVMPhysics,  # physics data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )
````

Evaluate reaction term vector \$f(u)=(f_1(u)\\dots f_n(u))\$
for given vector of nodal unknowns  \$u=(u_1\\dots u_n)\$.

"""
function prototype_reaction!(
    this::FVMPhysics,  # physics data
    f::AbstractArray,# result vector
    u::AbstractArray # vector of unknowns in point
    )

    for i=1:this.num_species
        f[i]=0
    end
end


###########################################################
"""
    Prototype boundary reaction term for [`FVMPhysics`](@ref).

````
function prototype_breaction!(
    this::FVMPhysics,    # physics data
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

When called, `this.bregion` contains the boundary region number,
and  `this.num_bspecies[this.bregion]` contains the number \$n_b\$ of
boundary species present in the boundary node.

"""
function prototype_breaction!(
    this::FVMPhysics,    # physics data
    f::AbstractArray,  # result vector for interior species
    bf::AbstractArray, # result vector for boundary species
    u::AbstractArray,  # interior unknowns
    bu::AbstractArray  # boundary unknowns
    )
    for i=1:this.num_species
        f[i]=0
    end
    for i=1:this.num_bspecies[this.bregion]
        bf[i]=0
    end
end

###########################################################
"""
    Prototype boundary storage term for [`FVMPhysics`](@ref).

````
function prototype_bstorage!(
    this::FVMPhysics, # physics data
    bf::AbstractArray, # result vector for boundary species
    bu::AbstractArray  # boundary unknowns 
    )
````

Evaluate storage rates at boundary: 
\$bf(u,bu)=(bf_1(u,bu)\\dots bf_{n_b}(u,bu))\$ where
\$u=(u_1\\dots u_n)\$ is the value (trace) 
of interior unknowns at the boundary and 
\$bu=(bu_1\\dots bu_{n_b})\$ is the value of the boundary unknowns.   

When called, `this.bregion` contains the boundary region number,
and  `this.num_bspecies[this.bregion]` contains the number \$n_b\$ of
boundary species present in the boundary node.


"""
function prototype_bstorage!(
    this::FVMPhysics, # physics data
    bf::AbstractArray, # result vector for boundary species
    bu::AbstractArray  # boundary unknowns 
    )

    for i=1:this.num_bspecies[this.bregion]
        bf[i]=bu[i]
    end
end




###########################################################
"""
Base class for physics data
"""
mutable struct FVMPhysicsBase <: FVMPhysics
    @AddFVMPhysicsBaseClassFields
    FVMPhysicsBase(nspec::Int) = FVMPhysicsBase(new(), nspec)
end

###########################################################
"""
Base class intialization for physics dats

````
function FVMPhysicsBase(this::FVMPhysics,nspec::Int)
````

"""
function FVMPhysicsBase(this::FVMPhysics,nspec::Int)
    this.num_species=nspec
    this.num_bspecies=Array{Int64,1}(undef,0)
    this.bregion=0
    this.flux=prototype_flux!
    this.reaction=prototype_reaction!
    this.breaction=prototype_breaction!
    this.bstorage=prototype_bstorage!
    this.source=prototype_source!
    this.storage=prototype_storage!
    return this
end



