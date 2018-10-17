
"""
    See 
    https://discourse.julialang.org/t/inheritance-in-julia/5416/3
    http://www.stochasticlifestyle.com/type-dispatch-design-post-object-oriented-programming-julia/
    
    We have no inheritance, but duck typing, traits and composition...

    This macro allows to allows to define a group of fields which
    are assumed to be those of the base type.
"""
macro define_field_group(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end




"""
   Abstract type for user problem data.
 
   Derived types must contain field number_of_species::Int64
"""
abstract type FVMParameters end


"""
Define "Base class" fields to be pasted into all structs
of type FVMParameters
"""

@define_field_group AddDefaultFVMParameters begin
    num_species::Int64
    num_bspecies::Array{Int64,1}
    bregion::Int64
    flux::Function
    reaction::Function
    breaction::Function
    source::Function
    storage::Function
end




"""
    Default source term
"""
function default_source!(this::FVMParameters, f,x)
    for i=1:this.num_species
        f[i]=0
    end
end

"""
    Default reaction term
"""
function default_reaction!(this::FVMParameters, f,u)
    for i=1:this.num_species
        f[i]=0
    end
end

"""
    Default flux term
"""
function default_flux!(this::FVMParameters, f,uk,ul)
    for i=1:this.num_species
        f[i]=uk[i]-ul[i]
    end
end



"""
    Default storage term
"""
function default_storage!(this::FVMParameters, f,u)
    for i=1:this.num_species
        f[i]=u[i]
    end
end

"""
               Problem data type for default parameter function
"""
mutable struct DefaultFVMParameters <: FVMParameters
    @AddDefaultFVMParameters
    DefaultFVMParameters(nspec::Int) = DefaultFVMParameters(new(), nspec)
end

function DefaultFVMParameters(this::FVMParameters,nspec::Int)
    this.num_species=nspec
    this.num_bspecies=Array{Int64,1}(undef,0)
    this.bregion=0
    this.flux=default_flux!
    this.reaction=default_reaction!
    this.breaction=default_reaction!
    this.source=default_source!
    this.storage=default_storage!
end

