
##########################################################
"""
   abstract type AbstractGrid

Abstract type for grid like datastructures.
"""
abstract type AbstractGrid end

##########################################################
"""
   abstract type Physics
    
Abstract type for user data record
"""
abstract type AbstractPhysics end



##########################################################
"""
       abstract type AbstractSystem
    
Abstract type for finite volume system structure
"""
abstract type AbstractSystem{Tv<:Number} end


abstract type Data end
