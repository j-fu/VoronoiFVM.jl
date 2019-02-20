
abstract type NodeBase end

"""
    mutable struct Node

Structure holding local node information.
Fields:

    index::Int32
    coord::Array{Float64,1}

"""
mutable struct Node <: NodeBase
    index::Int32
    coord::Array{Float64,1}
    Node()=new()
end

abstract type EdgeBase end

"""
    mutable struct Edge

Structure holding local edge information.

Fields:

    index::Int32
    nodeK::Int32
    nodeL::Int32
    coordK::Array{Float64,1}
    coordL::Array{Float64,1}


"""
mutable struct Edge <: EdgeBase
    index::Int32
    nodeK::Int32
    nodeL::Int32
    coordK::Array{Float64,1}
    coordL::Array{Float64,1}
    Edge()=new()
end

"""
   function edgelength(edge::Edge)
   
Calculate the length of an edge. 
"""
function edgelength(edge::Edge)
    l::Float64
    l=0.0
    for i=1:length(edge.coordK)
        d=edge.coordK[i]-edge.coordL[i]
        l=l+d*d
    end
    return l
end
