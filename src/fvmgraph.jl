# Structure for discrete geometry (for finite volume method)
# This can be easily generalized to higher dimensions 
# and unstructured grids
# 
# This is a graph embedded in R^d (d=1,2,3) with edge weights and node weights
# which can be generated from discretization grids in 1D/2D/3D
#
struct FVMGraph
    SpaceDim::Int32 # space dimension
    NumberOfNodes::Int64 # number of nodes
    NumberOfBoundaryRegions::Int64 # number of boundary regions
    Nodes::Array{Float64,2} # point coordinates.
    Edges::Array{Int32,2} # index array pointing into points list
    EdgeFactors::Array{Float64,1} # interface area / edge length
    NodeFactors::Array{Float64,1} # control volume sizes
    BoundaryNodes::Array{Int32,1} # indices of dirichlet points 
    BoundaryRegions::Array{Int32,1} # boundary region numbers
    BoundaryNodeFactors::Array{Float64,1} # boundary control volume sizes
    
    # Constructor for 1D discrete geometry
    function FVMGraph(X::Array{Float64,1})
    #  Primal Grid:
    #  o-----o-----o-----o-----o-----o-----o-----o-----o
    # Dual grid with control volumes
    #  |--|-----|-----|-----|-----|-----|-----|-----|--|
        n=length(X)
        pts=reshape(X,1,n)
        eds=zeros(Int32,2,n-1)
        efac=zeros(n-1)
        nfac=zeros(n)
        for i=1:n-1 # loop over edges
            eds[1,i]=i
            eds[2,i]=i+1
            h=pts[i+1]-pts[i]
            efac[i]=1.0/h
            nfac[i]+=0.5*h
            nfac[i+1]+=0.5*h
        end
        bpoints=zeros(Int32,2)
        bregions=zeros(Int32,2)
        bnodefac=zeros(Float64,2)
        bpoints[1]=1
        bpoints[2]=n
        bregions[1]=1
        bregions[2]=2
        bnodefac[1]=1
        bnodefac[2]=1
        new(1,n,2,pts,eds,efac,nfac,bpoints,bregions,bnodefac)
    end
end

