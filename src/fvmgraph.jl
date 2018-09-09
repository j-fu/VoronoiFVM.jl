
# Structure for discrete geometry (for finite volume method)
# This can be easily generalized to higher dimensions 
# and unstructured grids
# 
# This is a graph embedded in R^d (d=1,2,3) with edge weights and node weights
# which can be generated from discretization grids in 1D/2D/3D
#
struct FVMGraph
    SpaceDim::Int32
    NumberOfNodes::Int64
    Points::Array{Float64,2} # point coordinates.
    Edges::Array{Int32,2} # index array pointing into points list
    EdgeFactors::Array{Float64,1} # interface area / edge length
    NodeFactors::Array{Float64,1} # control volume sizes
    BPoints::Array{Int32,1} # indices of dirichlet points 

    # Constructor for 1D discrete geometry
    function FVMGraph(X::Array{Float64,1})
    #  Primal Grid:
    #  o-----o-----o-----o-----o-----o-----o-----o-----o
    # Dual grid with control volumes
    #  |--|-----|-----|-----|-----|-----|-----|-----|--|
        n=length(X)
        pts=reshape(X,n,1)
        eds=zeros(Int32,n-1,2)
        efac=zeros(n-1)
        nfac=zeros(n)
        for i=1:n-1 # loop over edges
            eds[i,1]=i
            eds[i,2]=i+1
            h=pts[i+1]-pts[i]
            efac[i]=1.0/h
            nfac[i]+=0.5*h
            nfac[i+1]+=0.5*h
        end
        bpoints=zeros(Int32,2)
        bpoints[1]=1
        bpoints[2]=n
        new(1,n,pts,eds,efac,nfac,bpoints)
    end
end
