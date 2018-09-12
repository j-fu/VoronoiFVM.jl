"""
    mutable struct FVMGraph

Fields:

    SpaceDim::Int32 # space dimension
    NumberOfNodes::Int64 # number of nodes
    NumberOfBoundaryRegions::Int64 # number of boundary regions
    Nodes::Array{Float64,2} # point coordinates.
    Edges::Array{Int32,2} # index array pointing into nodes list
    EdgeFactors::Array{Float64,1} # interface area / edge length
    NodeFactors::Array{Float64,1} # control volume sizes
    BoundaryNodes::Array{Int32,1} # indices of dirichlet nodes 
    BoundaryRegions::Array{Int32,1} # boundary region numbers


"""
mutable struct FVMGraph
    SpaceDim::Int32 # space dimension
    NumberOfNodes::Int64 # number of nodes
    NumberOfBoundaryRegions::Int64 # number of boundary regions
    Nodes::Array{Float64,2} # point coordinates.
    Edges::Array{Int32,2} # index array pointing into nodes list
    EdgeFactors::Array{Float64,1} # interface area / edge length
    NodeFactors::Array{Float64,1} # control volume sizes
    BoundaryNodes::Array{Int32,1} # indices of dirichlet nodes 
    BoundaryRegions::Array{Int32,1} # boundary region numbers
    BoundaryNodeFactors::Array{Float64,1} # boundary control volume sizes
    
    FVMGraph(X::Array{Float64,1})= FVMGraph(new(),X)
    FVMGraph(X::Array{Float64,1},Y::Array{Float64,1})= FVMGraph(new(),X,Y)
end
        
"""
   FVMGraph(X::Array{Float64,1})
   
Constructor for 1D finite volume graph
"""
function FVMGraph(self, X::Array{Float64,1})
    #  Primal Grid:
    #  o-----o-----o-----o-----o-----o-----o-----o-----o
    # Dual grid with control volumes
    #  |--|-----|-----|-----|-----|-----|-----|-----|--|
    n=length(X)
    nodes=reshape(X,1,n)
    edges=zeros(Int32,2,n-1)
    efac=zeros(n-1)
    nfac=zeros(n)
    for i=1:n-1 # loop over edges
        edges[1,i]=i
        edges[2,i]=i+1
        h=nodes[i+1]-nodes[i]
        efac[i]=1.0/h
        nfac[i]+=0.5*h
        nfac[i+1]+=0.5*h
    end
    bnodes=zeros(Int32,2)
    bregions=zeros(Int32,2)
    bnodefac=zeros(Float64,2)
    bnodes[1]=1
    bnodes[2]=n
    bregions[1]=1
    bregions[2]=2
    bnodefac[1]=1
    bnodefac[2]=1
    self.SpaceDim=1
    self.NumberOfNodes=n
    self.NumberOfBoundaryRegions=2
    self.Nodes=nodes
    self.Edges=edges
    self.EdgeFactors=efac
    self.NodeFactors=nfac
    self.BoundaryNodes=bnodes
    self.BoundaryRegions=bregions
    self.BoundaryNodeFactors=bnodefac
    return self
end


        

"""
   FVMGraph(X::Array{Float64,1},Y::Array{Float64,1})
   
Constructor for 2D finite volume graph
"""
function FVMGraph(self,X::Array{Float64,1},Y::Array{Float64,1})
    #  Primal Grid:
    #  o-----o-----o-----o-----o-----o-----o-----o-----o
    # Dual grid with control volumes
    #  |--|-----|-----|-----|-----|-----|-----|-----|--|
    nx=length(X)
    ny=length(Y)
    n=nx*ny
    nodes=Array{Float64,2}(undef,2,nx*ny)
    edges=Array{Int64,2}(undef,2,nx*(ny-1)+ny*(nx-1))
    efac=zeros(size(edges,2))
    nfac=zeros(size(nodes,2))
    edges.=0
    bnodes=Array{Int32,1}(undef,2*nx+2*ny)
    bregions=Array{Int32,1}(undef,2*nx+2*ny)
    bnodefac=Array{Float64,1}(undef,2*nx+2*ny)
    
    
    inode=1
    iedge=1
    ibnode=1
    for iy=1:ny
        for ix=1:nx
            dy=0
            dy=0
            nodes[1,inode]=X[ix]
            nodes[2,inode]=Y[iy]
            iedge0=iedge
            if iy==1 || iy==ny
                bnodes[ibnode]=inode
                fac=0.0
                if ix>1
                    fac+=0.5*X[ix]-X[ix-1]
                end
                if ix<nx
                    fac+=0.5*X[ix+1]-X[ix]
                end
                bnodefac[ibnode]=fac
                if iy==1
                    bregions[ibnode]=1
                end
                if iy==ny
                    bregions[ibnode]=3
                end
                ibnode+=1
            end
            if ix==1 || ix==nx
                bnodes[ibnode]=inode
                fac=0.0
                if iy>1
                    fac+=0.5*Y[iy]-Y[iy-1]
                end
                if iy<ny
                    fac+=0.5*Y[iy+1]-Y[iy]
                end
                bnodefac[ibnode]=fac
                if ix==1
                    bregions[ibnode]=4
                end
                if ix==nx
                    bregions[ibnode]=2
                end
                ibnode+=1
            end
            if ix<nx
                edges[1,iedge]=inode
                edges[2,iedge]=inode+1
                dx=X[ix+1]-X[ix]
                dy=0.0
                if iy>1
                    dy+=0.5*(Y[iy]-Y[iy-1])
                end
                if iy<ny
                    dy+=0.5*(Y[iy+1]-Y[iy])
                end
                efac[iedge]+=dy/dx
                iedge+=1
            end
            if iy<ny
                edges[1,iedge]=inode
                edges[2,iedge]=inode+nx
                dy=Y[iy+1]-Y[iy]
                dx=0.0
                if ix>1
                    dx+=0.5*(X[ix]-X[ix-1])
                end
                if ix<nx
                    dx+=0.5*(X[ix+1]-X[ix])
                end
                efac[iedge]+=dx/dy
                iedge+=1
            end
            if ix<nx && iy<ny
                dx=X[ix+1]-X[ix]
                dy=Y[iy+1]-Y[iy]
                nvol=0.25*dx*dy
                nfac[inode]+=nvol
                nfac[inode+1]+=nvol
                nfac[inode+nx]+=nvol
                nfac[inode+nx+1]+=nvol
            end
            inode+=1
        end
    end

    self.SpaceDim=1
    self.NumberOfNodes=n
    self.NumberOfBoundaryRegions=4
    self.Nodes=nodes
    self.Edges=edges
    self.EdgeFactors=efac
    self.NodeFactors=nfac
    self.BoundaryNodes=bnodes
    self.BoundaryRegions=bregions
    self.BoundaryNodeFactors=bnodefac
    return self
end
