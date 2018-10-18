"""
    mutable struct FVMGraph

Fields:

    dim_space::Int32 # space dimension
    num_nodes::Int64 # number of nodes
    num_bregions::Int64 # number of boundary regions
    node_coordinates::Array{Float64,2} # point coordinates
    edge_nodes::Array{Int64,2} # index array pointing into nodes list
    edge_factors::Array{Float64,1} # interface area / edge length
    node_factors::Array{Float64,1} # control volume sizes
    bnode_nodes::Array{Int64,1} # indices of dirichlet nodes 
    bnode_regions::Array{Int32,1} # boundary region numbers
    bnode_factors::Array{Float64,1} # boundary control volume sizes
    num_bregion_nodes::Array{Int64,1} # number of nodes per boundary region
    bnode_index::Array{Int64,1} # index of boundary node in boundary region 


"""
mutable struct FVMGraph
    dim_space::Int32 # space dimension
    num_nodes::Int64 # number of nodes
    num_bregions::Int64 # number of boundary regions
    node_coordinates::Array{Float64,2} # point coordinates
    edge_nodes::Array{Int64,2} # index array pointing into nodes list
    edge_factors::Array{Float64,1} # interface area / edge length
    node_factors::Array{Float64,1} # control volume sizes
    bnode_nodes::Array{Int64,1} # indices of dirichlet nodes 
    bnode_regions::Array{Int32,1} # boundary region numbers
    bnode_factors::Array{Float64,1} # boundary control volume sizes
    num_bregion_nodes::Array{Int64,1} # number of nodes per boundary region
    bnode_index::Array{Int64,1} # index of boundary node in boundary region 
    
    FVMGraph(X::Array{Float64,1})= FVMGraph(new(),X)
    FVMGraph(X::Array{Float64,1},Y::Array{Float64,1})= FVMGraph(new(),X,Y)
end
        
"""
   FVMGraph(X::Array{Float64,1})
   
Constructor for 1D finite volume graph
"""
function FVMGraph(this::FVMGraph, X::Array{Float64,1})
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
    bnode_nodes=zeros(Int32,2)
    bnode_regions=zeros(Int32,2)
    bnodefac=zeros(Float64,2)
    bnode_nodes[1]=1
    bnode_nodes[2]=n
    bnode_regions[1]=1
    bnode_regions[2]=2
    bnodefac[1]=1
    bnodefac[2]=1
    this.dim_space=1
    this.num_nodes=n
    this.num_bregions=2
    this.node_coordinates=nodes
    this.edge_nodes=edges
    this.edge_factors=efac
    this.node_factors=nfac
    this.bnode_nodes=bnode_nodes
    this.bnode_regions=bnode_regions
    this.bnode_factors=bnodefac

    finalize(this)
    return this
end


        

"""
   FVMGraph(X::Array{Float64,1},Y::Array{Float64,1})
   
Constructor for 2D finite volume graph
"""
function FVMGraph(this::FVMGraph,X::Array{Float64,1},Y::Array{Float64,1})
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
    bnode_nodes=Array{Int32,1}(undef,2*nx+2*ny)
    bnode_regions=Array{Int32,1}(undef,2*nx+2*ny)
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
                bnode_nodes[ibnode]=inode
                fac=0.0
                if ix>1
                    fac+=0.5*(X[ix]-X[ix-1])
                end
                if ix<nx
                    fac+=0.5*(X[ix+1]-X[ix])
                end
                bnodefac[ibnode]=fac
                if iy==1
                    bnode_regions[ibnode]=1
                end
                if iy==ny
                    bnode_regions[ibnode]=3
                end
                ibnode+=1
            end
            if ix==1 || ix==nx
                bnode_nodes[ibnode]=inode
                fac=0.0
                if iy>1
                    fac+=0.5*(Y[iy]-Y[iy-1])
                end
                if iy<ny
                    fac+=0.5*(Y[iy+1]-Y[iy])
                end
                bnodefac[ibnode]=fac
                if ix==1
                    bnode_regions[ibnode]=4
                end
                if ix==nx
                    bnode_regions[ibnode]=2
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

    this.dim_space=1
    this.num_nodes=n
    this.num_bregions=4
    this.node_coordinates=nodes
    this.edge_nodes=edges
    this.edge_factors=efac
    this.node_factors=nfac
    this.bnode_nodes=bnode_nodes
    this.bnode_regions=bnode_regions
    this.bnode_factors=bnodefac

    finalize(this)
    return this
end

"""
   Finalize construction of graph
"""
function finalize(this::FVMGraph)
    this.num_bregion_nodes=zeros(Int64,this.num_bregions)
    this.bnode_index=zeros(Int64,length(this.bnode_nodes))
    
    for i=1:length(this.bnode_nodes)
        ireg=this.bnode_regions[i]
        this.num_bregion_nodes[ireg]+=1
        this.bnode_index[i]=this.num_bregion_nodes[ireg]
    end
end

