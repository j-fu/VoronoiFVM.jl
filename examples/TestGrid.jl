module TestGrid

using Printf
using TwoPointFluxFVM
using PyPlot

const Node=TwoPointFluxFVM.Node
const Edge=TwoPointFluxFVM.Edge

mutable struct Physics  <: TwoPointFluxFVM.Physics
    k::Float64
    eps::Float64 
    Physics()=new()
end


function breaction!(this::Physics,node::Node,f,bf,u,bu)
    if  this.bregion==2
        f[1]=this.k*(u[1]-bu[1])
        bf[1]=this.k*(bu[1]-u[1])+ this.k*(bu[1]-u[2])
        f[2]=this.k*(u[2]-bu[1])
    else
        f[1]=0        
        f[2]=0
    end
end

function bstorage!(this::Physics,node::Node,bf,bu)
    if  this.bregion==2
        bf[1]=bu[1]
    else
        bf[1]=0        
    end
end


function flux!(this::Physics,edge::Edge,f,uk,ul)
    f[1]=this.eps*(uk[1]-ul[1])
    f[2]=this.eps*(uk[2]-ul[2])
end 

function source!(this::Physics,node::Node,f)
    x1=node.coord[1]-0.5
    f[1]=exp(-20*x1^2)
end 



function Physics(this)
    TwoPointFluxFVM.PhysicsBase(this,2)
    this.num_bspecies=[ 0, 1, 0, 0]
    this.eps=1
    this.k=1.0
    this.breaction=breaction!
    this.bstorage=bstorage!
    this.flux=flux!
    this.source=source!
    return this
end




function main(;n=10,pyplot=true,verbose=false)
    h=1.0/convert(Float64,n-1)
    X=collect(0.0:h:1.0)
    Y=collect(0.0:h:1.0)
    grid1=TwoPointFluxFVM.Grid(X)
    grid2=TwoPointFluxFVM.Grid(X,Y)

    # print(summary(grid1))
    # println("")
    # show(grid1)
    # println("")
    # print(grid1)
    # println("")
    # dump(grid1)

    
    sys=TwoPointFluxFVM.SparseSystem(grid2,Physics(),2)
    TwoPointFluxFVM.add_species(sys,1,1)
    TwoPointFluxFVM.add_boundary_species(sys,2,3)
    # TwoPointFluxFVM.add_boundary_species(sys,3,2)

#    TwoPointFluxFVM.add_species(sys,2,1)


    println(length(sys.node_species.nzval))
    a=unknowns(sys)
    
    for i=1:nnodes(grid2)
        a[1,i]=i
    end

    for i=1:nbfaces(grid2)
        if grid2.bfaceregions[i]==3
            a[2,grid2.bfacenodes[1,i]]=111
            a[2,grid2.bfacenodes[2,i]]=111
        end
    end

    b=similar(a)
    println(typeof(b))
    bv=TwoPointFluxFVM
    
    println(a)
    println(a[1,:])
    println(a[2,:])

    println(typeof(a))
    println(typeof(a[1,:]))

    ax=a[1,:]
    println(typeof(ax))
    ax[5]=1221
    println(ax)
    println(a)
    
    ay=view(a,1,:)
    ay.=ax
    println(typeof(ay))
    ay[6]=2332
    println(ay)
    println(a)
    
    # println(typeof(a[1:2,:]))
    # println(a[1:2,:])

    if isdefined(sys, :matrix)
        println(sys.matrix)
    end
end




end
