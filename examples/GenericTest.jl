module GenericTest

# TODO:
# More realistic test - use boundary flux as functional
#

using VoronoiFVM,ExtendableGrids
using GridVisualize
using ExtendableSparse
using ForwardDiff,DiffResults
using Plots

function g!(f,u,edge)
    f[1]=u[1,1]^2-u[1,2]^2
end

function f(X)

    x=X[1]

    valuetype=typeof(x)
    n=50
    nspecies=1 
    ispec=1
    
    function r!(f,u,edge)
        f[1]=x*u[1]^5
    end

    function bc!(f,u,node)
        boundary_dirichlet!(f,u,node,ispec,1,0.0)
        boundary_dirichlet!(f,u,node,ispec,3,x)
    end
    
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)
    sys=VoronoiFVM.System(grid;valuetype, species=[1], flux=g!, reaction=r!, bcondition=bc!)
    sol=solve(sys;inival=0.5)
    [sum(sol)]
end

df(x)=ForwardDiff.derivative(f,x)


function runf()
    X=0.1:0.05:2
    dresult=DiffResults.JacobianResult(ones(1))

    
    F=zeros(0)
    DF=zeros(0)
    @time for x ∈ X
        ForwardDiff.jacobian!(dresult,f,[x])
        push!(F,DiffResults.value(dresult)[1])
        push!(DF,DiffResults.jacobian(dresult)[1])
    end
    p=plot()
    plot!(p, X,F)
    plot!(p, X,DF)
    p
end




@Base.kwdef mutable struct MyData{Tv}
    x::Tv=1.0
end



function Base.merge(vdata::MyData{Tv}, params::AbstractVector{Tv}) where Tv
    vdata.x=params[1]
    vdata
end

function Base.merge(data::MyData{Tu}, params::AbstractVector{Tv}) where {Tu, Tv}
    vdata=copy!(MyData{Tv}(),data)
    vdata.x=params[1]
    vdata
end




function gx!(f,u,edge,data)
    f[1]=u[1,1]^2-u[1,2]^2
end

function rx!(f,u,edge,data)
    f[1]=data.x*u[1]^5
end

function bcx!(f,u,node,data)
    boundary_dirichlet!(f,u,node,1,1,0.0)
    boundary_dirichlet!(f,u,node,1,3,data.x)
end




function runfx()
    n=50
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)

    # ugly but simple. By KISS we should first provide this way.
    sys=nothing
    data=nothing
    function fx(X)
        Tv=eltype(X)
        if isnothing(sys)
            data=MyData{Tv}(1.0)
            sys=VoronoiFVM.System(grid;valuetype=Tv,species=[1], flux=gx!, reaction=rx!, bcondition=bcx!, data, unknown_storage=:dense)
        end
        data.x=X[1]
        sol=solve(sys,inival=0.5)
        [sum(sol)]
    end

    dresult=DiffResults.JacobianResult(ones(1))

    
    X=0.1:0.05:2
    F=zeros(0)
    DF=zeros(0)
    @time for x ∈ X
        ForwardDiff.jacobian!(dresult,fx,[x])
        push!(F,DiffResults.value(dresult)[1])
        push!(DF,DiffResults.jacobian(dresult)[1])
    end

    p=plot()
    plot!(p, X,F)
    plot!(p, X,DF)
    p
end



function gp!(f,u,edge)
    f[1]=u[1,1]^2-u[1,2]^2
end

function rp!(f,u,edge)
    p=parameters(u)
    f[1]=p[1]*u[1]^5
end

function bcp!(f,u,node)
    p=parameters(u)
    boundary_dirichlet!(f,u,node,1,1,0.0)
    boundary_dirichlet!(f,u,node,1,3,p[1])
end

#=
With parameters(u,node), we probably can cover both cases.
But first we can tag a new version with the current variant:

- rename dudp dfdp
- Parameter derivative would be sufficient to be calculated with solution, not
  during iterative processs.
- But this would mean have another assembly loop: 
   eval_dudp(system, sol, u), and we would pass parameters extra with given solution and be consistent
   with the global approach. Check with fvsys, probably we did the same there

  This is the way it _should be done_. How to pass this through ? We just need to evaluate with parameters
  as duals. Easy split of eval_and_assemble ?
  node.params 
-

=#
function runfp()
    n=50
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)
    sys=VoronoiFVM.System(grid;species=[1], flux=gp!, reaction=rp!, bcondition=bcp!, nparams=1, unknown_storage=:dense)

    # f(u,p)=0
    # was ist du/dp ? Newton ist statinär
    # also

    # u  =     df/du^-1f(u,p)
    # df/du u = f(u,p)
    # df/du u_p = df/dp(u,p)
    
    function fp(X)
        sol=solve(sys,inival=0.5,params=X)
        sum(sol), -sum(sys.matrix\vec(sys.dudp[1]))
    end
    
    X=0.1:0.05:2
    F=zeros(0)
    DF=zeros(0)
    @time for x ∈ X
        f,df=fp([x])
        push!(F,f)
        push!(DF,df)
    end

    p=plot()
    plot!(p, X,F)
    plot!(p, X,DF)
    p
end



function main(;Plotter=nothing,n=5, Tv=Float64)
    nspecies=1 
    ispec=1    
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)
    physics=VoronoiFVM.Physics(flux=g!)
    sys=VoronoiFVM.System(grid,physics;Tv)
    enable_species!(sys,ispec,[1])
    boundary_dirichlet!(sys,ispec,1,0.0)
    boundary_dirichlet!(sys,ispec,3,1.0)
    inival=unknowns(sys,inival=0.5)
    solve(sys;inival,factorization=SparspakLU(valuetype=Tv))
    # vis=GridVisualizer(Plotter=Plotter)
    # scalarplot!(vis,grid,solution[1,:],clear=true,colormap=:summer)
    # reveal(vis)
end

## Called by unit test

function test()
    main() ≈ 0.2
end

end

