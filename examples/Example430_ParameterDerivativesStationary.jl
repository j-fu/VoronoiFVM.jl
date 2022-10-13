#=

# 430: Parameter Derviatives (stationary)
 ([source code](SOURCE_URL))

This is still experimental.
=#
module Example430_ParameterDerivativesStationary


using VoronoiFVM,ExtendableGrids
using GridVisualize
using ExtendableSparse
using ForwardDiff,DiffResults


"""
    f(P)

Parameter dependent function which creates system and solves it
"""
function f(P)

    p=P[1]

    valuetype=typeof(p)
    n=10
    nspecies=1 
    ispec=1
    
    function flux!(f,u,edge)
        f[1]=u[1,1]^2-u[1,2]^2
    end
    
    function r!(f,u,edge)
        f[1]=p*u[1]^5
    end

    function bc!(f,u,node)
        boundary_dirichlet!(f,u,node,ispec,1,0.0)
        boundary_dirichlet!(f,u,node,ispec,3,p)
    end
    
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)
    sys=VoronoiFVM.System(grid;valuetype, species=[1], flux=flux!, reaction=r!, bcondition=bc!)
    tff=VoronoiFVM.TestFunctionFactory(sys)
    tfc=testfunction(tff,[1],[3])
    inival=unknowns(sys,)
    sol=solve(sys; inival=0.5)
    [integrate(sys,tfc,sol)[1]]
end

df(p)=ForwardDiff.derivative(f,p)


"""
    runf(;Plotter)

Plot f(p), df(p)
"""
function runf(;Plotter=nothing)
    P=0.1:0.05:2
    dresult=DiffResults.JacobianResult(ones(1))
    F=zeros(0)
    DF=zeros(0)
    @time for p ∈ P
        ForwardDiff.jacobian!(dresult,f,[p])
        push!(F,DiffResults.value(dresult)[1])
        push!(DF,DiffResults.jacobian(dresult)[1])
    end
    vis=GridVisualizer(;Plotter, legend=:lt)
    scalarplot!(vis,P,F,color=:red,label="f")
    scalarplot!(vis,P,DF,color=:blue,label="df",clear=false,show=true)
    sum(DF)
end





function fluxg!(f,u,edge,data)
    f[1]=u[1,1]^2-u[1,2]^2
end

function rg!(f,u,edge,data)
    f[1]=data.p*u[1]^5
end

function bcg!(f,u,node,data)
    boundary_dirichlet!(f,u,node,1,1,0.0)
    boundary_dirichlet!(f,u,node,1,3,data.p)
end

@Base.kwdef mutable struct MyData{Tv}
    p::Tv=1.0
end

"""
    rung(;Plotter)

Same as runf, but use DiffResults to calculate function at once. Also, pass parameter via data
"""
function rung(;Plotter=nothing, iteration=:cg, factorization=:default)
    n=10
    X=collect(0:1.0/n:1)
    grid=VoronoiFVM.Grid(X,X)

    # ugly but simple. By KISS we should first provide this way.

    sys=nothing
    data=nothing
    tfc=nothing

    function g(P)
        Tv=eltype(P)
        if isnothing(sys)
            data=MyData(one(Tv))
            sys=VoronoiFVM.System(grid;valuetype=Tv,species=[1], flux=fluxg!, reaction=rg!, bcondition=bcg!, data, unknown_storage=:dense)
            tff=VoronoiFVM.TestFunctionFactory(sys)
            tfc=testfunction(tff,[1],[3])
        end
        data.p=P[1]
        sol=solve(sys;inival=0.5,factorization,iteration)
        [integrate(sys,tfc,sol)[1]]
    end

    dresult=DiffResults.JacobianResult(ones(1))

    
    P=0.1:0.05:2
    G=zeros(0)
    DG=zeros(0)
    @time for p ∈ P
        ForwardDiff.jacobian!(dresult,g,[p])
        push!(G,DiffResults.value(dresult)[1])
        push!(DG,DiffResults.jacobian(dresult)[1])
    end

    vis=GridVisualizer(;Plotter, legend=:lt)
    scalarplot!(vis,P,G,color=:red,label="g")
    scalarplot!(vis,P,DG,color=:blue,label="dg",clear=false,show=true)
    sum(DG)
end




function test()
    testval=270.21652747842154
    xf=runf()
    xg=rung()
    xg2=rung(factorization=:ilu0, iteration=:cg)
    xf≈ testval &&
        xg≈ testval &&
         xg2≈ testval
    
end

end

