module test_legacy
# Tests for keeping legacy API for versions <0.15

using VoronoiFVM
using Test
using Printf
using ExtendableGrids
using GridVisualize


function main101()

    ## Flux function which describes the flux
    ## between neighboring control volumes
    function g!(f,u,edge)
        f[1]=u[1,1]-u[1,2]
    end

    nspecies=1 ## Number of species
    ispec=1    ## Index of species we are working with

    ## Create a one dimensional discretization grid
    ## Each grid cell belongs to a region marked by a region number
    ## By default, there is only one region numbered with 1
    X=collect(0:0.2:1)
    grid=VoronoiFVM.Grid(X)

    ## Create a physics structure
    physics=VoronoiFVM.Physics(flux=g!)

    ## Create a finite volume system 
    sys=VoronoiFVM.System(grid,physics)

    ## Enable species 1 in region 1
    enable_species!(sys,ispec,[1])

    ## Set boundary conditions at boundary regions 1 and 2
    boundary_dirichlet!(sys,ispec,1,0.0)
    boundary_dirichlet!(sys,ispec,2,1.0)
    
    ## Create & initialize array for solution and initial value
    inival=unknowns(sys,inival=0)
    solution=unknowns(sys)

    ## Solve stationary problem
    solve!(solution,inival,sys)

    ## Return test value
    return solution[3]
end




function main102(;n=10,Plotter=nothing,verbose=false,D=0.01,v=1.0)

    function central_flux!(f,u,edge,data)
        f_diff=data.D*(u[1,1]-u[1,2])
        vh=project(edge,data.v)
        f[1]=f_diff+vh*(u[1,1]+u[1,2])/2
    end
    
    function upwind_flux!(f,u,edge,data)
        fdiff=data.D*(u[1,]-u[1,2])
        vh=project(edge,data.v)
        if vh>0
            f[1]=fdiff+vh*u[1,1]
        else
            f[1]=fdiff+vh*u[1,2]
        end
    end
    
    function bernoulli(x)
        if abs(x)<nextfloat(eps(typeof(x)))
            return 1
        end
        return x/(exp(x)-1)
    end
    
    function exponential_flux!(f,u,edge,data)
        vh=project(edge,data.v)
        Bplus= data.D*bernoulli(vh/data.D)
        Bminus=data.D*bernoulli(-vh/data.D)
        f[1]=Bminus*u[1,1]-Bplus*u[1,2]
    end
    

    function calculate(grid,data,flux,verbose)
        
        sys=VoronoiFVM.System(grid,VoronoiFVM.Physics(flux=flux, data=data))
        
        ## Add species 1 to region 1
        enable_species!(sys,1,[1])
        
        ## Set boundary conditions
        boundary_dirichlet!(sys,1,1,0.0)
        boundary_dirichlet!(sys,1,2,1.0)
        
        ## Create a solution array
        inival=unknowns(sys,inival=0.5)
        solution=unknowns(sys)
        
        ## Create solver control info
        control=VoronoiFVM.NewtonControl()
        control.verbose=verbose
        
        ## Stationary solution of the problem
        solve!(solution,inival,sys, control=control)
        return solution
    end
    
    ## Create a one-dimensional discretization
    h=1.0/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:1))
    
    data=(v=[v],D=D)
    
    # Calculate three stationary solutions with different ways to calculate flux
    solution_exponential=calculate(grid,data,exponential_flux!,verbose)
    solution_upwind=calculate(grid,data,upwind_flux!,verbose)
    solution_central=calculate(grid,data,central_flux!,verbose)
    
    # Visualize solutions using GridVisualize.jl
    p=GridVisualizer(Plotter=Plotter,layout=(3,1))
    scalarplot!(p[1,1],grid,solution_exponential[1,:],title="exponential")
    scalarplot!(p[2,1],grid,solution_upwind[1,:],title="upwind")
    scalarplot!(p[3,1],grid,solution_central[1,:],title="centered",show=true)
    
    ## Return test value
    return sum(solution_exponential)+sum(solution_upwind)+sum(solution_central)
end






function main103(;n=10,Plotter=nothing,D=0.01,v=1.0,tend=100)
    ## Bernoulli function used in the exponential fitting discretization
    function bernoulli(x)
        if abs(x)<nextfloat(eps(typeof(x)))
            return 1
        end
        return x/(exp(x)-1)
    end
    
    function exponential_flux!(f,u,edge,data)
        vh=project(edge,data.v)
        Bplus= data.D*bernoulli(vh/data.D)
        Bminus=data.D*bernoulli(-vh/data.D)
        f[1]=Bminus*u[1,1]-Bplus*u[1,2]
    end
    
    function outflow!(f,u,node,data)
        if node.region==2
            f[1]=data.v[1]*u[1]
        end
    end 
    
    
    ## Create a one-dimensional discretization
    h=1.0/n
    grid=VoronoiFVM.Grid(0:h:1)
    
    data=(v=[v],D=D)
    
    sys=VoronoiFVM.System(grid,VoronoiFVM.Physics(flux=exponential_flux!, data=data, breaction=outflow!))
    
    ## Add species 1 to region 1
    enable_species!(sys,1,[1])
    
    ## Set boundary conditions
    boundary_neumann!(sys,1,1,0.0)
    
    ## Create a solution array
    inival=unknowns(sys)
    inival[1,:].=map(x->1-2x,grid)
    
    ## Transient solution of the problem
    control=VoronoiFVM.NewtonControl()
    control.Δt=0.01*h
    control.Δt_min=0.01*h
    control.Δt_max=0.1*tend
    tsol=solve(inival,sys,[0,tend],control=control)
    
    vis=GridVisualizer(Plotter=Plotter)
    for i=1:length(tsol.t)
        scalarplot!(vis[1,1],grid,tsol[1,:,i],flimits=(0,1),title="t=$(tsol.t[i])",show=true)
        sleep(0.01)
    end
    maximum(tsol)<=1.0 && maximum(tsol[end])<1.0e-20
end


function main105(;n=10,Plotter=nothing,verbose=false, unknown_storage=:sparse)
    
    ## Create a one-dimensional discretization
    h=1.0/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:1))

    ## A parameter which is "passed" to the flux function via scope
    ϵ=1.0e-3

    ## Flux function which describes the flux
    ## between neighboring control volumes
    function flux!(f,u,edge)
        f[1]=ϵ*(u[1,1]-u[1,2])
    end

    ## Source term
    function source!(f,node)
        if node[1]<=0.5
            f[1]=1
        else
            f[1]=-1
        end
    end
    
    ## Reaction term
    function reaction!(f,u,node)
        f[1]=exp(u[1]) - exp(-u[1]) 
    end
    
    ## Create a physics structure
    physics=VoronoiFVM.Physics(
        flux=flux!,
        source=source!,
        reaction=reaction!)
    

    ## Create a finite volume system - either
    ## in the dense or  the sparse version.
    ## The difference is in the way the solution object
    ## is stored - as dense or as sparse matrix

    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)

    ## Add species 1 to region 1
    enable_species!(sys,1,[1])

    ## Set boundary conditions
    boundary_dirichlet!(sys,1,1,0.0)
    boundary_dirichlet!(sys,1,2,1.0)

    ## Create a solution array
    inival=unknowns(sys,inival=0.5)
    solution=unknowns(sys)

    ## Create solver control info
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose

    ## Stationary solution of the problem
    solve!(solution,inival,sys, control=control)

    scalarplot(grid,solution[1,:],title="Nonlinear Poisson", Plotter=Plotter)

    return sum(solution)
end



function test()
    @test main101() ≈ 0.4
    @test main102() ≈ 2.523569744561089
    @test main103()

    testval=1.5247901344230088
    @test main105(unknown_storage=:sparse) ≈ testval
    @test main105(unknown_storage=:dense) ≈ testval
    true
end

end
