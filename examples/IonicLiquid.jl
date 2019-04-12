module IonicLiquid

using Printf
if isinteractive()
    using PyPlot
end

using VoronoiFVM

mutable struct Physics <: VoronoiFVM.Physics
    flux::Function
    storage::Function
    reaction::Function
    eps::Float64 
    z::Float64
    ic::Int32
    iphi::Int32
    Physics()=new()
end



function plot_solution(sys,U0)
    physics=sys.physics
    iphi=physics.iphi
    ic=physics.ic
    PyPlot.clf()
    @views begin
        PyPlot.plot(sys.grid.coord[1,:],U0[iphi,:], label="Potential", color="g")
        PyPlot.plot(sys.grid.coord[1,:],U0[ic,:], label="c-", color="b")
    end
    PyPlot.grid()
    PyPlot.legend(loc="upper right")
    PyPlot.pause(1.0e-10)
end


function main(;n=20,pyplot=false,dlcap=false,verbose=false,dense=false)
    
    h=1.0/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:1))
    
    physics=Physics()
    physics.eps=1.0e-4
    physics.z=-1
    physics.iphi=1
    physics.ic=2

    ic=physics.ic
    iphi=physics.iphi

    classflux=function(physics,edge,f,u)
        nspecies=2
        uk=VoronoiFVM.UK(u,nspecies)
        ul=VoronoiFVM.UL(u,nspecies)
        ic=physics.ic
        iphi=physics.iphi
        f[iphi]=physics.eps*(uk[iphi]-ul[iphi])
        arg=uk[iphi]-ul[iphi]
        bp,bm=fbernoulli_pm(uk[iphi]-ul[iphi])
        f[ic]=bm*uk[ic]-bp*ul[ic]
    end 


    physics.flux=function(physics,edge,f,u)
        nspecies=2
        uk=VoronoiFVM.UK(u,nspecies)
        ul=VoronoiFVM.UL(u,nspecies)
        ic=physics.ic
        iphi=physics.iphi
        f[iphi]=physics.eps*(uk[iphi]-ul[iphi])
        muk=-log(1-uk[ic])
        mul=-log(1-ul[ic])
        bp,bm=fbernoulli_pm(physics.z*2*(uk[iphi]-ul[iphi])+(muk-mul))
        f[ic]=bm*uk[ic]-bp*ul[ic]
    end 


    physics.storage=function(physics,node, f,u)
        ic=physics.ic
        iphi=physics.iphi
        f[iphi]=0
        f[ic]=u[ic]
    end

    physics.reaction=function(physics,node, f,u)
        ic=physics.ic
        iphi=physics.iphi
        f[iphi]=physics.z*(1-2*u[ic])
        f[ic]=0
    end

    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics,2)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics,2)
    end

    add_species(sys,1,[1])
    add_species(sys,2,[1])

    
    sys.boundary_values[iphi,1]=5
    sys.boundary_values[iphi,2]=0.0
    
    sys.boundary_factors[iphi,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[iphi,2]=VoronoiFVM.Dirichlet

    sys.boundary_values[ic,2]=0.5
    sys.boundary_factors[ic,2]=VoronoiFVM.Dirichlet
    
    inival=unknowns(sys)
    @views inival[iphi,:].=2
    @views inival[ic,:].=0.5
    U=unknowns(sys)


    if pyplot
        plot_solution(sys,inival)
    end
    
    physics.eps=1.0e-3
    control=VoronoiFVM.NewtonControl()
    control.verbose=verbose
    u1=0
    if !dlcap
        control.damp_initial=0.5
        t=0.0
        tend=1.0
        tstep=1.0e-4
        while t<tend
            t=t+tstep
            solve(sys,inival,U,control=control,tstep=tstep)
            inival.=U
            u1=U[2]
            if verbose
                @printf("time=%g\n",t)
            end
            if pyplot
                plot_solution(sys,U)
            end
            tstep*=1.4
        end
        return u1
    else
        delta=1.0e-4
        @views inival[iphi,:].=0
        @views inival[ic,:].=0.5
        sys.boundary_values[iphi,1]=0
        
        dphi=1.0e-1
        phimax=5
        delta=1.0e-4
        vplus=zeros(0)
        cdlplus=zeros(0)
        vminus=zeros(0)
        cdlminus=zeros(0)
        cdl=0
        for dir in [1,-1]
            phi=0.0
            while phi<phimax
                sys.boundary_values[iphi,1]=dir*phi
                solve(sys,inival,U,control=control)
                inival.=U
                Q=integrate(sys,physics.reaction,U)
                sys.boundary_values[iphi,1]=dir*phi+delta
                solve(sys,inival,U,control=control)
                inival.=U
                if pyplot
                    plot_solution(sys,U)
                end
                Qdelta=integrate(sys,physics.reaction,U)
                cdl=(Qdelta[iphi]-Q[iphi])/delta
                if dir==1
                    push!(vplus,dir*phi)
                    push!(cdlplus,cdl)
                else
                    push!(vminus,dir*phi)
                    push!(cdlminus,cdl)
                end
                phi+=dphi
            end
        end
        if pyplot
            PyPlot.clf()
            PyPlot.plot(vplus,cdlplus,color="g")
            PyPlot.plot(vminus,cdlminus,color="g")
            PyPlot.grid()
            PyPlot.legend(loc="upper right")
            PyPlot.pause(1.0e-10)
        end
        return cdl
    end
end

end
