# # 160: Unipolar degenerate drift-diffusion
# ([source code](SOURCE_URL))
module Example160_UnipolarDriftDiffusion1D

using Printf

using VoronoiFVM
using ExtendableGrids
using GridVisualize

mutable struct Data
    eps::Float64
    z::Float64
    ic::Int32
    iphi::Int32
    Data()=new()
end

function plot_solution(p,sys,U0,data)
    scalarplot!(p,sys.grid,U0[data.iphi,:], label="ψ", color=:green)
    scalarplot!(p,sys.grid,U0[data.ic,:], label="c-", color=:blue,show=true,clear=false)
end




function classflux!(f,u0,edge,data)
    u=unknowns(edge,u0)
    ic=data.ic
    iphi=data.iphi
    f[iphi]=data.eps*(u[iphi,1]-u[iphi,2])
    bp,bm=fbernoulli_pm(u[iphi,1]-u[iphi,2])
    f[ic]=bm*u[ic,1]-bp*u[ic,2]
end


function storage!(f,u,node,data)
    ic=data.ic
    iphi=data.iphi
    f[iphi]=0
    f[ic]=u[ic]
end

function reaction!(f,u,node,data)
    ic=data.ic
    iphi=data.iphi
    f[iphi]=data.z*(1-2*u[ic])
    f[ic]=0
end

function sedanflux!(f,u0,edge,data)
    u=unknowns(edge,u0)
    ic=data.ic
    iphi=data.iphi
    f[iphi]=data.eps*(u[iphi,1]-u[iphi,2])
    mu1=-log(1-u[ic,1])
    mu2=-log(1-u[ic,2])
    bp,bm=fbernoulli_pm(data.z*2*(u[iphi,1]-u[iphi,2])+(mu1-mu2))
    f[ic]=bm*u[ic,1]-bp*u[ic,2]
end


function main(;n=20,Plotter=nothing,dlcap=false,verbose=false,unknown_storage=:sparse)

    h=1.0/convert(Float64,n)
    grid=VoronoiFVM.Grid(collect(0:h:1))

    data=Data()
    data.eps=1.0e-4
    data.z=-1
    data.iphi=1
    data.ic=2

    ic=data.ic
    iphi=data.iphi


    physics=VoronoiFVM.Physics(data=data,
                               num_species=2,
                               flux=sedanflux!,
                               reaction=reaction!,
                               storage=storage!
                               )
    sys=VoronoiFVM.System(grid,physics,unknown_storage=unknown_storage)
    enable_species!(sys,1,[1])
    enable_species!(sys,2,[1])

    boundary_dirichlet!(sys,iphi,1,5.0)
    boundary_dirichlet!(sys,iphi,2,0.0)
    boundary_dirichlet!(sys,ic,2,0.5)


    inival=unknowns(sys)
    @views inival[iphi,:].=2
    @views inival[ic,:].=0.5
    U=unknowns(sys)

    p=GridVisualizer(Plotter=Plotter)
    plot_solution(p,sys,inival,data)


    data.eps=1.0e-3
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
            solve!(U,inival,sys,control=control,tstep=tstep)
            inival.=U
            u1=U[2]
            if verbose
                @printf("time=%g\n",t)
            end
            plot_solution(p,sys,U,data)
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
                solve!(U,inival,sys,control=control)
                inival.=U
                Q=integrate(sys,physics.reaction,U)
                sys.boundary_values[iphi,1]=dir*phi+delta
                solve!(U,inival,sys,control=control)
                inival.=U
                plot_solution(p,sys,U,data)
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

        px=GridVisualizer(Plotter=Plotter)
        scalarplot!(px,vplus,cdlplus,color=:green,clear=true)
        scalarplot!(px,vminus,cdlminus,color=:green,clear=false,show=true)
        return cdl
    end
end

function test()
        main(unknown_storage=:sparse) ≈ 0.9999546021312723 &&
            main(unknown_storage=:dense) ≈ 0.9999546021312723 &&
            main(dlcap=true) ≈ .010759276468375045 &&
            main(dlcap=true,unknown_storage=:dense) ≈ .010759276468375045
end
end
