module ImpedanceTest
#
#  Impedance calculation for 
#
#    C u_t - (D u_x)_x + Ru = 0   in (0,1)
#      u(0,t)=1 + exp(iωt)
#      u(1,t)=0
#
#    Measurement: I(t)= D u_x(1,t)
#     
#    Steady state: 
#    - (D u0_x)_x + Ru0 = 0
#    u0(0,t)=1 
#    u0(1,t)=0 
#
#    Small signal ansatz for ω
#
#    u(x,t)= u0(x)+ ua(x) exp(iωt)
#
#    iωC ua - (D ua_x)_x + R u_a =0
#      ua(0)=1
#      ua(1)=0
# 
#        

using Printf
using VoronoiFVM

if isinteractive()
    using PyPlot
end

# Structure containing  userdata information
mutable struct Data  <: VoronoiFVM.AbstractData
    D::Real           
    C::Real
    R::Real
    Data()=new()
end


function main(;nref=0,pyplot=false,verbose=false)
    L=1.0

    # Create array which is refined close to 0
    h0=0.1/2.0^nref
    h1=0.5/2.0^nref
    X=VoronoiFVM.geomspace(0.0,L,h0,h1)

    # Create discretitzation grid
    grid=VoronoiFVM.Grid(X)

    # Create and fill data 
    data=Data()
    data.R=1
    data.D=1
    data.C=2

    # Declare constitutive functions
    flux=function(f,u,edge,data)
        f[1]=data.D*(u[1]-u[2])
    end

    storage=function(f,u,node,data)
        f[1]=data.C*u[1]
    end

    reaction=function(f,u,node,data)
        f[1]=data.R*u[1]
    end

    # Create physics struct
    physics=VoronoiFVM.Physics(data=data,
                               flux=flux,
                               storage=storage,
                               reaction=reaction
                               )
    # Create discrete system and enabe species
    sys=VoronoiFVM.DenseSystem(grid,physics)
    enable_species!(sys,1,[1])

    # Create test functions for current measurement
    factory=VoronoiFVM.TestFunctionFactory(sys)
    tf0=testfunction(factory,[2],[1])
    tfL=testfunction(factory,[1],[2])

    # Solve steady state problem
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0
    
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet

    inival=unknowns(sys)
    steadystate=unknowns(sys)
    inival.=0.0
    solve!(steadystate,inival,sys)

    # Create Impeadancs system from steady state
    excited_spec=1
    excited_bc=1
    isys=VoronoiFVM.ImpedanceSystem(sys,steadystate,excited_spec, excited_bc)

    # Prepare recording of impedance results
    allomega=zeros(0)

    # for calculated data
    allI0=zeros(Complex{Float64},0)
    allIL=zeros(Complex{Float64},0)

    # for exact data
    allIx0=zeros(Complex{Float64},0)
    allIxL=zeros(Complex{Float64},0)

    ω0=0.5
    ω1=1.0e4
    ω=ω0

    testval=0.0
    UZ=unknowns(isys)
    while ω<ω1
        
        iω=1im*ω

        
        # solve impedance system
        solve!(UZ,isys,ω)

        # calculate aproximate solution
        I0=integrate(isys,tf0,ω,UZ)[1]
        IL=integrate(isys,tfL,ω,UZ)[1]

        # record approximate solution
        push!(allomega, ω)
        push!(allI0,I0)
        push!(allIL,IL)

        # record exact solution
        z=sqrt(iω*data.C/data.D+data.R/data.D);
        eplus=exp(z*L);
        eminus=exp(-z*L);
        Ix0=-data.D*z*(eminus+eplus)/(eminus-eplus); 
        IxL=2.0*data.D*z/(eminus-eplus);
        push!(allIx0,Ix0)
        push!(allIxL,IxL)

        if pyplot
            PyPlot.clf()
            PyPlot.grid()
            plot(grid.coord[1,:],real(UZ[1,:]),label="Re")
            plot(grid.coord[1,:],imag(UZ[1,:]),label="Im")
            pause(1.0e-10)
        end    

        # increase omega
        ω=ω*1.2

    end
    
    if pyplot
        PyPlot.clf()
        PyPlot.grid()
        plot(real(allI0),imag(allI0),label="calc")
        plot(real(allIx0),imag(allIx0),label="exact")
        PyPlot.legend(loc="upper left")
        pause(1.0e-10)
        waitforbuttonpress()

        PyPlot.clf()
        PyPlot.grid()
        plot(real(allIL),imag(allIL),label="calc")
        plot(real(allIxL),imag(allIxL),label="exact")
        PyPlot.legend(loc="upper left")
        pause(1.0e-10)
        waitforbuttonpress()
    end
    #return test value
    return  imag(allIL[5])
end
end

