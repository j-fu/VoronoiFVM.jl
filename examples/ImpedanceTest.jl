module ImpedanceTest

# This gives us he @printf macro (c-like output)
using Printf

# That's the thing we want to do
using VoronoiFVM

# Allow plotting
if isinteractive()
    using PyPlot
end



#
# Structure containing  userdata information
#
# We choose a mutable struct which allows to overwrite
# fields later.
mutable struct Data  <: VoronoiFVM.Data
    D::Real           # Example for "user data" passed to the callback
    C::Real
    R::Real
    Data()=new() # Provide inner constructor resulting in uninitialized struct
end



# Main function for user interaction from REPL and
# for testimg. Default physics need to generate correct
# test value.


function main(;n=11,pyplot=false,verbose=false, dense=false)
    L=1

    
    h=1.0/convert(Float64,n-1)
    grid=VoronoiFVM.Grid(collect(0:h:L))
    
    data=Data()
    data.R=1
    data.D=1
    data.C=2

    physics=VoronoiFVM.Physics(
        data=data,
        
    flux=function(f,u,edge,data)
        f[1]=data.D*(u[1]-u[2])
    end,


    storage=function(f,u,node,data)
        f[1]=data.C*u[1]
    end,

    reaction=function(f,u,node,data)
        f[1]=data.R*u[1]
    end)

    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end
    enable_species!(sys,1,[1])

    factory=VoronoiFVM.TestFunctionFactory(sys)
    tf0=testfunction(factory,[2],[1])
    tfL=testfunction(factory,[1],[2])


    
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0
    
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet
    
    inival=unknowns(sys)
    U=unknowns(sys)
    inival.=0.0
    
    solve!(U,inival,sys)

    ω0=0.5
    ω1=1.0e4

    ω=ω0

    excited_spec=1
    excited_bc=1
    isys=VoronoiFVM.ImpedanceSystem(sys,U,excited_spec, excited_bc)

    allomega=zeros(0)
    allI0=zeros(Complex{Float64},0)
    allIL=zeros(Complex{Float64},0)
    allIx0=zeros(Complex{Float64},0)
    allIxL=zeros(Complex{Float64},0)

    testval=0.0
    UZ=unknowns(isys)
    while ω<ω1
        
        iω=1im*ω
        
        z=sqrt(iω*data.C/data.D+data.R/data.D);
        eplus=exp(z*L);
        eminus=exp(-z*L);
        Ix0=-data.D*z*(eminus+eplus)/(eminus-eplus); 
        IxL=2.0*data.D*z/(eminus-eplus);
        
        
        solve!(UZ,isys,ω)

        if pyplot
            PyPlot.clf()
            PyPlot.grid()
            plot(grid.coord[1,:],real(UZ[1,:]),label="Re")
            plot(grid.coord[1,:],imag(UZ[1,:]),label="Im")
            pause(1.0e-10)
        end    
        
        I0=integrate(isys,tf0,ω,UZ)[1]
        IL=integrate(isys,tfL,ω,UZ)[1]
        push!(allomega, ω)
        push!(allI0,I0)
        push!(allIL,IL)
        push!(allIx0,Ix0)
        push!(allIxL,IxL)
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
    return  imag(allIL[5])
end
end

