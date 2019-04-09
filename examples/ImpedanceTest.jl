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
mutable struct Physics  <: VoronoiFVM.Physics
    flux::Function      # flux function, mandatory
    reaction::Function  # reaction term, optional
    storage::Function   # storage term, mandatory
    D::Real           # Example for "user data" passed to the callback
    C::Real
    R::Real
    Physics()=new() # Provide inner constructor resulting in uninitialized struct
end



# Main function for user interaction from REPL and
# for testimg. Default physics need to generate correct
# test value.


function main(;n=11,pyplot=false,verbose=false, dense=false)
    L=1

    
    h=1.0/convert(Float64,n-1)
    grid=VoronoiFVM.Grid(collect(0:h:L))
    
    
    physics=Physics()
    physics.R=1
    physics.D=1
    physics.C=2

    physics.flux=function(physics,edge,f,uk,ul)
        f[1]=physics.D*(uk[1]-ul[1])
    end 


    physics.storage=function(physics,node, f,u)
        f[1]=physics.C*u[1]
    end

    physics.reaction=function(physics,node,f,u)
        f[1]=physics.R*u[1]

    end

    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics,1)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics,1)
    end
    add_species(sys,1,[1])

    factory=VoronoiFVM.TestFunctionFactory(sys)
    tf0=testfunction(factory,[2],[1])
    tfL=testfunction(factory,[1],[2])


    
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0
    
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet
    
    inival=unknowns(sys)
    inival.=0.0
    
    U=solve(sys,inival)

    ω0=0.5
    ω1=1.0e4

    ω=ω0
    
    isys=VoronoiFVM.ImpedanceSystem(sys,U,1,1)

    allomega=zeros(0)
    allI0=zeros(Complex{Float64},0)
    allIL=zeros(Complex{Float64},0)
    allIx0=zeros(Complex{Float64},0)
    allIxL=zeros(Complex{Float64},0)

    testval=0.0
    while ω<ω1
        
        iω=1im*ω
        
        z=sqrt(iω*physics.C/physics.D+physics.R/physics.D);
        eplus=exp(z*L);
        eminus=exp(-z*L);
        Ix0=-physics.D*z*(eminus+eplus)/(eminus-eplus); 
        IxL=2.0*physics.D*z/(eminus-eplus);
        
        
        UZ=solve(isys,ω)

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

