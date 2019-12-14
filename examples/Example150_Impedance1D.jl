# # 150: Impedance calculation
# ([source code](SOURCE_URL))
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

module Example150_Impedance1D

using Printf
using VoronoiFVM

# Structure containing  userdata information
mutable struct Data  <: VoronoiFVM.AbstractData
    D::Float64           
    C::Float64
    R::Float64
    Data()=new()
end


function main(;nref=0,Plotter=nothing,verbose=false, dense=false)

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
    if dense
        sys=VoronoiFVM.DenseSystem(grid,physics)
    else
        sys=VoronoiFVM.SparseSystem(grid,physics)
    end
    enable_species!(sys,1,[1])

    # Create test functions for current measurement

    excited_bc=1
    excited_bcval=1.0
    excited_spec=1


    factory=VoronoiFVM.TestFunctionFactory(sys)
    measurement_testfunction=testfunction(factory,[1],[2])

    # Solve steady state problem
    sys.boundary_values[1,1]=1.0
    sys.boundary_values[1,2]=0
    
    sys.boundary_factors[1,1]=VoronoiFVM.Dirichlet
    sys.boundary_factors[1,2]=VoronoiFVM.Dirichlet

    inival=unknowns(sys)
    steadystate=unknowns(sys)
    inival.=0.0
    solve!(steadystate,inival,sys)

    function meas_stdy(meas,U)
        u=reshape(U,sys)
        meas[1]=VoronoiFVM.integrate_stdy(sys,measurement_testfunction,u)[1]
        nothing
    end

    function meas_tran(meas,U)
        u=reshape(U,sys)
        meas[1]=VoronoiFVM.integrate_tran(sys,measurement_testfunction,u)[1]
        nothing
    end
    

    dmeas_stdy=measurement_derivative(sys,meas_stdy,steadystate)
    dmeas_tran=measurement_derivative(sys,meas_tran,steadystate)



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
        # obtain measurement in frequency  domain
        IL=freqdomain_impedance(isys,ω,steadystate,excited_spec,excited_bc,excited_bcval,dmeas_stdy, dmeas_tran)

        # record approximate solution
        push!(allomega, ω)
        push!(allIL,IL)

        # record exact solution
        z=sqrt(iω*data.C/data.D+data.R/data.D);
        eplus=exp(z*L);
        eminus=exp(-z*L);
        IxL=2.0*data.D*z/(eminus-eplus);
        push!(allIxL,IxL)

        # increase omega
        ω=ω*1.2

    end
    
    if isplots(Plotter)
        p=Plotter.plot(grid=true)
        Plotter.plot!(p,real(allIL),imag(allIL),label="calc")
        Plotter.plot!(p,real(allIxL),imag(allIxL),label="exact")
        Plotter.gui(p)
    end
    #return test value
    return  imag(allIL[5])
end

function test()
    main(dense=true) ≈ 0.23106605162049176 &&  main(dense=false) ≈ 0.23106605162049176
end


end

